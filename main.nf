#!/usr/bin/env nextflow
/*
 * nf-demultiplex-matfet
 * Maternal/fetal demultiplexing of single-cell samples, in one of two modes.
 *
 * --mode pregnancy (default): pregnancy samples where mother and fetus both contribute
 * many cells, one pregnancy group at a time:
 *   SOUPORCELL_DEMUX (per library: remap + allele counts)
 *     -> GROUP_DEMUX (per group: group matrix, Souporcell discovery, composition gate)
 *     -> LABEL_MATERNAL_FETAL (per group: marker orientation, per-library split, report)
 * See docs/pipeline_design.md.
 *
 * --mode fetal: fetal tissue, one fetus per group, maternal cells rare and no
 * maternal sample:
 *   FETAL (Souporcell k=2 on the first library of each batch + cross-library heatmap check;
 *          per library: cellsnp-lite pileup; per batch: fetal genotype + cell scores)
 * See docs/fetal_mode_design.md.
 */

nextflow.enable.dsl = 2

include { SOUPORCELL_DEMUX     } from './subworkflows/local/souporcell_demux'
include { GROUP_DEMUX          } from './subworkflows/local/group_demux'
include { LABEL_MATERNAL_FETAL } from './subworkflows/local/label_maternal_fetal'
include { FETAL                } from './subworkflows/local/fetal'

def helpMessage() {
    log.info """
    nf-demultiplex-matfet — maternal/fetal demultiplexing

    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results -profile sanger [--mode fetal]

    Modes:
      pregnancy      (default) Pregnancy samples where mother and fetus both contribute
                     many cells: Souporcell group discovery, composition gate, marker
                     labelling. docs/pipeline_design.md
      fetal          Fetal tissue, maternal cells rare and no maternal sample: pooled fetal
                     genotype per batch of libraries, every cell scored against it. Valid
                     only when maternal cells are a small share (<= ~2%) of each batch.
                     Souporcell k=2 runs first on the checked libraries and their clusters are
                     compared across libraries (heatmap + warnings in the report; it never
                     changes the calls). docs/fetal_mode_design.md

    Samplesheet (CSV, header required; relative paths resolve against the
    samplesheet's own directory):
      group_id,sample_id,bam,barcodes,expression[,sample_role][,expected_genotypes][,souporcell_dir][,cellsnp_dir]

      group_id            pregnancy: libraries sharing one mother and one fetus.
                          fetal: libraries from one fetus.
                          A group of one library is allowed.
      sample_id           Unique library identifier across the whole samplesheet.
      bam                 CellRanger or STARsolo BAM with CB/UB tags (indexed, or index is created).
      barcodes            Filtered cell-barcode file (optionally .gz).
      expression          Filtered expression matrix directory (CellRanger
                          filtered_feature_bc_matrix or STARsolo GeneFull/filtered).
      sample_role         pregnancy only. Optional, default 'mixture'. 'maternal_reference'
                          is not supported yet.
      expected_genotypes  pregnancy only. Optional, default 2. Only 2 is supported.
      souporcell_dir      Optional. An existing per-library Souporcell output directory
                          (alt.mtx, ref.mtx, barcodes.tsv, common_variants_covered.vcf).
                          pregnancy: the library skips remapping and allele counting.
                          fetal: the library joins the Souporcell check without a new
                          run, whether or not it would be selected (needs alt.mtx, ref.mtx,
                          clusters.tsv, common_variants_covered.vcf).
      cellsnp_dir         fetal mode only. Optional. An existing cellsnp-lite output directory
                          for the library (cellSNP.base.vcf.gz, cellSNP.tag.DP.mtx,
                          cellSNP.tag.AD.mtx, cellSNP.samples.tsv) made with the same panel and
                          flags; the library then skips the pileup.

    Fetal mode batches: each group's libraries, in samplesheet order, are scored in batches
    of --fetal_batch_size (default 6, as validated).
    Fetal mode Souporcell check (--fetal_souporcell): 'per_batch' (default; the first library
    of each batch, as in the paper's range checks), 'all', or 'off'.
    """.stripIndent()
}

/*
 * Pre-flight: this pipeline runs every tool out of params.container
 * (souporcell_pipeline.py, samtools, numpy/scipy). With no container engine the task
 * scripts run on the bare host and fail with exit 127 hours into Souporcell, so refuse
 * to start instead.
 */
def checkContainerEngine() {
    if (!workflow.containerEngine) {
        error """
        No container engine is enabled, but every process in this pipeline runs its
        tools from a container. Re-run with Singularity enabled, e.g.:
          nextflow run main.nf --input <samplesheet.csv> -profile sanger
        """.stripIndent()
    }
    // Catch a mistyped or missing local image before any task is submitted.
    if (params.container ==~ /^\/.*\.(sif|img)$/ && !file(params.container).exists()) {
        error "Container image not found: ${params.container}"
    }
}

/*
 * Read and validate the whole samplesheet before any job is submitted.
 * Returns a list of [meta, barcodes, bam, expression, counts_dir or null]
 * in samplesheet order, with meta = [id, group, index, group_size]. counts_dir is the
 * library's existing allele counts: souporcell_dir in pregnancy mode, cellsnp_dir in
 * fetal mode (cellsnp_dir is refused in pregnancy mode). In fetal mode,
 * meta.soup_dir holds the library's existing Souporcell output for the check, or null.
 */
def resolvePath(sheet_dir, String p) {
    return p.startsWith('/') ? file(p) : sheet_dir.resolve(p)
}

def mustExist(f, what, row_id) {
    if (!f.exists()) { error "Samplesheet row '${row_id}': ${what} not found: ${f}" }
    return f
}

def parseSamplesheet(sheet_path, mode) {
    def sheet = file(sheet_path, checkIfExists: true)
    def sheet_dir = sheet.parent

    def raw_rows = sheet.splitCsv(header: true, strip: true)
    if (!raw_rows) { error "Samplesheet ${sheet} has no rows" }

    def rows = []
    raw_rows.eachWithIndex { row, index ->
        [ 'group_id', 'sample_id', 'bam', 'barcodes', 'expression' ].each { col ->
            if (!row[col]?.trim()) {
                error "Samplesheet row ${index + 1} is missing required column '${col}': ${row}"
            }
        }
        def id    = row.sample_id.trim()
        def group = row.group_id.trim()
        [ sample_id: id, group_id: group ].each { col, value ->
            if (!(value ==~ /^[A-Za-z0-9._-]+$/)) {
                error "Samplesheet ${col} '${value}' may only contain letters, digits, '.', '_' and '-'"
            }
        }

        // Columns that belong to the other mode are refused rather than ignored.
        def wrong_mode = mode == 'fetal' ? [ 'sample_role', 'expected_genotypes' ] : [ 'cellsnp_dir' ]
        wrong_mode.each { col ->
            if (row[col]?.trim()) {
                error "Samplesheet row '${id}': column '${col}' does not apply to --mode ${mode}"
            }
        }

        def role = row.sample_role?.trim() ?: 'mixture'
        if (role == 'maternal_reference') {
            error "Samplesheet row '${id}': sample_role 'maternal_reference' is not supported yet; " +
                  "maternal reference libraries have not been validated (docs/pipeline_design.md)"
        }
        if (role != 'mixture') {
            error "Samplesheet row '${id}': unknown sample_role '${role}' (expected 'mixture')"
        }

        def k = row.expected_genotypes?.trim() ?: '2'
        if (k != '2') {
            error "Samplesheet row '${id}': expected_genotypes must be 2 (one mother, one fetus); " +
                  "got '${k}'. Other values have not been validated."
        }

        def counts_dir = null
        def soup_dir = null
        if (row.souporcell_dir?.trim()) {
            soup_dir = mustExist(resolvePath(sheet_dir, row.souporcell_dir.trim()), 'souporcell_dir', id)
            // fetal needs only what the cross-library check reads (bin/souporcell_distance.py)
            def needed = mode == 'fetal'
                ? [ 'alt.mtx', 'ref.mtx', 'clusters.tsv', 'common_variants_covered.vcf' ]
                : [ 'alt.mtx', 'ref.mtx', 'barcodes.tsv', 'common_variants_covered.vcf' ]
            needed.each { name -> mustExist(soup_dir.resolve(name), "souporcell_dir file ${name}", id) }
            if (mode == 'pregnancy') { counts_dir = soup_dir }
        }
        if (row.cellsnp_dir?.trim()) {
            counts_dir = mustExist(resolvePath(sheet_dir, row.cellsnp_dir.trim()), 'cellsnp_dir', id)
            [ 'cellSNP.base.vcf.gz', 'cellSNP.tag.DP.mtx', 'cellSNP.tag.AD.mtx', 'cellSNP.samples.tsv' ].each { name ->
                mustExist(counts_dir.resolve(name), "cellsnp_dir file ${name}", id)
            }
        }

        def meta = [ id: id, group: group, index: index ]
        if (mode == 'fetal') { meta.soup_dir = soup_dir }
        rows << [
            meta,
            mustExist(resolvePath(sheet_dir, row.barcodes.trim()),   'barcodes',   id),
            mustExist(resolvePath(sheet_dir, row.bam.trim()),        'bam',        id),
            mustExist(resolvePath(sheet_dir, row.expression.trim()), 'expression', id),
            counts_dir,
        ]
    }

    def dups = rows.collect { it[0].id }.countBy { it }.findAll { _id, n -> n > 1 }.keySet()
    if (dups) { error "Duplicate sample_id(s) in samplesheet: ${dups.join(', ')}" }

    def sizes = rows.collect { it[0].group }.countBy { it }
    rows.each { it[0].group_size = sizes[it[0].group] }
    return rows
}

/*
 * Fetal mode: cut each group's libraries, in samplesheet order, into batches of
 * batch_size. Adds meta.batch (<group_id>_b<k>), meta.batch_size, meta.batch_pos (1 for
 * the batch's first library) and meta.group_batches (batches in the group) to every row.
 */
def assignBatches(rows, int batch_size) {
    rows.groupBy { it[0].group }.each { group, grows ->
        def batches = grows.sort { it[0].index }.collate(batch_size)
        batches.eachWithIndex { batch, k ->
            batch.eachWithIndex { row, pos ->
                row[0].batch         = "${group}_b${k + 1}".toString()
                row[0].batch_pos     = pos + 1
                row[0].batch_size    = batch.size()
                row[0].group_batches = batches.size()
            }
        }
    }
    return rows
}

/*
 * Fetal mode: the scoring settings were validated as a whole (D8 benchmark, embryo runs).
 * Changing one is allowed but warned about.
 */
def warnUnvalidatedFetalSettings() {
    def validated = [ fetal_batch_size: 6, fetal_min_sites: 1000, fetal_rhos: '0,0.1,0.25' ]
    validated.each { name, value ->
        if (params[name].toString() != value.toString()) {
            log.warn "--${name} ${params[name]} differs from the validated value ${value}; " +
                     "fetal mode's error rates were measured at ${value} (docs/fetal_mode_design.md)"
        }
    }
}

def completionMessage() {
    if (workflow.success && params.mode == 'fetal') {
        def check = params.fetal_souporcell == 'off' ? '' : 'fetal/souporcell_check/, '
        log.info "Done. Results in ${params.outdir} (fetal/<group_id>/, ${check}<sample_id>/fetal_mode_cells.tsv.gz, reports/<group_id>.html)."
    }
    else if (workflow.success) {
        log.info "Done. Results in ${params.outdir} (groups/<group_id>/, <sample_id>/, reports/<group_id>.html)."
    }
    else {
        log.error "Pipeline finished with errors (see above). Any partial results are in ${params.outdir}."
    }
}

workflow {
    if (params.help || !params.input) {
        helpMessage()
        if (!params.input && !params.help) {
            log.error "Missing required --input samplesheet. See --help."
            System.exit(1)
        }
        return
    }

    if (!(params.mode in [ 'pregnancy', 'fetal' ])) {
        error "Unknown --mode '${params.mode}' (expected 'pregnancy' or 'fetal'). See --help."
    }

    // Site files have no defaults outside -profile sanger (conf/sanger.config).
    if (!params.container) {
        error "--container is required: a Singularity image with Souporcell, cellsnp-lite, samtools and Python (-profile sanger sets the farm image)"
    }
    if (!params.soc_vcf) {
        error "--soc_vcf is required: the common-variant panel VCF (-profile sanger sets the farm panel)"
    }
    if (!params.soc_fasta && !(params.mode == 'fetal' && params.fetal_souporcell == 'off')) {
        error "--soc_fasta is required for Souporcell's remap: the genome FASTA the BAMs were aligned to (-profile sanger sets the farm FASTA)"
    }
    checkContainerEngine()

    // Souporcell always runs at k=2 (modules/local/souporcell.nf); other values were never validated.
    if (params.containsKey('souporcell_k') && params.souporcell_k.toString() != '2') {
        error "--souporcell_k is fixed at 2 in this pipeline (got '${params.souporcell_k}')"
    }

    workflow.onComplete = {
        completionMessage()
    }

    def rows = parseSamplesheet(params.input, params.mode)

    // ---- Fetal mode: pileup per library, then pooled fetal genotype and cell scores per batch ----
    if (params.mode == 'fetal') {
        if (!(params.fetal_batch_size instanceof Integer) || params.fetal_batch_size < 1) {
            error "--fetal_batch_size must be a positive integer, got '${params.fetal_batch_size}'"
        }
        if (!(params.fetal_souporcell in [ 'per_batch', 'all', 'off' ])) {
            error "--fetal_souporcell must be 'per_batch', 'all' or 'off', got '${params.fetal_souporcell}'"
        }
        warnUnvalidatedFetalSettings()
        assignBatches(rows, params.fetal_batch_size)
        def n_batches = rows.collect { it[0].batch }.unique().size()
        log.info "Fetal mode: ${rows.size()} libraries in ${n_batches} batches of up to ${params.fetal_batch_size}"

        // Souporcell check: the first library of each batch (or all), plus any library with a
        // souporcell_dir; none when switched off.
        def soc_rows = params.fetal_souporcell == 'off' ? [] : rows.findAll { meta, _b, _bam, _e, _c ->
            params.fetal_souporcell == 'all' || meta.batch_pos == 1 || meta.soup_dir != null
        }
        log.info "Fetal mode Souporcell check: ${soc_rows.size()} libraries (--fetal_souporcell ${params.fetal_souporcell}, " +
                 "${soc_rows.count { it[0].soup_dir == null }} to run)"
        ch_souporcell = Channel.fromList(soc_rows.collect { meta, barcodes, bam, _e, _c -> [ meta, barcodes, bam, meta.soup_dir ] })

        FETAL(Channel.fromList(rows), ch_souporcell, file(params.soc_vcf, checkIfExists: true))
        return
    }

    // One group meta per pregnancy group; libraries listed in samplesheet order.
    def gmetas = [:]
    rows.each { meta, _barcodes, _bam, _expression, _soup_dir ->
        gmetas.get(meta.group, [ id: meta.group, samples: [], k: 2 ]).samples << meta.id
    }

    def meta_by_id = rows.collectEntries { meta, _barcodes, _bam, _expression, _soup_dir -> [ meta.id, meta ] }

    ch_rows = Channel.fromList(rows)

    // ---- Per library: remap + allele counts, unless an existing Souporcell dir is given ----
    // Only the sample id goes in, so the cached remap depends on the library alone:
    // adding, removing or regrouping other samplesheet rows must not redo hours of
    // remapping on -resume. The group fields are reattached by id afterwards.
    ch_to_count = ch_rows
        .filter { _meta, _barcodes, _bam, _expression, soup_dir -> soup_dir == null }
        .map { meta, barcodes, bam, _expression, _soup_dir -> [ [ id: meta.id ], barcodes, bam ] }
    SOUPORCELL_DEMUX(ch_to_count)

    ch_counts = SOUPORCELL_DEMUX.out.dir.map { lib, dir -> [ meta_by_id[lib.id], dir ] }.mix(
        ch_rows
            .filter { _meta, _barcodes, _bam, _expression, soup_dir -> soup_dir != null }
            .map { meta, _barcodes, _bam, _expression, soup_dir -> [ meta, soup_dir ] }
    )

    // Each group moves on as soon as its own libraries are counted; within a group,
    // libraries are put back in samplesheet order so the group matrix is reproducible.
    ch_groups = ch_counts
        .map { meta, dir -> [ groupKey(meta.group, meta.group_size), meta.index, dir ] }
        .groupTuple()
        .map { key, indices, dirs ->
            def ordered = [ indices, dirs ].transpose().sort { it[0] }.collect { it[1] }
            [ gmetas[key.toString()], ordered ]
        }

    // ---- Per group: group matrix, Souporcell discovery, composition gate ----
    GROUP_DEMUX(ch_groups, file(params.soc_vcf, checkIfExists: true))

    // ---- Per group: maternal/fetal orientation, per-library split, report ----
    ch_expression = Channel.fromList(
        gmetas.values().collect { gmeta ->
            [ gmeta, rows.findAll { it[0].group == gmeta.id }.collect { it[3] } ]
        }
    )

    LABEL_MATERNAL_FETAL(
        GROUP_DEMUX.out.discovery,
        GROUP_DEMUX.out.provenance,
        GROUP_DEMUX.out.gate,
        ch_expression,
        file(params.marker_panel, checkIfExists: true)
    )
}
