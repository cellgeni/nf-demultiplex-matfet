#!/usr/bin/env nextflow
/*
 * nf-demultiplex-matfet
 * Maternal/fetal demultiplexing of single-cell pregnancy samples, one pregnancy group
 * at a time:
 *   SOUPORCELL_DEMUX (per library: remap + allele counts)
 *     -> GROUP_DEMUX (per group: group matrix, Souporcell discovery, composition gate)
 *     -> LABEL_MATERNAL_FETAL (per group: marker orientation, per-library split, report)
 * See docs/pipeline_design.md.
 */

nextflow.enable.dsl = 2

include { SOUPORCELL_DEMUX     } from './subworkflows/local/souporcell_demux'
include { GROUP_DEMUX          } from './subworkflows/local/group_demux'
include { LABEL_MATERNAL_FETAL } from './subworkflows/local/label_maternal_fetal'

def helpMessage() {
    log.info """
    nf-demultiplex-matfet — group-aware maternal/fetal demultiplexing

    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results -profile sanger

    Samplesheet (CSV, header required; relative paths resolve against the
    samplesheet's own directory):
      group_id,sample_id,bam,barcodes,expression[,sample_role][,expected_genotypes][,souporcell_dir]

      group_id            Pregnancy group: libraries sharing one mother and one fetus.
                          A group of one library is allowed.
      sample_id           Unique library identifier across the whole samplesheet.
      bam                 CellRanger BAM with cell barcodes (indexed, or index is created).
      barcodes            Filtered cell-barcode file (optionally .gz).
      expression          CellRanger filtered_feature_bc_matrix directory.
      sample_role         Optional, default 'mixture'. 'maternal_reference' is not
                          supported yet.
      expected_genotypes  Optional, default 2. Only 2 is supported.
      souporcell_dir      Optional. An existing per-library Souporcell output directory
                          (alt.mtx, ref.mtx, barcodes.tsv, common_variants_covered.vcf);
                          the library then skips remapping and allele counting.

    See docs/pipeline_design.md for the design, the composition gate and the outputs.
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
 * Returns a list of [meta, barcodes, bam, expression, souporcell_dir or null]
 * in samplesheet order, with meta = [id, group, index, group_size].
 */
def resolvePath(sheet_dir, String p) {
    return p.startsWith('/') ? file(p) : sheet_dir.resolve(p)
}

def mustExist(f, what, row_id) {
    if (!f.exists()) { error "Samplesheet row '${row_id}': ${what} not found: ${f}" }
    return f
}

def parseSamplesheet(sheet_path) {
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

        def soup_dir = null
        if (row.souporcell_dir?.trim()) {
            soup_dir = mustExist(resolvePath(sheet_dir, row.souporcell_dir.trim()), 'souporcell_dir', id)
            [ 'alt.mtx', 'ref.mtx', 'barcodes.tsv', 'common_variants_covered.vcf' ].each { name ->
                mustExist(soup_dir.resolve(name), "souporcell_dir file ${name}", id)
            }
        }

        rows << [
            [ id: id, group: group, index: index ],
            mustExist(resolvePath(sheet_dir, row.barcodes.trim()),   'barcodes',   id),
            mustExist(resolvePath(sheet_dir, row.bam.trim()),        'bam',        id),
            mustExist(resolvePath(sheet_dir, row.expression.trim()), 'expression', id),
            soup_dir,
        ]
    }

    def dups = rows.collect { it[0].id }.countBy { it }.findAll { _id, n -> n > 1 }.keySet()
    if (dups) { error "Duplicate sample_id(s) in samplesheet: ${dups.join(', ')}" }

    def sizes = rows.collect { it[0].group }.countBy { it }
    rows.each { it[0].group_size = sizes[it[0].group] }
    return rows
}

def completionMessage() {
    if (workflow.success) {
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

    checkContainerEngine()

    def rows = parseSamplesheet(params.input)

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

    workflow.onComplete = {
        completionMessage()
    }
}
