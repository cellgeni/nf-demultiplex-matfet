#!/usr/bin/env nextflow
/*
 * nf-demux-preg
 * Maternal/fetal demultiplexing of single-cell pregnancy samples.
 *   SOUPORCELL_DEMUX (K=2)  ->  LABEL_MATERNAL_FETAL (marker panel)
 * See ROADMAP.md.
 */

nextflow.enable.dsl = 2

include { SOUPORCELL_DEMUX     } from './subworkflows/local/souporcell_demux'
include { LABEL_MATERNAL_FETAL } from './subworkflows/local/label_maternal_fetal'

def helpMessage() {
    log.info """
    nf-demux-preg — maternal/fetal demultiplexing (Souporcell K=2 + marker labeling)

    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results -profile sanger,singularity

    Samplesheet (CSV, header required):
      sample_id,bam,barcodes,expression

      sample_id    Unique library identifier (no spaces).
      bam          CellRanger BAM with cell barcodes (indexed, or index is created).
      barcodes     Filtered cell-barcode file (optionally .gz).
      expression   CellRanger filtered_feature_bc_matrix directory.

    K is fixed at 2 (one maternal + one fetal genotype); it is not a samplesheet column.
    See ROADMAP.md §7 for the input contract and §8 for the decision logic.
    """.stripIndent()
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

    // ---- Parse and light-validate the samplesheet ----
    ch_rows = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            [ 'sample_id', 'bam', 'barcodes', 'expression' ].each { col ->
                if (!row.containsKey(col) || !row[col]?.trim()) {
                    error "Samplesheet row is missing required column '${col}': ${row}"
                }
            }
            def meta       = [ id: row.sample_id.trim() ]
            def bam        = file(row.bam.trim(),        checkIfExists: true)
            def barcodes   = file(row.barcodes.trim(),   checkIfExists: true)
            def expression = file(row.expression.trim(), checkIfExists: true)
            [ meta, barcodes, bam, expression ]
        }

    // Fail fast on duplicate sample_ids.
    ch_rows
        .map { meta, _barcodes, _bam, _expression -> meta.id }
        .toList()
        .map { ids ->
            def dups = ids.countBy { it }.findAll { _k, v -> v > 1 }.keySet()
            if (dups) { error "Duplicate sample_id(s) in samplesheet: ${dups.join(', ')}" }
        }

    // ---- Subworkflow 1: Souporcell demux (K=2) ----
    ch_soup_in = ch_rows.map { meta, barcodes, bam, _expression -> [ meta, barcodes, bam ] }
    SOUPORCELL_DEMUX(ch_soup_in)

    // ---- Subworkflow 2: maternal/fetal labeling ----
    ch_expr = ch_rows.map { meta, _barcodes, _bam, expression -> [ meta, expression ] }
    ch_label_in = SOUPORCELL_DEMUX.out.dir
        .join(ch_expr)                                    // (meta, soup_dir, expression)
        .map { meta, soup_dir, expression -> [ meta, soup_dir, expression ] }

    LABEL_MATERNAL_FETAL(ch_label_in, file(params.marker_panel, checkIfExists: true))

    workflow.onComplete = {
        log.info "Done. Results in ${params.outdir} (per-sample tables + reports/<sample_id>.html)."
    }
}

/*
 * Standalone entry: run ONLY the labeling subworkflow against existing Souporcell output.
 * This is the "LABEL_MATERNAL_FETAL is independently runnable/testable" design rule (ROADMAP §4).
 *
 *   nextflow run main.nf -entry LABEL_ONLY --input label_samplesheet.csv -profile test
 *
 * Samplesheet columns: sample_id,souporcell_dir,expression
 *   souporcell_dir must contain clusters.tsv and ambient_rna.txt (and qc.txt for the report).
 */
workflow LABEL_ONLY {
    if (params.help || !params.input) {
        helpMessage()
        if (!params.input) { log.error "Missing required --input samplesheet."; System.exit(1) }
        return
    }

    ch_label_in = Channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .map { row ->
            [ 'sample_id', 'souporcell_dir', 'expression' ].each { col ->
                if (!row.containsKey(col) || !row[col]?.trim()) {
                    error "LABEL_ONLY samplesheet row missing required column '${col}': ${row}"
                }
            }
            // Resolve relative paths against projectDir so the shipped fixture runs
            // from any launch directory.
            def resolve = { p -> p.startsWith('/') ? file(p, checkIfExists: true)
                                                    : file("${projectDir}/${p}", checkIfExists: true) }
            def meta       = [ id: row.sample_id.trim() ]
            def soup_dir   = resolve(row.souporcell_dir.trim())
            def expression = resolve(row.expression.trim())
            [ meta, soup_dir, expression ]
        }

    LABEL_MATERNAL_FETAL(ch_label_in, file(params.marker_panel, checkIfExists: true))

    workflow.onComplete = {
        log.info "Done. Results in ${params.outdir} (per-sample tables + reports/<sample_id>.html)."
    }
}
