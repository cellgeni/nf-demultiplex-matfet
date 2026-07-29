// LABEL_MATERNAL_FETAL — label the two genotype clusters as maternal/fetal/unresolved.
//
// Stable contract (for standalone runs): its input interface is the *Souporcell output
// files* (clusters.tsv, ambient_rna.txt) + the expression matrix — NOT nf-demultiplex
// internals. Anyone with existing Souporcell output can run just this subworkflow.

include { MARKER_SUMMARY } from '../../modules/local/marker_summary'
include { LABEL_CLUSTERS } from '../../modules/local/label_clusters'
include { BUILD_REPORT   } from '../../modules/local/report'

workflow LABEL_MATERNAL_FETAL {
    take:
    ch_label_in    // tuple(meta, souporcell_dir, expression)
    marker_panel   // path

    main:
    MARKER_SUMMARY(ch_label_in, marker_panel)

    LABEL_CLUSTERS(
        MARKER_SUMMARY.out.summary.map { meta, summary, soup_dir -> [meta, summary, soup_dir] }
    )

    // Assemble report inputs by joining on meta.
    ch_report = MARKER_SUMMARY.out.summary
        .map { meta, summary, soup_dir -> [meta, summary, soup_dir] }
        .join(LABEL_CLUSTERS.out.labels)          // + cluster_labels, qc_summary
        .map { meta, summary, soup_dir, labels, qc -> [meta, summary, labels, qc, soup_dir] }

    BUILD_REPORT(ch_report)

    emit:
    cluster_labels = LABEL_CLUSTERS.out.labels
    assignments    = LABEL_CLUSTERS.out.assignments
    report         = BUILD_REPORT.out.report
}
