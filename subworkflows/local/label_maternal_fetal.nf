// LABEL_MATERNAL_FETAL — orient a group's genotype clusters as maternal/fetal/unresolved
// from pooled marker expression, then split the labelled cells back into libraries.
//   MARKER_SUMMARY -> LABEL_CLUSTERS -> SPLIT_BY_LIBRARY, BUILD_REPORT
// Clusters are only oriented when the composition gate verdict is two_profiles.

include { MARKER_SUMMARY   } from '../../modules/local/marker_summary'
include { LABEL_CLUSTERS   } from '../../modules/local/label_clusters'
include { SPLIT_BY_LIBRARY } from '../../modules/local/split_by_library'
include { BUILD_REPORT     } from '../../modules/local/report'

workflow LABEL_MATERNAL_FETAL {
    take:
    ch_discovery   // tuple(gmeta, clusters, ambient, profiles)
    ch_provenance  // tuple(gmeta, barcode_provenance.tsv)
    ch_gate        // tuple(gmeta, gate.tsv, gate_libraries.tsv)
    ch_expression  // tuple(gmeta, [expression_dir, ...]) in samplesheet order
    marker_panel   // path

    main:
    ch_marker_in = ch_discovery
        .join(ch_provenance)
        .join(ch_expression)
        .map { gmeta, clusters, _ambient, _profiles, provenance, expression ->
            [ gmeta, clusters, provenance, expression ] }
    MARKER_SUMMARY(ch_marker_in, marker_panel)

    ch_label_in = MARKER_SUMMARY.out.summary
        .join(ch_discovery)
        .join(ch_gate)
        .map { gmeta, summary, clusters, ambient, _profiles, gate, _gate_libraries ->
            [ gmeta, summary, clusters, ambient, gate ] }
    LABEL_CLUSTERS(ch_label_in)

    SPLIT_BY_LIBRARY(
        LABEL_CLUSTERS.out.assignments
            .join(ch_provenance)
            .join(ch_gate)
            .map { gmeta, assignments, provenance, gate, _gate_libraries ->
                [ gmeta, assignments, provenance, gate ] }
    )

    BUILD_REPORT(
        MARKER_SUMMARY.out.summary
            .join(LABEL_CLUSTERS.out.labels)
            .join(ch_gate)
    )

    emit:
    cluster_labels = LABEL_CLUSTERS.out.labels
    libraries      = SPLIT_BY_LIBRARY.out.libraries
    report         = BUILD_REPORT.out.report
}
