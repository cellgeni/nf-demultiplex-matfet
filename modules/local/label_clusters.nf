// Apply the two-tier maternal/fetal decision engine to a group's clusters, only when
// the composition gate accepted two profiles, and propagate labels to barcodes.

process LABEL_CLUSTERS {
    tag "${gmeta.id}"
    publishDir "${params.outdir}/groups/${gmeta.id}", mode: 'copy',
        pattern: '{cluster_labels.tsv,group_qc_summary.tsv}'

    input:
    tuple val(gmeta), path(marker_summary), path(clusters), path(ambient), path(gate)

    output:
    tuple val(gmeta), path("cluster_labels.tsv"), path("group_qc_summary.tsv"), emit: labels
    tuple val(gmeta), path("cell_assignments.tsv.gz"),                          emit: assignments

    script:
    """
    label_clusters.py \\
        --group-id ${gmeta.id} \\
        --marker-summary ${marker_summary} \\
        --clusters ${clusters} \\
        --ambient ${ambient} \\
        --gate ${gate} \\
        --cluster-labels-out cluster_labels.tsv \\
        --cell-assignments-out cell_assignments.tsv.gz \\
        --qc-out group_qc_summary.tsv \\
        --min-cluster-cells ${params.min_cluster_cells} \\
        --min-marker-detected-fraction ${params.min_marker_detected_fraction} \\
        --min-marker-fraction-delta ${params.min_marker_fraction_delta} \\
        --min-marker-logfc ${params.min_marker_logfc} \\
        --rps4y1-fraction-floor ${params.rps4y1_fraction_floor} \\
        --rps4y1-min-fold ${params.rps4y1_min_fold} \\
        --min-vote-margin ${params.min_vote_margin} \\
        --imbalance-warn-ratio ${params.imbalance_warn_ratio}
    """
}
