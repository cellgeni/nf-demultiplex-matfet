// Apply the two-tier maternal/fetal decision engine and propagate labels to barcodes.

process LABEL_CLUSTERS {
    tag "${meta.id}"
    publishDir "${params.outdir}/${meta.id}", mode: 'copy'

    input:
    tuple val(meta), path(marker_summary), path(souporcell_dir)

    output:
    tuple val(meta), path("cluster_labels.tsv"), path("sample_qc_summary.tsv"), emit: labels
    tuple val(meta), path("cell_assignments.tsv.gz"),                           emit: assignments

    script:
    """
    label_clusters.py \\
        --sample-id ${meta.id} \\
        --marker-summary ${marker_summary} \\
        --clusters ${souporcell_dir}/clusters.tsv \\
        --ambient ${souporcell_dir}/ambient_rna.txt \\
        --cluster-labels-out cluster_labels.tsv \\
        --cell-assignments-out cell_assignments.tsv.gz \\
        --qc-out sample_qc_summary.tsv \\
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
