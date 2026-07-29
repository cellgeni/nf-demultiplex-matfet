// Aggregate panel-marker expression per Souporcell cluster (singlets only).

process MARKER_SUMMARY {
    tag "${meta.id}"
    // Publish only the summary; the souporcell dir flows through the channel for
    // downstream processes but is already published by SOUPORCELL (avoid duplicating it).
    publishDir "${params.outdir}/${meta.id}", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(meta), path(souporcell_dir), path(expression)
    path marker_panel

    output:
    tuple val(meta), path("cluster_marker_summary.tsv"), path("${souporcell_dir}"), emit: summary
    path "gene_resolution.tsv",                                                     emit: resolution

    script:
    """
    aggregate_markers.py \\
        --sample-id ${meta.id} \\
        --clusters ${souporcell_dir}/clusters.tsv \\
        --matrix ${expression} \\
        --panel ${marker_panel} \\
        --out cluster_marker_summary.tsv \\
        --resolution-out gene_resolution.tsv
    """
}
