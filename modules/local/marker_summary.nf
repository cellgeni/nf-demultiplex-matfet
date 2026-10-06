// Aggregate panel-marker expression per group genotype cluster (singlets only),
// pooled across every library in the group.

process MARKER_SUMMARY {
    tag "${gmeta.id}"
    publishDir "${params.outdir}/groups/${gmeta.id}", mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(gmeta), path(clusters), path(provenance), path(expression, stageAs: 'expr*/*')
    path marker_panel

    output:
    tuple val(gmeta), path("cluster_marker_summary.tsv"), emit: summary
    path "gene_resolution.tsv",                           emit: resolution

    script:
    // Expression directories arrive in samplesheet order, matching gmeta.samples.
    def rows = [gmeta.samples, expression instanceof List ? expression : [expression]].transpose()
        .collect { sample_id, dir -> "'${sample_id}\t${dir}'" }
        .join(' ')
    """
    printf '%s\\n' ${rows} > matrices.tsv

    aggregate_markers.py \\
        --group-id ${gmeta.id} \\
        --clusters ${clusters} \\
        --provenance ${provenance} \\
        --matrices matrices.tsv \\
        --panel ${marker_panel} \\
        --out cluster_marker_summary.tsv \\
        --resolution-out gene_resolution.tsv
    """
}
