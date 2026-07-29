// Human-readable per-sample HTML report (verdict + evidence + demux QC + provenance).

process BUILD_REPORT {
    tag "${meta.id}"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(meta), path(marker_summary), path(cluster_labels), path(qc_summary), path(souporcell_dir)

    output:
    tuple val(meta), path("${meta.id}.html"), emit: report

    script:
    """
    make_report.py \\
        --sample-id ${meta.id} \\
        --marker-summary ${marker_summary} \\
        --cluster-labels ${cluster_labels} \\
        --qc-summary ${qc_summary} \\
        --souporcell-qc ${souporcell_dir}/qc.txt \\
        --out ${meta.id}.html
    """
}
