// Human-readable per-group HTML report (verdict + gate + per-library counts + evidence).

process BUILD_REPORT {
    tag "${gmeta.id}"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(gmeta), path(marker_summary), path(cluster_labels), path(qc_summary), path(gate), path(gate_libraries)

    output:
    tuple val(gmeta), path("${gmeta.id}.html"), emit: report

    script:
    """
    make_report.py \\
        --group-id ${gmeta.id} \\
        --marker-summary ${marker_summary} \\
        --cluster-labels ${cluster_labels} \\
        --qc-summary ${qc_summary} \\
        --gate ${gate} \\
        --gate-libraries ${gate_libraries} \\
        --out ${gmeta.id}.html
    """
}
