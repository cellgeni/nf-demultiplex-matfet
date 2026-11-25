// Fetal mode, per fetus: human-readable HTML report (bin/make_fetal_report.py), with the run's
// Souporcell cross-library check when it ran (souporcell_check/, empty input when switched off).

process FETAL_REPORT {
    tag "${gid}"
    publishDir "${params.outdir}/reports", mode: 'copy'

    input:
    tuple val(gid), path(summary), path(sex_check), path(warnings), path(flagged), path(soc_check)

    output:
    tuple val(gid), path("${gid}.html"), emit: report

    script:
    def check = soc_check ? "--souporcell-check ${soc_check}" : ''
    """
    make_fetal_report.py \\
        --group-id ${gid} \\
        --summary ${summary} \\
        --sex-check ${sex_check} \\
        --warnings ${warnings} \\
        --flagged ${flagged} \\
        ${check} \\
        --out ${gid}.html
    """
}
