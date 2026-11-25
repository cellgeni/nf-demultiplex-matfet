// Fetal mode, whole run: genotype distance between every pair of Souporcell clusters of the checked
// libraries, per fetus and across fetuses, with heatmaps for the reports (bin/souporcell_distance.py).
// A report-only check: it warns for review and never changes fetal mode's calls.

process SOUPORCELL_DISTANCE {
    publishDir "${params.outdir}/fetal", mode: 'copy'

    input:
    tuple val(groups), val(libraries), path(dirs, stageAs: 'soc*/*')

    output:
    path "souporcell_check", emit: dir

    script:
    // One manifest row per library: group, library, staged Souporcell directory.
    def rows = [ groups, libraries, [ dirs ].flatten() ].transpose().collect { g, lib, d -> "'${g}' '${lib}' '${d}'" }.join(' ')
    """
    printf '%s\\t%s\\t%s\\n' group library dir ${rows} > manifest.tsv
    souporcell_distance.py \\
        --manifest manifest.tsv \\
        --min-singlets ${params.fetal_souporcell_min_singlets} \\
        --warn-msd ${params.fetal_souporcell_warn_msd} \\
        --outdir souporcell_check
    """
}
