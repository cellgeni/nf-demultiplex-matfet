// Fetal mode, per fetus: fetal_mode per cell, the pooled bound, the in-range checks and
// warnings (bin/fetal_summary.py). Publishes the fetus's tables to <outdir>/fetal/<group_id>/
// and each library's calls to <outdir>/<sample_id>/fetal_mode_cells.tsv.gz.

process FETAL_SUMMARY {
    tag "${gid}"
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> fn.startsWith('summary/per_library/') ? fn.replaceFirst(/^summary\/per_library\//, '')
                                                              : fn.replaceFirst(/^summary\//, "fetal/${gid}/") }

    input:
    tuple val(gid), path(cells), path(sex_counts)

    output:
    tuple val(gid), path("summary/fetal_summary.tsv"), path("summary/sex_check.tsv"), path("summary/warnings.txt"),
          path("summary/flagged_cells.tsv"), emit: summary
    path "summary/cells.tsv.gz",              emit: cells
    path "summary/per_library/*/*",           emit: libraries

    script:
    """
    fetal_summary.py \\
        --group-id ${gid} \\
        --cells ${cells} \\
        --sex-counts ${sex_counts} \\
        --warn-flagged-fraction ${params.fetal_warn_flagged_fraction} \\
        --max-maternal-share ${params.fetal_max_maternal_share} \\
        --outdir summary
    """
}
