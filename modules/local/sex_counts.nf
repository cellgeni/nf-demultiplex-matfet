// Fetal mode, per library: total, XIST and Y-gene UMIs per cell from the expression matrix,
// for the per-fetus sex call and the genotype-free XIST/Y check (bin/sex_counts.py).

process SEX_COUNTS {
    tag "${meta.id}"

    input:
    tuple val(meta), path(expression)

    output:
    tuple val(meta), path("${meta.id}.sex_counts.tsv.gz"), emit: counts

    script:
    """
    sex_counts.py --library ${meta.id} --expression ${expression} --prefix ${meta.id}
    """
}
