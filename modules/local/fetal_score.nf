// Fetal mode, per batch of up to params.fetal_batch_size libraries of one fetus: pooled
// fetal genotype and per-cell maternal score (bin/fetal_score.py, the frozen prototype v3).
// Publishes to <outdir>/fetal/<group_id>/batches/.

process FETAL_SCORE {
    tag "${bmeta.id}"
    publishDir "${params.outdir}/fetal/${bmeta.group}/batches", mode: 'copy', pattern: "${bmeta.id}.*"

    input:
    tuple val(bmeta), val(libraries), path(cellsnp, stageAs: 'pileup*/*'), path(expression, stageAs: 'expression*/*')
    path panel_vcf

    output:
    tuple val(bmeta), path("${bmeta.id}.cells.tsv.gz"), path("${bmeta.id}.summary.tsv"),
          path("${bmeta.id}.removed_sites.tsv"), emit: scores
    path "versions.yml",                       emit: versions

    script:
    """
    # Staged as pileup1/<dir>, pileup2/<dir>, ... (expression likewise) in batch order.
    pileups=( \$(ls -d pileup*/* | sort -V) )
    exprs=( \$(ls -d expression*/* | sort -V) )
    libs=( ${libraries.join(' ')} )
    args=()
    for i in "\${!libs[@]}"; do
        args+=( --library "\${libs[\$i]}=\${pileups[\$i]}" --expression "\${libs[\$i]}=\${exprs[\$i]}" )
    done

    fetal_score.py "\${args[@]}" \\
        --panel ${panel_vcf} \\
        --min-sites ${params.fetal_min_sites} \\
        --rhos ${params.fetal_rhos} \\
        --prefix ${bmeta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fetal_min_sites: "${params.fetal_min_sites}"
        fetal_rhos: "${params.fetal_rhos}"
        python: "\$(python3 --version | sed 's/Python //')"
        numpy: "\$(python3 -c 'import numpy; print(numpy.__version__)')"
        scipy: "\$(python3 -c 'import scipy; print(scipy.__version__)')"
    END_VERSIONS
    """
}
