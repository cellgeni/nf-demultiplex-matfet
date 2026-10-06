// Souporcell genotype discovery on a whole pregnancy group's matrix:
// clustering -> troublet (doublets) -> consensus.py (profiles + ambient RNA).
// The invocation reproduces the validation runs (docs/pipeline_design.md): seed 1,
// 4 threads, the binaries' own defaults for every other option.

process SOUPORCELL_DISCOVERY {
    tag "${gmeta.id}"
    publishDir "${params.outdir}/groups/${gmeta.id}/discovery", mode: 'copy',
        pattern: '{clusters.tsv,clusters_tmp.tsv,ambient_rna.txt,cluster_genotypes.vcf,*.err}'

    input:
    tuple val(gmeta), path(matrix)

    output:
    tuple val(gmeta), path("clusters.tsv"), path("ambient_rna.txt"), path("cluster_genotypes.vcf"), emit: results
    path "{clusters_tmp.tsv,*.err}", emit: logs
    path "versions.yml",             emit: versions

    script:
    // The pinned image has no `souporcell` executable on PATH (that name is a
    // directory); these are the real binaries inside params.container.
    def SOUPORCELL_BIN = '/opt/souporcell/souporcell/target/release/souporcell'
    def TROUBLET_BIN   = '/opt/souporcell/troublet/target/release/troublet'
    def CONSENSUS_BIN  = '/opt/souporcell/consensus.py'
    """
    ${SOUPORCELL_BIN} \\
        --alt_matrix ${matrix}/alt.mtx \\
        --ref_matrix ${matrix}/ref.mtx \\
        --barcodes ${matrix}/barcodes.tsv \\
        --num_clusters ${gmeta.k} \\
        --threads ${task.cpus} \\
        --seed ${params.discovery_seed} \\
        > clusters_tmp.tsv 2> clustering.err

    ${TROUBLET_BIN} \\
        --alts ${matrix}/alt.mtx \\
        --refs ${matrix}/ref.mtx \\
        --clusters clusters_tmp.tsv \\
        > clusters.tsv 2> troublet.err

    python ${CONSENSUS_BIN} \\
        -c clusters.tsv \\
        -a ${matrix}/alt.mtx \\
        -r ${matrix}/ref.mtx \\
        -v ${matrix}/variants.vcf \\
        --soup_out ambient_rna.txt \\
        --vcf_out cluster_genotypes.vcf \\
        --output_dir . \\
        > consensus.err 2>&1

    for f in clusters_tmp.tsv clusters.tsv ambient_rna.txt cluster_genotypes.vcf; do
        [ -s "\$f" ] || { echo "ERROR: \$f is empty; see the .err logs" >&2; exit 1; }
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        souporcell: "\$(${SOUPORCELL_BIN} --version 2>&1 | head -n1)"
        num_clusters: "${gmeta.k}"
        seed: "${params.discovery_seed}"
    END_VERSIONS
    """
}
