// Thin wrapper over the pinned nf-demultiplex Souporcell path (git submodule).
// Heavy logic stays single-source in modules/nf-demultiplex/bin/. Only K is fixed here.

process SOUPORCELL {
    tag "${meta.id}"
    publishDir "${params.outdir}/souporcell", mode: 'copy'

    input:
    tuple val(meta), path(barcodes), path(bam)

    output:
    tuple val(meta), path("${meta.id}"), emit: dir
    path "versions.yml",                 emit: versions

    script:
    def soup_bin = "${projectDir}/modules/nf-demultiplex/bin"
    """
    # nf-demultiplex souporcell.sh needs a coordinate-sorted, indexed BAM alongside.
    if [ ! -f "${bam}.bai" ] && [ ! -f "${bam.baseName}.bai" ]; then
        samtools index -@ ${task.cpus} ${bam}
    fi

    # Positional args (from nf-demultiplex/bin/souporcell.sh):
    #   sample_id barcodes bam known_genotypes soc_vcf soc_fasta K skip_remap no_umi
    ${soup_bin}/souporcell.sh \\
        ${meta.id} \\
        ${barcodes} \\
        ${bam} \\
        ${params.known_genotypes} \\
        ${params.soc_vcf} \\
        ${params.soc_fasta} \\
        ${params.souporcell_k} \\
        ${params.skip_remap} \\
        ${params.no_umi}

    # Souporcell QC (ambient RNA %, doublets, per-cluster cell counts)
    ${soup_bin}/soup_qc.sh ${meta.id}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        souporcell_k: "${params.souporcell_k}"
        souporcell_skip_remap: "${params.skip_remap}"
        souporcell_no_umi: "${params.no_umi}"
        nf_demultiplex_bin: "${soup_bin}"
        samtools: "\$(samtools --version | head -n1 | sed 's/samtools //')"
    END_VERSIONS
    """
}
