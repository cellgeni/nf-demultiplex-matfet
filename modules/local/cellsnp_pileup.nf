// Fetal mode, per library: allele counts at the common-variant panel with cellsnp-lite,
// read straight from the CellRanger/STARsolo BAM (no remapping). Flags as in the
// prototype runs that produced the embryo results and the D8 benchmark.

process CELLSNP_PILEUP {
    tag "${meta.id}"

    input:
    tuple val(meta), path(barcodes), path(bam)
    path panel_vcf

    output:
    tuple val(meta), path("cellsnp_${meta.id}"), emit: dir
    path "versions.yml",                         emit: versions

    script:
    """
    # Use the index next to the source BAM when there is one; index only if missing.
    src=\$(readlink -f ${bam})
    for idx in "\${src}.bai" "\${src%.bam}.bai"; do
        if [ -f "\${idx}" ]; then ln -sf "\${idx}" ${bam}.bai; break; fi
    done
    [ -f ${bam}.bai ] || samtools index ${bam}    # single-threaded: -@ hangs with samtools 1.9 in this image

    if [[ ${barcodes} == *.gz ]]; then zcat ${barcodes} > cells.tsv; else cat ${barcodes} > cells.tsv; fi

    cellsnp-lite \\
        -s ${bam} -b cells.tsv -R ${panel_vcf} -O cellsnp_${meta.id} \\
        -p ${task.cpus} --cellTAG CB --UMItag UB --minMAF 0 --minCOUNT 1 --gzip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cellsnp-lite: "\$(cellsnp-lite --version 2>&1 | grep -oE '[0-9]+\\.[0-9]+\\.[0-9]+' | head -n1)"
    END_VERSIONS
    """
}
