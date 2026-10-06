// Merge a group's per-library allele matrices into one group matrix, then run the
// independent integrity checks. Any failed check fails the task.
// The large .mtx/.vcf files stay in the work directory; only provenance and QC publish.
// Nextflow does not publish files declared inside an output directory, so the QC files
// are copied to the task root and published from there.

process BUILD_GROUP_MATRIX {
    tag "${gmeta.id}"
    publishDir "${params.outdir}/groups/${gmeta.id}/matrix", mode: 'copy',
        pattern: '*.{tsv,json,txt}'

    input:
    tuple val(gmeta), path(soup_dirs, stageAs: 'lib*/*')
    path panel_vcf

    output:
    tuple val(gmeta), path("${gmeta.id}"),                            emit: matrix
    tuple val(gmeta), path("${gmeta.id}/barcode_provenance.tsv"),     emit: provenance
    path "{barcode_provenance.tsv,library_offsets.tsv,matrix_qc.tsv,manifest.json}", emit: qc
    path "integrity_check.txt",                                       emit: check

    script:
    // Libraries arrive in samplesheet order; the builder concatenates cells in this order.
    def rows = [gmeta.samples, soup_dirs instanceof List ? soup_dirs : [soup_dirs]].transpose()
        .collect { sample_id, dir -> "'${gmeta.id},${sample_id},mixture,${gmeta.k},${dir}'" }
        .join(' ')
    """
    printf '%s\\n' 'group_id,sample_id,sample_role,expected_genotypes,souporcell_dir' ${rows} > groups.csv

    build_group_matrix.py \\
        --groups groups.csv \\
        --group-id ${gmeta.id} \\
        --outdir . \\
        --panel-vcf ${panel_vcf} \\
        --no-hash

    check_group_matrix.py \\
        --groups groups.csv \\
        --group-id ${gmeta.id} \\
        --outdir . \\
        --panel-vcf ${panel_vcf} \\
        > integrity_check.txt || { cat integrity_check.txt; exit 1; }
    cat integrity_check.txt

    cp ${gmeta.id}/{barcode_provenance.tsv,library_offsets.tsv,matrix_qc.tsv,manifest.json} .
    """
}
