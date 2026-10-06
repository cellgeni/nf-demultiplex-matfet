// Write each library's labelled cells and summary from the group result,
// joined through barcode_provenance.tsv. Publishes to <outdir>/<sample_id>/.

process SPLIT_BY_LIBRARY {
    tag "${gmeta.id}"
    publishDir "${params.outdir}", mode: 'copy',
        saveAs: { fn -> fn.replaceFirst(/^per_library\//, '') }

    input:
    tuple val(gmeta), path(cell_assignments), path(provenance), path(gate)

    output:
    tuple val(gmeta), path("per_library/*/*"), emit: libraries

    script:
    """
    split_by_library.py \\
        --group-id ${gmeta.id} \\
        --cell-assignments ${cell_assignments} \\
        --provenance ${provenance} \\
        --gate ${gate} \\
        --outdir per_library
    """
}
