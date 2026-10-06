// Decide whether the group's second Souporcell genotype is real (two_profiles),
// too small to use (one_profile), or not supported (unresolved).
// Thresholds and their evidence: docs/pipeline_design.md, "Composition gate".

process COMPOSITION_GATE {
    tag "${gmeta.id}"
    publishDir "${params.outdir}/groups/${gmeta.id}", mode: 'copy'

    input:
    tuple val(gmeta), path(clusters), path(ambient), path(profiles), path(provenance)

    output:
    tuple val(gmeta), path("gate.tsv"), path("gate_libraries.tsv"), emit: gate

    script:
    """
    composition_gate.py \\
        --group-id ${gmeta.id} \\
        --clusters ${clusters} \\
        --ambient ${ambient} \\
        --profiles ${profiles} \\
        --provenance ${provenance} \\
        --min-discordant-fraction ${params.gate_min_discordant_fraction} \\
        --min-minority-singlets ${params.gate_min_minority_singlets} \\
        --out gate.tsv \\
        --libraries-out gate_libraries.tsv
    """
}
