// GROUP_DEMUX — infer the genotypes once per pregnancy group and gate the result.
//   BUILD_GROUP_MATRIX -> SOUPORCELL_DISCOVERY -> COMPOSITION_GATE

include { BUILD_GROUP_MATRIX   } from '../../modules/local/build_group_matrix'
include { SOUPORCELL_DISCOVERY } from '../../modules/local/souporcell_discovery'
include { COMPOSITION_GATE     } from '../../modules/local/composition_gate'

workflow GROUP_DEMUX {
    take:
    ch_groups    // tuple(gmeta, [souporcell_dir, ...]) in samplesheet order
    panel_vcf    // path: reference panel; only its ##contig header is read

    main:
    BUILD_GROUP_MATRIX(ch_groups, panel_vcf)
    SOUPORCELL_DISCOVERY(BUILD_GROUP_MATRIX.out.matrix)

    COMPOSITION_GATE(
        SOUPORCELL_DISCOVERY.out.results.join(BUILD_GROUP_MATRIX.out.provenance)
    )

    emit:
    discovery  = SOUPORCELL_DISCOVERY.out.results       // tuple(gmeta, clusters, ambient, profiles)
    provenance = BUILD_GROUP_MATRIX.out.provenance      // tuple(gmeta, barcode_provenance.tsv)
    gate       = COMPOSITION_GATE.out.gate              // tuple(gmeta, gate.tsv, gate_libraries.tsv)
    versions   = SOUPORCELL_DISCOVERY.out.versions
}
