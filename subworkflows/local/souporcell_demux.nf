// SOUPORCELL_DEMUX — run Souporcell at K=2 and publish unmodified outputs.
// Single source of truth for the Souporcell invocation (wraps the nf-demultiplex submodule).

include { SOUPORCELL } from '../../modules/local/souporcell'

workflow SOUPORCELL_DEMUX {
    take:
    ch_samples   // tuple(meta, barcodes, bam)

    main:
    SOUPORCELL(ch_samples)

    emit:
    dir      = SOUPORCELL.out.dir       // tuple(meta, souporcell_dir)
    versions = SOUPORCELL.out.versions
}
