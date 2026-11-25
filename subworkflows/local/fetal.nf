// FETAL — fetal mode: flag maternal cells in libraries from one fetus, with no maternal
// sample (docs/fetal_mode_design.md).
//   SOUPORCELL (k=2, the checked libraries: by default the first of each batch)
//     -> SOUPORCELL_DISTANCE (whole run: cluster-vs-cluster genotype distance, heatmaps, warnings)
//   CELLSNP_PILEUP (per library, unless cellsnp_dir is given) -> FETAL_SCORE (per batch)
//   SEX_COUNTS (per library) -> FETAL_SUMMARY (per fetus) -> FETAL_REPORT (per fetus)
// The Souporcell check runs alongside scoring and only reports: it never changes fetal mode's calls.
// Batches are fixed in main.nf from the samplesheet: each group's libraries in samplesheet
// order, cut into runs of params.fetal_batch_size.

include { SOUPORCELL          } from '../../modules/local/souporcell'
include { SOUPORCELL_DISTANCE } from '../../modules/local/souporcell_distance'
include { CELLSNP_PILEUP      } from '../../modules/local/cellsnp_pileup'
include { FETAL_SCORE         } from '../../modules/local/fetal_score'
include { SEX_COUNTS          } from '../../modules/local/sex_counts'
include { FETAL_SUMMARY       } from '../../modules/local/fetal_summary'
include { FETAL_REPORT        } from '../../modules/local/fetal_report'

workflow FETAL {
    take:
    ch_libraries  // tuple(meta, barcodes, bam, expression, cellsnp_dir or null)
                  // meta: id, group, index, group_size, batch, batch_size, group_batches
    ch_souporcell // tuple(meta, barcodes, bam, souporcell_dir or null): libraries for the Souporcell check
    panel_vcf     // path: common-variant panel with population allele frequencies

    main:
    // ---- Souporcell k=2 on the checked libraries, then one cross-library comparison for the run ----
    // Only the sample id goes in, as in pregnancy mode, so the cached run depends on the library alone.
    SOUPORCELL(
        ch_souporcell
            .filter { _meta, _barcodes, _bam, soup_dir -> soup_dir == null }
            .map { meta, barcodes, bam, _soup_dir -> [ [ id: meta.id ], barcodes, bam ] }
    )
    ch_soup_group = ch_souporcell.map { meta, _barcodes, _bam, _soup_dir -> [ meta.id, meta.group, meta.index ] }
    ch_soup_dirs = SOUPORCELL.out.dir
        .map { lib, dir -> [ lib.id, dir ] }
        .mix(
            ch_souporcell
                .filter { _meta, _barcodes, _bam, soup_dir -> soup_dir != null }
                .map { meta, _barcodes, _bam, soup_dir -> [ meta.id, soup_dir ] }
        )
        .join(ch_soup_group)
        .toSortedList { a, b -> a[3] <=> b[3] }                  // samplesheet order
        .filter { it }
        .map { rows -> [ rows.collect { it[2] }, rows.collect { it[0] }, rows.collect { it[1] } ] }
    SOUPORCELL_DISTANCE(ch_soup_dirs)

    // Only the sample id goes into the pileup, so its cache depends on the library alone:
    // changing the batching or the groups does not redo hours of pileup on -resume.
    CELLSNP_PILEUP(
        ch_libraries
            .filter { _meta, _barcodes, _bam, _expression, cellsnp_dir -> cellsnp_dir == null }
            .map { meta, barcodes, bam, _expression, _cellsnp_dir -> [ [ id: meta.id ], barcodes, bam ] },
        panel_vcf
    )

    ch_pileups = CELLSNP_PILEUP.out.dir
        .map { lib, dir -> [ lib.id, dir ] }
        .mix(
            ch_libraries
                .filter { _meta, _barcodes, _bam, _expression, cellsnp_dir -> cellsnp_dir != null }
                .map { meta, _barcodes, _bam, _expression, cellsnp_dir -> [ meta.id, cellsnp_dir ] }
        )

    ch_meta = ch_libraries.map { meta, _barcodes, _bam, expression, _cellsnp_dir -> [ meta.id, meta, expression ] }

    // Each batch moves on as soon as its own libraries are piled up, in samplesheet order.
    ch_batches = ch_pileups
        .join(ch_meta)
        .map { id, dir, meta, expression ->
            [ groupKey(meta.batch, meta.batch_size), meta.index, id, dir, expression, meta ] }
        .groupTuple()
        .map { key, indices, ids, dirs, expressions, metas ->
            def o = [ indices, ids, dirs, expressions ].transpose().sort { it[0] }
            [ [ id: key.toString(), group: metas[0].group, group_batches: metas[0].group_batches ],
              o.collect { it[1] }, o.collect { it[2] }, o.collect { it[3] } ]
        }

    FETAL_SCORE(ch_batches, panel_vcf)

    SEX_COUNTS(ch_libraries.map { meta, _barcodes, _bam, expression, _cellsnp_dir -> [ [ id: meta.id ], expression ] })

    // Per fetus: wait for all of its batches and all of its libraries' sex counts.
    ch_fetus_scores = FETAL_SCORE.out.scores
        .map { bmeta, cells, _summary, _removed -> [ groupKey(bmeta.group, bmeta.group_batches), bmeta.id, cells ] }
        .groupTuple()
        .map { key, ids, cells -> [ key.toString(), [ ids, cells ].transpose().sort { it[0] }.collect { it[1] } ] }

    ch_fetus_sex = SEX_COUNTS.out.counts
        .map { lib, counts -> [ lib.id, counts ] }
        .join(ch_meta)
        .map { _id, counts, meta, _expression -> [ groupKey(meta.group, meta.group_size), meta.index, counts ] }
        .groupTuple()
        .map { key, indices, counts -> [ key.toString(), [ indices, counts ].transpose().sort { it[0] }.collect { it[1] } ] }

    FETAL_SUMMARY(ch_fetus_scores.join(ch_fetus_sex))
    // Every report waits for the run's Souporcell check, or gets none when it is switched off.
    ch_report = params.fetal_souporcell == 'off'
        ? FETAL_SUMMARY.out.summary.map { it + [ [] ] }
        : FETAL_SUMMARY.out.summary.combine(SOUPORCELL_DISTANCE.out.dir)
    FETAL_REPORT(ch_report)

    emit:
    scores   = FETAL_SCORE.out.scores        // tuple(bmeta, cells.tsv.gz, summary.tsv, removed_sites.tsv)
    summary  = FETAL_SUMMARY.out.summary     // tuple(group_id, fetal_summary, sex_check, warnings, flagged_cells)
    check    = SOUPORCELL_DISTANCE.out.dir   // path: souporcell_check/
    report   = FETAL_REPORT.out.report
    versions = CELLSNP_PILEUP.out.versions.mix(FETAL_SCORE.out.versions, SOUPORCELL.out.versions)
}
