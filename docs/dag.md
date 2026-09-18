# Pipeline DAG

Two views of the same pipeline: a curated diagram that names the files moving between
processes, and the raw graph Nextflow emits from the actual channel topology.

## Regenerating

Nextflow builds the DAG from the workflow definition without running anything
(`-preview`), so no real data is needed — only paths that exist. Formats are chosen by
the output extension (`.mmd`, `.html`, `.dot`, `.svg`, `.png`; the last two need `dot`
on PATH).

```bash
module load cellgen/nextflow

# Full pipeline (SOUPORCELL + labeling). Needs a samplesheet whose paths exist;
# the files are never read under -preview, so empty stand-ins are fine.
nextflow run . -preview --input samplesheet.csv -with-dag docs/dag.svg

# Labeling-only entry, against the shipped fixture
nextflow run . -entry LABEL_ONLY -profile test -preview -with-dag docs/dag_label.mmd
```

`nextflow.config` also has `dag.enabled = true`, so every real run drops its own
`${params.outdir}/pipeline_info/dag.html`.

## Curated view — full pipeline

```mermaid
flowchart TB
    sheet[/"samplesheet.csv<br/>sample_id,bam,barcodes,expression"/]
    panel[/"assets/marker_panel.tsv"/]

    sheet -->|splitCsv + validate| rows(("ch_rows<br/>meta, barcodes, bam, expression"))
    rows -->|"toList"| dup{{"duplicate sample_id check"}}

    subgraph SD["SOUPORCELL_DEMUX"]
        SOUP(["SOUPORCELL<br/>K=2, samtools index + soup_qc.sh"])
    end

    subgraph LMF["LABEL_MATERNAL_FETAL"]
        MS(["MARKER_SUMMARY<br/>aggregate_markers.py"])
        LC(["LABEL_CLUSTERS<br/>label_clusters.py"])
        BR(["BUILD_REPORT<br/>make_report.py"])
    end

    rows -->|"meta, barcodes, bam"| SOUP
    SOUP -->|"versions.yml"| vers[/"versions.yml"/]
    SOUP -->|"souporcell/&lt;id&gt;/<br/>clusters.tsv, ambient_rna.txt, qc.txt"| joinA(("join on meta"))
    rows -->|"meta, expression"| joinA
    panel --> MS
    joinA -->|ch_label_in| MS

    MS -->|gene_resolution.tsv| res[/"gene_resolution.tsv"/]
    MS -->|"cluster_marker_summary.tsv + soup dir"| LC
    LC -->|cluster_labels.tsv, sample_qc_summary.tsv| joinB(("join on meta"))
    MS --> joinB
    LC -->|cell_assignments.tsv.gz| assign[/"cell_assignments.tsv.gz"/]
    joinB --> BR
    BR --> rep[/"reports/&lt;sample_id&gt;.html"/]
```

## Curated view — `-entry LABEL_ONLY`

Skips Souporcell entirely and takes existing Souporcell output directories from the
samplesheet, which is what makes `LABEL_MATERNAL_FETAL` independently runnable
(ROADMAP §4).

```mermaid
flowchart TB
    sheet[/"samplesheet.csv<br/>sample_id,souporcell_dir,expression"/]
    panel[/"assets/marker_panel.tsv"/]
    sheet -->|splitCsv + validate| MS(["MARKER_SUMMARY"])
    panel --> MS
    MS --> LC(["LABEL_CLUSTERS"])
    MS --> BR(["BUILD_REPORT"])
    LC --> BR
    LC --> assign[/"cell_assignments.tsv.gz"/]
    BR --> rep[/"reports/&lt;sample_id&gt;.html"/]
```

## Nextflow-generated view

Rendered from the real channel topology, so it also shows the operator nodes
(`map`, `join`, `toList`) between processes: [dag.svg](dag.svg) / [dag.png](dag.png).

![Nextflow DAG](dag.png)
