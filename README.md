# SCRATCH-TCRmerge

> **Not to be confused with `SCRATCH-TCR`**, which analyses TCR data. SCRATCH-TCRmerge
> joins TCRtoolkit's per-cell TCR output onto a SCRATCH GEX object (`.h5ad`) by cell barcode.

Joins per-cell **TCR** data onto an annotated **GEX** AnnData by cell barcode, and emits a
merged `.h5ad` plus the joined tables as `.tsv` and `.rds`.

TCR data is per-cell *metadata*, not a second expression matrix, so this does not concatenate
two objects — it adds TCR columns to the GEX object's `.obs`.

## Why this is a separate tool

TCRtoolkit's single-cell route needs the GEX object as a Seurat `.rds`. R's sparse matrices
(`dgCMatrix`) index non-zeros with a 32-bit signed integer, so a single matrix caps out at
**2,147,483,647 non-zeros** — roughly 430K–1.1M cells at typical density. Large cohorts cannot
become a Seurat object at all.

This tool sidesteps that entirely: **the expression matrix is never read**. The input `.h5ad` is
copied byte-for-byte and only its `/obs` group is rewritten, so runtime and memory are independent
of cell count. The GEX side stays in Python; TCRtoolkit keeps taking `.rds` exactly as before.

## Inputs

| Input | What it is |
|---|---|
| GEX `.h5ad` | SCRATCH-QC's per-run object, e.g. `GBM_DFCI1_CSF_singlet.h5ad` — selected by the `gex_h5ad_pattern` glob (default `*_singlet.h5ad`) since the name changes every run |
| TCR table(s) | TCRtoolkit's `bridge/merged_vdj_object/post_qc_cells.tsv` (one row per cell) |

**Which GEX object?** SCRATCH-QC's `*_singlet.h5ad` carries QC metadata and embeddings but **not**
cell-type labels — those come from SCRATCH-Annotation (celltypist writes
`*_celltypist_annotation_object.h5ad` with `celltypist_cell_label_coarse`). If you need subset
classification (CD8 effector, CD4 Treg…) alongside the TCR data, point `gex_h5ad_pattern` at the
annotation object instead.

## Outputs

| Output | Contents |
|---|---|
| `merged.h5ad` | The GEX object with TCR columns added to `.obs` (prefixed `tcr_` by default) |
| `tables/merged_obs.tsv` | The full joined per-cell metadata |
| `tables/barcode_join_report.tsv` | Every barcode convention tried, with its match rate |
| `tables/per_sample_join_summary.tsv` | Cells matched per sample |
| `tables/tcr_rows_without_gex_cell.tsv` | TCR cells with no matching GEX cell |
| `tables/merge_summary.tsv` | One-line summary of the join |
| `tables_rds/*.rds` | The same tables as R data.frames (`emit_rds`, default on) |
| `tcr_source/*` | The rest of the TCRtoolkit bundle, copied verbatim (see below) |
| `merged.tcr_only.h5ad` | Optional: only cells carrying TCR data (`subset_to_tcr`) |

### What gets joined, and what is only carried over

TCRtoolkit's `bridge/merged_vdj_object/` holds eight files. Only one is joinable:

| File | Treatment |
|---|---|
| `post_qc_cells.tsv` | **Joined** onto `.obs` — one row per cell |
| `pre_qc_cells.tsv` | Carried over. Same barcodes as post-QC, so joining both would duplicate every cell |
| `pre/post_qc_summary.tsv` | Carried over. Per-**sample** (one row per sample), not per-cell |
| `pre/post_qc_seurat.rds`, `pre/post_qc_combineTCR.rds` | Carried over. R objects the Python merger cannot read |

Carried-over files land in `tcr_source/` so the merged dataset holds the whole TCR bundle
alongside the merged object. Set "Other TCR files to carry over" to empty to skip them.

The `.rds` files are **data.frames, not Seurat objects** — deliberately. Rebuilding a Seurat
object here would reintroduce the matrix-size ceiling this tool exists to avoid.

## Barcode reconciliation

VDJ and GEX barcode spellings routinely disagree. TCRtoolkit emits a bare 16-mer
(`AGAGAATGTACTACAA`) with no `-1` suffix; GEX objects usually carry `-1`, a sample prefix, or both.

Rather than assuming a convention, every candidate key is scored against the GEX `obs_names` and
the best match rate wins:

```
barcode              barcode-1            sample_barcode
sample_barcode-1     sample-barcode       cell_id            cell_id-1
```

All candidates and their rates land in `tables/barcode_join_report.tsv`. If the best candidate
matches fewer than `min_match_rate` of TCR rows (default 1%), the run **fails with both formats
printed** rather than silently writing an object with no TCR data.

## Running locally

```bash
nextflow run . \
  --gex_h5ad   /path/to/annotated_gex.h5ad \
  --tcr_tables /path/to/post_qc_cells.tsv \
  --outdir     results
```

Multiple TCR tables: `--tcr_tables 'a_cells.tsv,b_cells.tsv'`.

Without Nextflow, the script stands alone:

```bash
python bin/merge_tcr_gex.py \
  --gex-h5ad annotated_gex.h5ad \
  --tcr-table post_qc_cells.tsv \
  --out-h5ad merged.h5ad \
  --tables-dir tables
```

## Running on Cirro

Register this repo as a Cirro process pointing at `.cirro/`, then select **both** datasets under
"Datasets to use": the TCRtoolkit output dataset (for `post_qc_cells.tsv`) and the GEX dataset
(for the `.h5ad`).

Two ways to name the inputs, in priority order:

1. **Pick them explicitly** — "GEX object (.h5ad)" and "TCR per-cell table (.tsv)" are file
   pickers that browse the selected datasets. Recommended: you see exactly which file is used.
2. **Leave the pickers empty** and `.cirro/preprocess.py` discovers them by filename pattern
   (`*_singlet.h5ad` and `post_qc_cells.tsv` by default), logging what it chose.

Either way, keep the TCR dataset selected: the carried-over files in `tcr_source/` are found by
scanning the selected datasets, not through a picker.

The process must have **`uses_sample_sheet` disabled** in its Cirro definition. This pipeline
consumes files, not samples, and SCRATCH-QC output datasets carry no sample metadata — leaving it
enabled makes Cirro demand a sample selection that cannot be satisfied.

Form options: TCR table filename, the sample/barcode/cell-id column names, the `.obs` prefix,
the minimum match rate, whether to write a TCR-only object, and whether to write `.rds` tables.

## Options

| Param | Default | Purpose |
|---|---|---|
| `gex_h5ad_pattern` | `*_singlet.h5ad` | Glob picking the GEX object (Cirro only; locally use `--gex_h5ad`) |
| `tcr_table_pattern` | `post_qc_cells.tsv` | Which TCR table to take (Cirro only; locally use `--tcr_tables`) |
| `--sample_col` | `sample` | Sample column in the TCR table |
| `--barcode_col` | `barcode` | Barcode column in the TCR table |
| `--cell_id_col` | `cell_id` | Composite id; synthesised if absent |
| `--obs_prefix` | `tcr_` | Prefix for added `.obs` columns |
| `--min_match_rate` | `0.01` | Abort below this barcode match rate |
| `--subset_to_tcr` | `false` | Also write a TCR-only `.h5ad` |
| `--emit_rds` | `true` | Also write `.rds` tables |
| `--container` | `syedsazaidi/scratch-tcr:latest` | Ships anndata, h5py, pandas, R/Seurat |

## Checking the result

```python
import anndata as ad
a = ad.read_h5ad("merged.h5ad", backed="r")
print(a.obs["tcr_has_tcr_match"].sum(), "of", a.n_obs, "cells have TCR data")
print(a.obs.filter(like="tcr_").columns.tolist())
```

Then cross-check `tables/barcode_join_report.tsv`: if the winning candidate's rate is far below
what you expect, the barcode convention is wrong, not the data.
