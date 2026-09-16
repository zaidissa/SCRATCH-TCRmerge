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
| GEX `.h5ad` | Your annotated object (cell-type labels, embeddings, counts) |
| TCR table(s) | TCRtoolkit's `bridge/merged_vdj_object/post_qc_cells.tsv` (one row per cell) |

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
| `merged.tcr_only.h5ad` | Optional: only cells carrying TCR data (`subset_to_tcr`) |

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
(for the `.h5ad`). `.cirro/preprocess.py` finds one `.h5ad` and the TCR table(s) and wires them in.

Form options: TCR table filename, the sample/barcode/cell-id column names, the `.obs` prefix,
the minimum match rate, whether to write a TCR-only object, and whether to write `.rds` tables.

## Options

| Param | Default | Purpose |
|---|---|---|
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
