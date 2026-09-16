#!/usr/bin/env python3
"""
SCRATCH-TCRmerge: join per-cell TCR tables onto a GEX AnnData by cell barcode.

TCR data is per-cell METADATA, not a second expression matrix, so this does not
concatenate two objects: it adds TCR columns to the GEX object's .obs.

Scale note
----------
The expression matrix is never read. The input .h5ad is copied byte-for-byte and
only its /obs group is rewritten, so cost is independent of cell count and there
is no in-memory matrix (and therefore none of R's ~2.1e9 non-zero dgCMatrix cap,
which is why the join happens here and not in Seurat).

Barcode reconciliation
----------------------
VDJ and GEX barcode spellings routinely disagree ("AGAGAATGTACTACAA" vs
"AGAGAATGTACTACAA-1" vs "<sample>_AGAGAATGTACTACAA"). Rather than trusting one
convention, every candidate key is scored against the GEX obs_names and the best
match rate wins -- the same strategy TCRtoolkit uses internally. All candidates
and their rates are written to the join report so a bad join is visible, not silent.
"""

import argparse
import os
import shutil
import sys

import h5py
import numpy as np
import pandas as pd

try:                                     # anndata >= 0.11
    from anndata.io import read_elem, write_elem
except ImportError:                      # older anndata
    from anndata.experimental import read_elem, write_elem


# ── barcode key candidates ────────────────────────────────────────────────────
# Each returns a pd.Series of candidate join keys built from the TCR table.

def _bare(df, bc, sm):
    return df[bc].astype(str)


def _bare_minus1(df, bc, sm):
    return df[bc].astype(str) + "-1"


def _sample_bare(df, bc, sm):
    return df[sm].astype(str) + "_" + df[bc].astype(str)


def _sample_minus1(df, bc, sm):
    return df[sm].astype(str) + "_" + df[bc].astype(str) + "-1"


def _sample_dash(df, bc, sm):
    return df[sm].astype(str) + "-" + df[bc].astype(str)


def _cellid(df, bc, sm):
    return df["__cell_id__"].astype(str)


def _cellid_minus1(df, bc, sm):
    return df["__cell_id__"].astype(str) + "-1"


CANDIDATES = [
    ("barcode", _bare),
    ("barcode-1", _bare_minus1),
    ("sample_barcode", _sample_bare),
    ("sample_barcode-1", _sample_minus1),
    ("sample-barcode", _sample_dash),
    ("cell_id", _cellid),
    ("cell_id-1", _cellid_minus1),
]


def read_obs(h5ad_path):
    """Read only the /obs group - never touches X."""
    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            sys.exit(f"ERROR: {h5ad_path} has no /obs group - is it a valid .h5ad?")
        return read_elem(f["obs"])


def load_tcr_tables(paths, sample_col, barcode_col, cell_id_col):
    frames = []
    for p in paths:
        sep = "," if p.lower().endswith(".csv") else "\t"
        df = pd.read_csv(p, sep=sep, low_memory=False)
        df["__source_file__"] = os.path.basename(p)
        frames.append(df)
        print(f"[read] {os.path.basename(p)}: {len(df)} rows, {df.shape[1] - 1} columns")
    tcr = pd.concat(frames, ignore_index=True)

    for col, name in ((barcode_col, "barcode"), (sample_col, "sample")):
        if col not in tcr.columns:
            sys.exit(f"ERROR: TCR table has no '{col}' column (--{name}-col). "
                     f"Available: {list(tcr.columns)}")
    # cell_id is optional; synthesise so the cell_id candidates always work
    tcr["__cell_id__"] = (tcr[cell_id_col] if cell_id_col in tcr.columns
                          else tcr[sample_col].astype(str) + "_" + tcr[barcode_col].astype(str))
    return tcr


def score_candidates(tcr, obs_names, sample_col, barcode_col):
    obs_set = set(map(str, obs_names))
    rows = []
    for name, fn in CANDIDATES:
        try:
            keys = fn(tcr, barcode_col, sample_col)
        except Exception as exc:                      # a candidate may not apply
            rows.append({"candidate": name, "matched": 0, "total": len(tcr),
                         "match_rate": 0.0, "note": f"skipped: {exc}"})
            continue
        matched = int(sum(1 for k in keys if k in obs_set))
        rows.append({"candidate": name, "matched": matched, "total": len(tcr),
                     "match_rate": round(matched / max(len(tcr), 1), 6), "note": ""})
    report = pd.DataFrame(rows).sort_values("matched", ascending=False).reset_index(drop=True)
    return report


def main():
    ap = argparse.ArgumentParser(
        description="Join per-cell TCR tables onto a GEX .h5ad by barcode.")
    ap.add_argument("--gex-h5ad", required=True, help="GEX AnnData (.h5ad)")
    ap.add_argument("--tcr-table", action="append", default=[],
                    help="per-cell TCR table (.tsv/.csv); repeatable")
    ap.add_argument("--tcr-tables", default="",
                    help="comma-separated TCR tables (alternative to --tcr-table)")
    ap.add_argument("--out-h5ad", default="merged.h5ad")
    ap.add_argument("--tables-dir", default="tables")
    ap.add_argument("--sample-col", default="sample")
    ap.add_argument("--barcode-col", default="barcode")
    ap.add_argument("--cell-id-col", default="cell_id")
    ap.add_argument("--prefix", default="tcr_",
                    help="prefix for added obs columns (avoids clobbering GEX columns)")
    ap.add_argument("--min-match-rate", type=float, default=0.01,
                    help="fail if the best candidate matches fewer than this fraction "
                         "of TCR rows (default 0.01) - catches a wrong barcode convention")
    ap.add_argument("--subset-to-tcr", action="store_true",
                    help="ALSO write <out>.tcr_only.h5ad containing only TCR+ cells "
                         "(this one does read X for the subset)")
    args = ap.parse_args()

    paths = list(args.tcr_table)
    if args.tcr_tables:
        paths += [p for p in args.tcr_tables.split(",") if p.strip()]
    paths = [p for p in paths if p and os.path.exists(p) and os.path.getsize(p) > 0]
    if not paths:
        sys.exit("ERROR: no non-empty TCR tables given (--tcr-table/--tcr-tables)")

    os.makedirs(args.tables_dir, exist_ok=True)

    # ── 1. inputs ────────────────────────────────────────────────────────────
    tcr = load_tcr_tables(paths, args.sample_col, args.barcode_col, args.cell_id_col)
    obs = read_obs(args.gex_h5ad)
    print(f"[read] {os.path.basename(args.gex_h5ad)}: {len(obs)} cells, "
          f"{obs.shape[1]} obs columns (matrix NOT loaded)")

    # ── 2. pick the barcode convention by match rate ──────────────────────────
    report = score_candidates(tcr, obs.index, args.sample_col, args.barcode_col)
    print("\n[join] candidate barcode keys:")
    for _, r in report.iterrows():
        print(f"   {r['candidate']:<20} {r['matched']:>7} / {r['total']:<7} "
              f"({r['match_rate'] * 100:.2f}%) {r['note']}")
    report.to_csv(os.path.join(args.tables_dir, "barcode_join_report.tsv"),
                  sep="\t", index=False)

    best = report.iloc[0]
    if best["matched"] == 0 or best["match_rate"] < args.min_match_rate:
        sys.exit(
            f"\nERROR: best candidate '{best['candidate']}' matched only "
            f"{best['matched']}/{best['total']} TCR rows "
            f"({best['match_rate'] * 100:.2f}% < {args.min_match_rate * 100:.2f}%).\n"
            f"The barcode conventions do not line up. GEX obs_names look like: "
            f"{list(map(str, obs.index[:3]))}\n"
            f"TCR barcodes look like: {list(tcr[args.barcode_col].astype(str)[:3])}\n"
            f"See {args.tables_dir}/barcode_join_report.tsv for every candidate.")

    key_fn = dict(CANDIDATES)[best["candidate"]]
    tcr = tcr.copy()
    tcr["__key__"] = key_fn(tcr, args.barcode_col, args.sample_col)
    print(f"\n[join] using '{best['candidate']}' "
          f"({best['matched']}/{best['total']} TCR rows match a GEX cell)")

    # ── 3. join (one TCR row per cell; duplicates reported, first kept) ───────
    dup = int(tcr["__key__"].duplicated().sum())
    if dup:
        print(f"[warn] {dup} duplicate TCR keys; keeping the first of each")
        tcr.drop_duplicates(subset="__key__", keep="first", inplace=True)

    drop = {"__key__", "__cell_id__"}
    payload = tcr.set_index("__key__")[[c for c in tcr.columns if c not in drop]]
    payload = payload.add_prefix(args.prefix)

    new_obs = obs.join(payload, how="left")
    assert len(new_obs) == len(obs), "join changed the cell count - aborting"
    new_obs = new_obs.loc[obs.index]                       # preserve original order

    matched_col = f"{args.prefix}has_tcr_match"
    probe = f"{args.prefix}{args.barcode_col}"
    new_obs[matched_col] = new_obs[probe].notna() if probe in new_obs else False
    n_matched = int(new_obs[matched_col].sum())
    print(f"[join] {n_matched} of {len(new_obs)} GEX cells carry TCR data "
          f"({n_matched / max(len(new_obs), 1) * 100:.2f}%)")

    # h5-safe dtypes: object -> str, NaN -> ""
    for c in new_obs.columns:
        if new_obs[c].dtype == object:
            new_obs[c] = new_obs[c].astype(str).replace({"nan": "", "None": "", "<NA>": ""})

    # ── 4. write merged h5ad by copying the file and replacing /obs ───────────
    shutil.copyfile(args.gex_h5ad, args.out_h5ad)
    with h5py.File(args.out_h5ad, "r+") as f:
        del f["obs"]
        write_elem(f, "obs", new_obs)
    print(f"[write] {args.out_h5ad} (X copied untouched)")

    # ── 5. side tables, for looking at things outside Python ─────────────────
    td = args.tables_dir
    new_obs.reset_index(names="cell_id").to_csv(
        os.path.join(td, "merged_obs.tsv"), sep="\t", index=False)

    unmatched = tcr[~tcr["__key__"].isin(set(map(str, obs.index)))]
    unmatched.drop(columns=["__key__", "__cell_id__"], errors="ignore").to_csv(
        os.path.join(td, "tcr_rows_without_gex_cell.tsv"), sep="\t", index=False)

    if args.sample_col in tcr.columns:
        per_sample = (tcr.assign(matched=tcr["__key__"].isin(set(map(str, obs.index))))
                        .groupby(args.sample_col)["matched"]
                        .agg(tcr_cells="size", matched_to_gex="sum")
                        .reset_index())
        per_sample["match_rate"] = (per_sample["matched_to_gex"] /
                                    per_sample["tcr_cells"].clip(lower=1)).round(6)
        per_sample.to_csv(os.path.join(td, "per_sample_join_summary.tsv"),
                          sep="\t", index=False)
        print("\n[summary] per sample:")
        for _, r in per_sample.iterrows():
            print(f"   {str(r[args.sample_col]):<14} {int(r['matched_to_gex']):>6}"
                  f" / {int(r['tcr_cells']):<6} ({r['match_rate'] * 100:.1f}%)")

    pd.DataFrame([{
        "gex_h5ad": os.path.basename(args.gex_h5ad),
        "gex_cells": len(obs),
        "tcr_tables": ";".join(os.path.basename(p) for p in paths),
        "tcr_rows": len(tcr),
        "barcode_key": best["candidate"],
        "cells_with_tcr": n_matched,
        "tcr_rows_unmatched": len(unmatched),
    }]).to_csv(os.path.join(td, "merge_summary.tsv"), sep="\t", index=False)

    # ── 6. optional TCR-only subset (reads X for the kept cells only) ────────
    if args.subset_to_tcr:
        import anndata as ad
        sub_path = args.out_h5ad.replace(".h5ad", "") + ".tcr_only.h5ad"
        adata = ad.read_h5ad(args.out_h5ad, backed="r")
        keep = np.asarray(new_obs[matched_col].values, dtype=bool)
        adata[keep].to_memory().write_h5ad(sub_path)
        print(f"[write] {sub_path} ({int(keep.sum())} TCR+ cells)")

    print("\nDone.")


if __name__ == "__main__":
    main()
