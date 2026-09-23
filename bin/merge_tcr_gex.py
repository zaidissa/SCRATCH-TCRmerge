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
match rate wins. All candidates and their rates are written to the join report so
a bad join is visible, not silent.

Library-type sample names
-------------------------
Paired 5' GEX and VDJ libraries from the same GEM well share cell barcodes but are
usually named after their library type: BTC-GBM-001-001-GEX vs BTC-GBM-001-001-TCR,
or ...-S3-GEX-LIB vs ...-S3-TCR-LIB. Exact sample-prefixed keys can never match
those. So the library words (--library-tokens, default GEX,TCR,VDJ,BCR,ADT) are
dropped as whole tokens from the sample part of BOTH sides, and the result is tried
as extra candidates. Matches are mapped back to the ORIGINAL GEX cell names, exact
candidates win ties, and normalised keys that would collide are excluded rather
than guessed. Barcodes are never joined without their sample: 10x barcodes collide
across libraries by chance, so a bare-barcode join would attach one sample's TCR to
another sample's cell.
"""

import argparse
import os
import re
import shutil
import sys

import h5py
import numpy as np
import pandas as pd

try:                                     # anndata >= 0.11
    from anndata.io import read_elem, write_elem
except ImportError:                      # older anndata
    from anndata.experimental import read_elem, write_elem


# ── sample-name normalisation ─────────────────────────────────────────────────

def make_normaliser(tokens_csv):
    """Return a function dropping library-type words from a sample name, or None."""
    tokens = {t.strip().upper() for t in (tokens_csv or "").split(",") if t.strip()}
    if not tokens:
        return None

    def normalise(sample):
        parts = re.split(r"[-_]", str(sample))
        kept = [p for p in parts if p.upper() not in tokens]
        # A name made only of library words would normalise to nothing - keep it as is.
        return "-".join(kept) if kept else str(sample)

    return normalise


def build_normalised_index(obs_names, normalise):
    """
    Map a normalised GEX key -> original obs_name. obs_names are split at the LAST
    underscore into <sample>_<barcode...>; only the sample part is normalised. Keys
    two different cells would share are dropped, never resolved by guessing.
    """
    index, clashes = {}, set()
    for name in map(str, obs_names):
        if "_" not in name:
            continue
        sample, tail = name.rsplit("_", 1)
        key = f"{normalise(sample)}_{tail}"
        if key in index and index[key] != name:
            clashes.add(key)
        else:
            index[key] = name
    for key in clashes:
        index.pop(key, None)
    return index, len(clashes)


# ── barcode key candidates ────────────────────────────────────────────────────
# (name, builder, target). builder(df, barcode_col, sample_col) -> pd.Series of keys.
# target "exact" keys are matched against obs_names as written; "normalised" keys
# against the library-normalised GEX index.

def build_candidates(normalise):
    cands = [
        ("barcode",          lambda d, bc, sm: d[bc].astype(str), "exact"),
        ("barcode-1",        lambda d, bc, sm: d[bc].astype(str) + "-1", "exact"),
        ("sample_barcode",   lambda d, bc, sm: d[sm].astype(str) + "_" + d[bc].astype(str), "exact"),
        ("sample_barcode-1", lambda d, bc, sm: d[sm].astype(str) + "_" + d[bc].astype(str) + "-1", "exact"),
        ("sample-barcode",   lambda d, bc, sm: d[sm].astype(str) + "-" + d[bc].astype(str), "exact"),
        ("cell_id",          lambda d, bc, sm: d["__cell_id__"].astype(str), "exact"),
        ("cell_id-1",        lambda d, bc, sm: d["__cell_id__"].astype(str) + "-1", "exact"),
    ]
    if normalise:
        cands += [
            ("libnorm:sample_barcode-1",
             lambda d, bc, sm: d[sm].astype(str).map(normalise) + "_" + d[bc].astype(str) + "-1",
             "normalised"),
            ("libnorm:sample_barcode",
             lambda d, bc, sm: d[sm].astype(str).map(normalise) + "_" + d[bc].astype(str),
             "normalised"),
        ]
    return cands


def read_obs(h5ad_path):
    """Read only the /obs group - never touches X."""
    with h5py.File(h5ad_path, "r") as f:
        if "obs" not in f:
            sys.exit(f"ERROR: {h5ad_path} has no /obs group - is it a valid .h5ad?")
        obs = read_elem(f["obs"])

    # Duplicate cell names would silently inflate the merge: a left join can match a
    # duplicated label more than once, and reindexing by a duplicated label multiplies
    # rows. Refuse rather than write an object with more cells than the input.
    dup = obs.index.duplicated()
    if dup.any():
        examples = list(map(str, obs.index[dup][:3]))
        sys.exit(f"ERROR: {os.path.basename(h5ad_path)} has {int(dup.sum())} duplicate "
                 f"obs_names (e.g. {examples}). Cell names must be unique to join on them - "
                 f"run adata.obs_names_make_unique() before merging.")
    return obs


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


def resolve(keys, target, obs_set, norm_index):
    """Map candidate keys to original obs_names (NaN where there is no match)."""
    if target == "exact":
        return keys.where(keys.isin(obs_set))
    return keys.map(norm_index)


def score_candidates(tcr, candidates, obs_set, norm_index, sample_col, barcode_col, n_clash):
    rows = []
    for name, fn, target in candidates:
        try:
            keys = fn(tcr, barcode_col, sample_col)
        except Exception as exc:                      # a candidate may not apply
            rows.append({"candidate": name, "target": target, "matched": 0,
                         "total": len(tcr), "match_rate": 0.0, "note": f"skipped: {exc}"})
            continue
        matched = int(resolve(keys, target, obs_set, norm_index).notna().sum())
        note = ""
        if target == "normalised":
            note = "library words dropped from sample names"
            if n_clash:
                note += f"; {n_clash} ambiguous GEX keys excluded"
        rows.append({"candidate": name, "target": target, "matched": matched,
                     "total": len(tcr), "match_rate": round(matched / max(len(tcr), 1), 6),
                     "note": note})
    report = pd.DataFrame(rows)
    # most matches first; on a tie prefer an exact convention over a normalised one
    report["_exact_first"] = (report["target"] != "exact").astype(int)
    report = (report.sort_values(["matched", "_exact_first"], ascending=[False, True])
                    .drop(columns="_exact_first").reset_index(drop=True))
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
    ap.add_argument("--library-tokens", default="GEX,TCR,VDJ,BCR,ADT",
                    help="comma-separated library-type words dropped from sample names before "
                         "matching (e.g. BTC-GBM-001-001-GEX == BTC-GBM-001-001-TCR). "
                         "Empty string requires identical sample names.")
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

    obs_set = set(map(str, obs.index))
    normalise = make_normaliser(args.library_tokens)
    norm_index, n_clash = ({}, 0) if normalise is None else build_normalised_index(obs.index, normalise)
    if normalise:
        print(f"[join] library words ignored in sample names: {args.library_tokens}"
              + (f" ({n_clash} ambiguous GEX keys excluded)" if n_clash else ""))

    # ── 2. pick the barcode convention by match rate ──────────────────────────
    candidates = build_candidates(normalise)
    report = score_candidates(tcr, candidates, obs_set, norm_index,
                              args.sample_col, args.barcode_col, n_clash)
    print("\n[join] candidate barcode keys:")
    for _, r in report.iterrows():
        print(f"   {r['candidate']:<26} {r['matched']:>7} / {r['total']:<7} "
              f"({r['match_rate'] * 100:.2f}%) {r['note']}")
    report.to_csv(os.path.join(args.tables_dir, "barcode_join_report.tsv"),
                  sep="\t", index=False)

    best = report.iloc[0]
    if best["matched"] == 0 or best["match_rate"] < args.min_match_rate:
        gex_samples = sorted({n.rsplit("_", 1)[0] for n in obs_set if "_" in n})[:5]
        tcr_samples = sorted(tcr[args.sample_col].astype(str).unique())[:5]
        sys.exit(
            f"\nERROR: best candidate '{best['candidate']}' matched only "
            f"{best['matched']}/{best['total']} TCR rows "
            f"({best['match_rate'] * 100:.2f}% < {args.min_match_rate * 100:.2f}%).\n"
            f"GEX obs_names look like: {list(map(str, obs.index[:3]))}\n"
            f"TCR barcodes look like:  {list(tcr[args.barcode_col].astype(str)[:3])}\n"
            f"GEX samples: {gex_samples}\n"
            f"TCR samples: {tcr_samples}\n"
            f"If these are the same specimens under different names, the sample names differ "
            f"by more than the library words in --library-tokens ({args.library_tokens!r}). "
            f"If they are different specimens, pick the GEX object that holds these samples.\n"
            f"See {args.tables_dir}/barcode_join_report.tsv for every candidate.")

    fn, target = next((f, t) for n, f, t in candidates if n == best["candidate"])
    tcr = tcr.copy()
    tcr["__key__"] = fn(tcr, args.barcode_col, args.sample_col)
    tcr["__obs__"] = resolve(tcr["__key__"], target, obs_set, norm_index)
    print(f"\n[join] using '{best['candidate']}' "
          f"({best['matched']}/{best['total']} TCR rows match a GEX cell)")

    # ── 3. join (one TCR row per GEX cell; duplicates reported, first kept) ───
    matched_rows = tcr[tcr["__obs__"].notna()]
    dup = int(matched_rows["__obs__"].duplicated().sum())
    if dup:
        print(f"[warn] {dup} TCR rows map to a GEX cell already matched; keeping the first of each")
    matched_rows = matched_rows.drop_duplicates(subset="__obs__", keep="first")

    drop = {"__key__", "__obs__", "__cell_id__"}
    payload = matched_rows.set_index("__obs__")[[c for c in tcr.columns if c not in drop]]
    payload = payload.add_prefix(args.prefix)

    new_obs = obs.join(payload, how="left")
    new_obs = new_obs.loc[obs.index]                       # preserve original order
    # Checked AFTER the reindex: that is where a duplicated label would multiply rows,
    # so a check before it would pass while the object silently grew.
    if len(new_obs) != len(obs):
        sys.exit(f"ERROR: the join changed the cell count ({len(obs)} -> {len(new_obs)}); "
                 f"refusing to write a mismatched object.")

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

    unmatched = tcr[tcr["__obs__"].isna()]
    unmatched.drop(columns=["__key__", "__obs__", "__cell_id__"], errors="ignore").to_csv(
        os.path.join(td, "tcr_rows_without_gex_cell.tsv"), sep="\t", index=False)

    if args.sample_col in tcr.columns:
        # Which GEX sample did each TCR sample land in? Makes a wrong pairing visible.
        gex_sample = tcr["__obs__"].map(lambda n: n.rsplit("_", 1)[0] if isinstance(n, str) and "_" in n else None)
        per_sample = (tcr.assign(matched=tcr["__obs__"].notna(), gex_sample=gex_sample)
                        .groupby(args.sample_col)
                        .agg(tcr_cells=("matched", "size"),
                             matched_to_gex=("matched", "sum"),
                             paired_gex_sample=("gex_sample",
                                                lambda s: ";".join(sorted({x for x in s if x}))))
                        .reset_index())
        per_sample["match_rate"] = (per_sample["matched_to_gex"] /
                                    per_sample["tcr_cells"].clip(lower=1)).round(6)
        per_sample.to_csv(os.path.join(td, "per_sample_join_summary.tsv"),
                          sep="\t", index=False)
        print("\n[summary] per sample (TCR sample -> GEX sample it paired with):")
        for _, r in per_sample.iterrows():
            print(f"   {str(r[args.sample_col]):<34} {int(r['matched_to_gex']):>6} / "
                  f"{int(r['tcr_cells']):<6} ({r['match_rate'] * 100:5.1f}%)  -> "
                  f"{r['paired_gex_sample'] or '(no GEX partner)'}")

    pd.DataFrame([{
        "gex_h5ad": os.path.basename(args.gex_h5ad),
        "gex_cells": len(obs),
        "tcr_tables": ";".join(os.path.basename(p) for p in paths),
        "tcr_rows": len(tcr),
        "barcode_key": best["candidate"],
        "sample_names_normalised": target == "normalised",
        "library_tokens": args.library_tokens,
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
