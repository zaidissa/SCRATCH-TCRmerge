#!/usr/bin/env python3
"""
Cirro preprocessing for SCRATCH-TCRmerge.

Finds the two inputs among the selected datasets' files and exposes them as
Nextflow params:

    gex_h5ad    the single .h5ad to annotate
    tcr_tables  comma-separated per-cell TCR tables (TCRtoolkit's
                bridge/merged_vdj_object/post_qc_cells.tsv by default)

Both are set with overwrite=True: process-input.json declares them (as empty
strings) so they already exist in ds.params before this runs, and add_param
asserts "already exists" without it.
"""

from fnmatch import fnmatch

from cirro.helpers.preprocess_dataset import PreprocessDataset


def pick_h5ad(ds, files, pattern):
    """
    SCRATCH-QC names its object per run (e.g. GBM_DFCI1_CSF_singlet.h5ad), so match on a
    glob pattern rather than a fixed name. Anything this tool produced itself is skipped,
    so re-running on an output dataset does not pick up the previous merge.
    """
    h5ads = [f for f in files if f.lower().endswith(".h5ad")]
    h5ads = [f for f in h5ads if not f.rsplit("/", 1)[-1].startswith("merged")] or h5ads

    if not h5ads:
        raise ValueError(
            "No .h5ad found in the selected datasets. Select the dataset holding the GEX "
            "object (e.g. SCRATCH-QC's *_singlet.h5ad)."
        )

    hits = [f for f in h5ads if fnmatch(f.rsplit("/", 1)[-1], pattern)]
    if not hits:
        ds.logger.warning(
            f"No .h5ad matched '{pattern}'; falling back to any .h5ad in the dataset. "
            f"Candidates: {[f.rsplit('/', 1)[-1] for f in h5ads]}"
        )
        hits = h5ads

    if len(hits) > 1:
        ds.logger.warning(
            f"{len(hits)} .h5ad files matched '{pattern}'; using the first: {hits[0]}. "
            f"Narrow 'GEX object filename' to disambiguate. "
            f"Matches: {[f.rsplit('/', 1)[-1] for f in hits]}"
        )
    return hits[0]


def pick_passthrough(ds, files, patterns, already_used):
    """
    Everything else in TCRtoolkit's merged_vdj_object bundle: per-sample summaries and
    the Seurat / combineTCR .rds objects, plus pre_qc_cells.tsv. None of it can be joined
    onto .obs - the summaries are per-sample, the .rds files are R objects, and
    pre_qc_cells.tsv repeats the same barcodes - so it is republished verbatim instead.

    Two ways of finding them, because ds.files is Cirro's *indexed* listing and shows only
    a fraction of a dataset (6 of 257 files in one real run), which left this empty:

      1. entries matched against whatever ds.files does expose;
      2. plain filenames (no wildcard) derived as siblings of the table being joined -
         they live in the same merged_vdj_object/ folder, so the path is predictable
         without listing anything.

    A derived path is a prediction, not an observation: main.nf skips any that turns out
    not to exist rather than failing the run.
    """
    globs = [p.strip() for p in patterns.split(",") if p.strip()]
    used = set(already_used)
    parents = {u.rsplit("/", 1)[0] for u in used}
    found = {}                                   # path -> how we found it

    # 1. whatever the dataset listing exposes
    listed = [
        f for f in files
        if f not in used and any(fnmatch(f.rsplit("/", 1)[-1], g) for g in globs)
    ]
    # Only carry files beside a table we are actually joining, so selecting a big
    # upstream dataset does not drag unrelated .rds files into the output. This filter
    # is strict: if nothing sits beside the table, the answer is nothing, not "fall back
    # to every match anywhere in the dataset".
    if parents:
        listed = [f for f in listed if f.rsplit("/", 1)[0] in parents]
    for f in listed:
        found[f] = "listed"

    # 2. siblings of the joined table(s), for entries that name a file outright
    for folder in sorted(parents):
        for g in globs:
            if any(ch in g for ch in "*?["):      # a glob cannot be derived, only matched
                continue
            candidate = f"{folder}/{g}"
            if candidate not in used and candidate not in found:
                found[candidate] = "derived"

    for path in sorted(found):
        ds.logger.info(f"TCR file (carried over, {found[path]}): {path}")
    return sorted(found)


def pick_tcr_tables(ds, files, pattern):
    exact = [f for f in files if f.rsplit("/", 1)[-1] == pattern]
    if exact:
        return exact

    ds.logger.warning(
        f"No file named '{pattern}' found; falling back to any *_cells.tsv."
    )
    cells = [f for f in files if f.rsplit("/", 1)[-1].endswith("_cells.tsv")]
    if cells:
        # pre_qc and post_qc both match; prefer post_qc (QC-filtered contigs).
        post = [f for f in cells if "post_qc" in f.rsplit("/", 1)[-1]]
        return post or cells

    raise ValueError(
        "No per-cell TCR table found. Expected TCRtoolkit's "
        "bridge/merged_vdj_object/post_qc_cells.tsv in one of the selected datasets."
    )


def derive_contig_tables(ds, tcr_paths, names):
    """
    TCRtoolkit writes its per-contig tables to VDJ_QC/VDJ_QC/tables/ in the SAME dataset
    as the per-cell table, at a fixed depth: .../data/bridge/merged_vdj_object/<table>.tsv
    and .../data/VDJ_QC/VDJ_QC/tables/<contigs>.tsv. Derive them from the picked table
    rather than scanning ds.files, which only lists part of a dataset.
    """
    out = {}
    if not tcr_paths:
        return out
    folder = tcr_paths[0].rsplit("/", 1)[0]              # .../data/bridge/merged_vdj_object
    data_root = folder.rsplit("/", 2)[0]                 # .../data
    for key, name in names.items():
        out[key] = f"{data_root}/VDJ_QC/VDJ_QC/tables/{name}"
        ds.logger.info(f"contig table ({key}, derived): {out[key]}")
    return out


def main():
    ds = PreprocessDataset.from_running()
    ds.logger.info("List of starting params")
    ds.logger.info(ds.params)

    files = list(ds.files["file"]) if len(ds.files) else []
    ds.logger.info(f"{len(files)} files visible across the selected datasets")

    # ds.params may be a plain dict or a params object depending on cirro version.
    try:
        params = dict(ds.params)
    except Exception:
        params = {}
    pattern = params.get("tcr_table_pattern") or "post_qc_cells.tsv"
    gex_pattern = params.get("gex_h5ad_pattern") or "*_singlet.h5ad"

    # A file picked explicitly in the form wins; the pattern scan is the fallback, so
    # selecting datasets and letting this discover the files still works unchanged.
    gex_picked = str(params.get("gex_h5ad_selected") or "").strip()
    if gex_picked:
        gex = gex_picked
        ds.logger.info(f"GEX object (picked in the form): {gex}")
    else:
        gex = pick_h5ad(ds, files, gex_pattern)
        ds.logger.info(f"GEX object (auto-discovered via '{gex_pattern}'): {gex}")
    ds.add_param("gex_h5ad", gex, overwrite=True)

    tcr_picked = str(params.get("tcr_table_selected") or "").strip()
    if tcr_picked:
        tcr = [tcr_picked]
        ds.logger.info(f"TCR table (picked in the form, joined): {tcr_picked}")
    else:
        tcr = pick_tcr_tables(ds, files, pattern)
        for t in tcr:
            ds.logger.info(f"TCR table (auto-discovered, joined): {t}")
    ds.add_param("tcr_tables", ",".join(tcr), overwrite=True)

    passthrough_patterns = params.get("tcr_passthrough_patterns") or (
        "pre_qc_cells.tsv,pre_qc_summary.tsv,post_qc_summary.tsv,"
        "pre_qc_seurat.rds,post_qc_seurat.rds,"
        "pre_qc_combineTCR.rds,post_qc_combineTCR.rds"
    )
    extra = pick_passthrough(ds, files, passthrough_patterns, tcr)
    ds.add_param("tcr_passthrough", ",".join(extra), overwrite=True)

    # ── per-contig side table ────────────────────────────────────────────────
    want_contigs = str(params.get("build_contig_table", True)).lower() not in ("false", "0", "no", "")
    contigs = contigs_qc = ""
    if want_contigs:
        names = {
            "contigs": params.get("contigs_filename") or "contigs_before_qc.tsv",
            "contigs_passed_qc": params.get("contigs_qc_filename") or "contigs_after_qc.tsv",
        }
        derived = derive_contig_tables(ds, tcr, names)
        contigs, contigs_qc = derived.get("contigs", ""), derived.get("contigs_passed_qc", "")
    ds.add_param("contigs", contigs, overwrite=True)
    ds.add_param("contigs_passed_qc", contigs_qc, overwrite=True)

    # Full contig sequences live in Cell Ranger's AIRR files, in the upstream align
    # dataset. Only used if such a dataset is among the selections and Cirro lists them.
    airr = sorted(f for f in files if f.rsplit("/", 1)[-1] == "airr_rearrangement.tsv")
    for a in airr:
        ds.logger.info(f"AIRR file (full contig sequences): {a}")
    if want_contigs and not airr:
        ds.logger.info(
            "No airr_rearrangement.tsv among the selected datasets: the contig table will "
            "carry the assembled V(D)J regions (~340 nt) rather than full contig sequences. "
            "Add the SCRATCH-align dataset to include them."
        )
    ds.add_param("airr_files", ",".join(airr), overwrite=True)

    ds.logger.info(
        f"Merging {len(tcr)} TCR table(s) onto 1 GEX object, carrying {len(extra)} "
        "further TCR file(s) into tcr_source/. The expression matrix is not read; "
        "only .obs is rewritten."
    )
    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()
