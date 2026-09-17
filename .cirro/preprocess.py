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
    """
    globs = [p.strip() for p in patterns.split(",") if p.strip()]
    used = set(already_used)
    hits = [
        f for f in files
        if f not in used and any(fnmatch(f.rsplit("/", 1)[-1], g) for g in globs)
    ]
    # Only carry files that sit beside a table we are actually joining, so selecting a
    # big upstream dataset does not drag unrelated .rds files into the output.
    if used:
        parents = {u.rsplit("/", 1)[0] for u in used}
        beside = [f for f in hits if f.rsplit("/", 1)[0] in parents]
        if beside:
            hits = beside
    return sorted(hits)


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
        "*_summary.tsv,*_seurat.rds,*_combineTCR.rds,pre_qc_cells.tsv"
    )
    extra = pick_passthrough(ds, files, passthrough_patterns, tcr)
    for e in extra:
        ds.logger.info(f"TCR file (carried over, not joined): {e}")
    ds.add_param("tcr_passthrough", ",".join(extra), overwrite=True)

    ds.logger.info(
        f"Merging {len(tcr)} TCR table(s) onto 1 GEX object, carrying {len(extra)} "
        "further TCR file(s) into tcr_source/. The expression matrix is not read; "
        "only .obs is rewritten."
    )
    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()
