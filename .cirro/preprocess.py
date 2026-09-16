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

from cirro.helpers.preprocess_dataset import PreprocessDataset


def pick_h5ad(ds, files):
    hits = [f for f in files if f.lower().endswith(".h5ad")]
    # Ignore anything this tool itself produced, so re-running on an output
    # dataset doesn't pick up the previous merge.
    hits = [f for f in hits if not f.rsplit("/", 1)[-1].startswith("merged")] or hits

    if not hits:
        raise ValueError(
            "No .h5ad found in the selected datasets. Select the dataset holding "
            "the annotated GEX object."
        )
    if len(hits) > 1:
        ds.logger.warning(
            f"{len(hits)} .h5ad files found; using the first: {hits[0]}. "
            "Select a single GEX dataset to remove the ambiguity."
        )
    return hits[0]


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

    try:
        pattern = dict(ds.params).get("tcr_table_pattern") or "post_qc_cells.tsv"
    except Exception:
        pattern = "post_qc_cells.tsv"

    gex = pick_h5ad(ds, files)
    ds.logger.info(f"GEX object: {gex}")
    ds.add_param("gex_h5ad", gex, overwrite=True)

    tcr = pick_tcr_tables(ds, files, pattern)
    for t in tcr:
        ds.logger.info(f"TCR table: {t}")
    ds.add_param("tcr_tables", ",".join(tcr), overwrite=True)

    ds.logger.info(
        f"Merging {len(tcr)} TCR table(s) onto 1 GEX object. The expression matrix "
        "is not read; only .obs is rewritten."
    )
    ds.logger.info(ds.params)


if __name__ == "__main__":
    main()
