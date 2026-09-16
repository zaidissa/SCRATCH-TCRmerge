#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * SCRATCH-TCRmerge
 *
 * Joins per-cell TCR tables (TCRtoolkit's bridge/merged_vdj_object/*_cells.tsv)
 * onto a GEX AnnData by cell barcode, and emits:
 *   - merged.h5ad          GEX object with TCR columns added to .obs
 *   - tables/*.tsv         the same joined metadata plus join diagnostics
 *   - tables/*.rds         optional R-side copies of those tables
 *
 * The expression matrix is never loaded: the .h5ad is copied and only its /obs
 * group is rewritten, so runtime and memory do not scale with cell count.
 */

include { MERGE_TCR_GEX } from './modules/local/merge_tcr_gex.nf'
include { TABLES_TO_RDS } from './modules/local/tables_to_rds.nf'

workflow {

    if (!params.gex_h5ad) {
        error "No GEX object given. Set --gex_h5ad /path/to/object.h5ad"
    }
    if (!params.tcr_tables) {
        error "No TCR tables given. Set --tcr_tables 'a.tsv,b.tsv' (or a glob)."
    }

    ch_gex = Channel.fromPath(params.gex_h5ad, checkIfExists: true)

    // Accept either a comma-separated list (what the Cirro preprocess writes)
    // or a glob. split(',', -1) keeps trailing empties out of the file list.
    def tcr_entries = params.tcr_tables.toString().contains(',')
        ? params.tcr_tables.toString().split(',', -1).findAll { it?.trim() }
        : [ params.tcr_tables.toString() ]

    ch_tcr = Channel
        .fromList(tcr_entries)
        .flatMap { entry -> file(entry.trim(), checkIfExists: true) }
        .collect()

    MERGE_TCR_GEX( ch_gex, ch_tcr )

    if (params.emit_rds == true) {
        TABLES_TO_RDS( MERGE_TCR_GEX.out.tables )
    }

    // Nextflow 26.04's parser rejects a top-level `workflow.onComplete` block
    // ("statements cannot be mixed with script declarations"), so the completion
    // handler is registered from inside the workflow body instead.
    workflow.onComplete = {
        log.info(workflow.success
            ? "SCRATCH-TCRmerge finished. Results: ${params.outdir}"
            : "SCRATCH-TCRmerge failed - see the error above.")
    }
}
