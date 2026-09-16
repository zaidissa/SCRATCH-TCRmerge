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

include { MERGE_TCR_GEX }      from './modules/local/merge_tcr_gex.nf'
include { TABLES_TO_RDS }      from './modules/local/tables_to_rds.nf'
include { PUBLISH_TCR_SOURCE } from './modules/local/publish_tcr_source.nf'

// Comma-separated list (what Cirro's preprocess.py writes) or a single path/glob.
// split(',', -1) keeps trailing empty entries out of the file list.
//
// Declared as a function, not `def x = { ... }`: Nextflow 26.04 rejects a top-level
// closure assignment as a statement mixed with script declarations.
def as_file_list(value) {
    value.toString().contains(',')
        ? value.toString().split(',', -1).findAll { it?.trim() }
        : [ value.toString() ]
}

workflow {

    if (!params.gex_h5ad) {
        error "No GEX object given. Set --gex_h5ad /path/to/object.h5ad"
    }
    if (!params.tcr_tables) {
        error "No TCR tables given. Set --tcr_tables 'a.tsv,b.tsv' (or a glob)."
    }

    ch_gex = Channel.fromPath(params.gex_h5ad, checkIfExists: true)

    ch_tcr = Channel
        .fromList( as_file_list(params.tcr_tables) )
        .flatMap { entry -> file(entry.trim(), checkIfExists: true) }
        .collect()

    MERGE_TCR_GEX( ch_gex, ch_tcr )

    if (params.emit_rds == true) {
        TABLES_TO_RDS( MERGE_TCR_GEX.out.tables )
    }

    // Carry the rest of the TCR bundle (per-sample summaries, the Seurat/combineTCR
    // .rds objects, pre_qc_cells.tsv) into the merged dataset without joining it.
    if (params.tcr_passthrough) {
        ch_passthrough = Channel
            .fromList( as_file_list(params.tcr_passthrough) )
            .flatMap { entry -> file(entry.trim(), checkIfExists: true) }
            .collect()

        PUBLISH_TCR_SOURCE( ch_passthrough )
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
