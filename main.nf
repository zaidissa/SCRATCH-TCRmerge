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

// Cirro sends real JSON booleans, but `--flag true` on the command line arrives as the
// STRING "true" - which `== true` rejects, while `?:` would accept even the string
// "false". Compare the lower-cased text so both routes agree.
def truthy(value) {
    value?.toString()?.toLowerCase() in ['true', '1', 'yes', 'y', 'on']
}

// Optional inputs: resolve to the files that actually exist, dropping the rest. Declared
// as a function rather than `def x = { ... }` - Nextflow 26.04 will not let a closure
// variable be called like a function from the workflow body.
def existing_files(value) {
    if (!value) return []
    as_file_list(value).collect { file(it.trim(), checkIfExists: false) }
                       .findAll { it.exists() }
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

    // Optional inputs are staged as a placeholder when absent, so the process always has a
    // file to stage and the script decides whether it is real. The placeholders must have
    // DISTINCT names: Nextflow refuses to stage two input files with the same name, so a
    // single shared NO_FILE fails with "input file name collision".
    def no_contigs = file("${projectDir}/assets/NO_CONTIGS",    checkIfExists: true)
    def no_qc      = file("${projectDir}/assets/NO_CONTIGS_QC", checkIfExists: true)
    def no_airr    = file("${projectDir}/assets/NO_AIRR",       checkIfExists: true)
    def contig_f   = existing_files(params.contigs)
    def qc_f       = existing_files(params.contigs_passed_qc)
    def airr_f     = existing_files(params.airr_files)

    ch_contigs  = Channel.value( contig_f ? contig_f[0] : no_contigs )
    ch_qcpassed = Channel.value( qc_f     ? qc_f[0]     : no_qc )
    ch_airr     = Channel.value( airr_f ?: [no_airr] )

    MERGE_TCR_GEX( ch_gex, ch_tcr, ch_contigs, ch_qcpassed, ch_airr )

    if (truthy(params.emit_rds)) {
        TABLES_TO_RDS( MERGE_TCR_GEX.out.tables )
    }

    // Carry the rest of the TCR bundle (per-sample summaries, the Seurat/combineTCR
    // .rds objects, pre_qc_cells.tsv) into the merged dataset without joining it.
    if (params.tcr_passthrough) {
        // These paths are partly DERIVED from the joined table's folder (Cirro's file
        // listing is incomplete), so some may not exist. Skip those with a warning
        // rather than failing a merge that is otherwise fine.
        ch_passthrough = Channel
            .fromList( as_file_list(params.tcr_passthrough) )
            .flatMap { entry -> file(entry.trim(), checkIfExists: false) }
            .filter { f ->
                def present = f.exists()
                if (!present) log.warn "tcr_source: skipping missing file ${f}"
                present
            }
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
