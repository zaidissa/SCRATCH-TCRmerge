process MERGE_TCR_GEX {

    tag "${gex_h5ad.simpleName}"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path gex_h5ad
      path tcr_tables

    output:
      path "merged.h5ad",            emit: merged_h5ad
      path "merged.tcr_only.h5ad",   emit: tcr_only_h5ad, optional: true
      path "tables",                 emit: tables
      path "merge_log.txt",          emit: log

    script:
    def tables_arg = tcr_tables.collect { "--tcr-table ${it}" }.join(' ')
    def subset_arg = (params.subset_to_tcr == true) ? '--subset-to-tcr' : ''
    // Nextflow runs task scripts under `set -e`, so a plain `cmd > log; cat log` aborts
    // BEFORE the cat when cmd fails - the traceback lands in merge_log.txt and never
    // reaches the task output. `|| { cat ...; exit 1; }` prints it first, then fails.
    """
    merge_tcr_gex.py \\
        --gex-h5ad ${gex_h5ad} \\
        ${tables_arg} \\
        --out-h5ad merged.h5ad \\
        --tables-dir tables \\
        --sample-col '${params.sample_col}' \\
        --barcode-col '${params.barcode_col}' \\
        --cell-id-col '${params.cell_id_col}' \\
        --prefix '${params.obs_prefix}' \\
        --min-match-rate ${params.min_match_rate} \\
        ${subset_arg} > merge_log.txt 2>&1 || { cat merge_log.txt; exit 1; }
    cat merge_log.txt
    """
}
