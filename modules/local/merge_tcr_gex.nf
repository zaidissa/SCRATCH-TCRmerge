process MERGE_TCR_GEX {

    tag "${gex_h5ad.simpleName}"
    label 'process_medium'
    container "${params.container}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path gex_h5ad
      path tcr_tables
      path contigs           // assets/NO_FILE when not supplied
      path contigs_passed_qc // assets/NO_FILE when not supplied
      path airr_files        // assets/NO_FILE when not supplied

    output:
      path "merged.h5ad",            emit: merged_h5ad
      path "merged.tcr_only.h5ad",   emit: tcr_only_h5ad, optional: true
      path "tables",                 emit: tables
      path "merge_log.txt",          emit: log

    script:
    def tables_arg = tcr_tables.collect { "--tcr-table ${it}" }.join(' ')
    // Cirro sends a real JSON boolean, but `--subset_to_tcr true` on the command line
    // arrives as the STRING "true": `== true` is false for that, and `?:` is true for
    // the string "false". Compare the lower-cased text instead so both agree.
    def subset_arg = (params.subset_to_tcr?.toString()?.toLowerCase() in
                      ['true', '1', 'yes', 'y', 'on']) ? '--subset-to-tcr' : ''
    // Optional inputs arrive as a NO_* placeholder when absent - one distinct name each,
    // because Nextflow refuses to stage two input files with the same name. Tested inline:
    // Nextflow 26.04 rejects calling a named closure variable like a function (a closure
    // passed to findAll is fine - it is the `name(args)` call on a variable that is not).
    def contig_arg = (contigs && !contigs.name.startsWith('NO_')) ? "--contigs ${contigs}" : ''
    def qc_arg     = (contigs_passed_qc && !contigs_passed_qc.name.startsWith('NO_'))
                     ? "--contigs-passed-qc ${contigs_passed_qc}" : ''
    def airr_list  = (airr_files instanceof List ? airr_files : [airr_files])
                     .findAll { it && !it.name.startsWith('NO_') }
    def airr_arg   = airr_list ? "--airr ${airr_list.join(',')}" : ''
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
        --library-tokens '${params.library_tokens}' \\
        ${contig_arg} ${qc_arg} ${airr_arg} \\
        --min-match-rate ${params.min_match_rate} \\
        ${subset_arg} > merge_log.txt 2>&1 || { cat merge_log.txt; exit 1; }
    cat merge_log.txt
    """
}
