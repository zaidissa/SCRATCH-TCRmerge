process TABLES_TO_RDS {

    tag "tables -> rds"
    label 'process_low'
    container "${params.container}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path tables

    output:
      path "tables_rds/*.rds", emit: rds

    script:
    """
    #!/usr/bin/env Rscript
    dir.create("tables_rds", showWarnings = FALSE)
    tsvs <- list.files("${tables}", pattern = "\\\\.tsv\$", full.names = TRUE)
    if (length(tsvs) == 0) stop("No .tsv files found in ${tables}")
    for (f in tsvs) {
        df  <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
        out <- file.path("tables_rds", sub("\\\\.tsv\$", ".rds", basename(f)))
        saveRDS(df, out)
        cat(sprintf("%-40s %d rows x %d cols -> %s\\n",
                    basename(f), nrow(df), ncol(df), out))
    }
    """
}
