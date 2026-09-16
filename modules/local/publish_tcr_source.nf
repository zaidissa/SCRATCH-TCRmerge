/*
 * Republish the parts of TCRtoolkit's bridge/merged_vdj_object/ bundle that are not
 * joined onto the GEX object, so the merged dataset is self-contained:
 *
 *   *_summary.tsv        per-SAMPLE counts - cannot be joined to cells by barcode
 *   *_seurat.rds         Seurat object, for TCR-side work in R
 *   *_combineTCR.rds     scRepertoire object, likewise
 *   pre_qc_cells.tsv     per-cell but the same barcodes as post_qc_cells.tsv, so
 *                        joining it too would double every cell
 *
 * Files are staged into tcr_source/ and published verbatim - nothing is parsed here.
 */

process PUBLISH_TCR_SOURCE {

    tag "${tcr_files.size()} file(s)"
    label 'process_low'
    container "${params.container}"

    publishDir "${params.outdir}", mode: 'copy', overwrite: true

    input:
      path tcr_files

    output:
      path "tcr_source/*", emit: files

    script:
    """
    mkdir -p tcr_source
    for f in ${tcr_files.join(' ')}; do
        cp -L "\$f" "tcr_source/\$(basename \$f)"
    done
    echo "carried over \$(ls tcr_source | wc -l) file(s) into tcr_source/"
    ls -la tcr_source
    """
}
