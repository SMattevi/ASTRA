#' Read and merge vcf containing all and only (Shapeit4) phased variants
#'
#' @param salmon_path A character vector with quant.sf path.
#' @param txdb txdb for salmon import
#' @param counts Filtering value on salmon counts
#'
#' @import tximport
#'
#' @return A dataframe.
#' @export
#'
read_salmon<- function(salmon_path,txdb,counts){
  k <- keys(txdb, keytype = "TXID")
  tx2gene <- select(txdb, k, "GENEID", "TXID")

  txi.salmon <- tximport::tximport(salmon_path, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=T)
  salmon_counts <- as.data.frame(txi.salmon$counts)
  colnames(salmon_counts) <- "counts"

  return(salmon_counts)
}
