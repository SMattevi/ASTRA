#' Read and merge vcf containing all and only (Shapeit4) phased variants
#'
#' @param dataset_choice A string with mart dataset.
#' @param window Window size up and downstream gene size.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' biomart_df("hsapiens_gene_ensembl")
biomart_df<- function(dataset_choice="hsapiens_gene_ensembl",window=2.5){
  mart = biomaRt::useMart("ensembl", dataset=dataset_choice)
  BM<-biomaRt::getBM(attributes=c("chromosome_name", "start_position",
                                  "end_position", "ensembl_gene_id",
                                  "gene_biotype","hgnc_symbol","description"),
                     mart = mart,useCache = F)
  BM$length<-BM$end_position-BM$start_position
  BM$start_position_l<-BM$start_position-(BM$length/window)
  BM$end_position_l<-BM$end_position+(BM$length/window)

  return(BM)
}
