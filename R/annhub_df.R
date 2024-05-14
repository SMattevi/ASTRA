#' Extract annotation db from AnnotationHub
#'
#' @param edb_id ensembl database id from AnnotationHub i.e. "AH113665"
#' @param window Window size up and downstream gene size.
#'
#' @import AnnotationHub
#'
#' @return A dataframe.
#' @export
#'
annotationhub_df<- function(edb_id="AH113665",window=2.5){
  ah <- AnnotationHub()
  edb <- ah[[edb_id]]
  BM <- as.data.frame(genes(edb))
  BM<-BM[,c("seqnames","start","end","gene_id","gene_biotype","symbol","description")]
  colnames(BM)<-c("chromosome_name", "start_position", "end_position",
                  "ensembl_gene_id","gene_biotype","hgnc_symbol",
                  "description")
  BM$length<-BM$end_position-BM$start_position
  BM$start_position_l<-BM$start_position-(BM$length/window)
  BM$end_position_l<-BM$end_position+(BM$length/window)

  return(BM)
}
