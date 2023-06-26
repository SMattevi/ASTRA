#' Flip given alleles
#'
#' @param df_toswap A dataframe
#' @param col1 Col name first allele to swap
#' @param col2 Col name second allele to swap
#' @param rownum Row num to swap
#'
#' @return A dataframe.
#'
flip_alleles<- function(df_toswap,col1, col2,rownum){
  tmp<-df_toswap[[col1]][[rownum]]
  df_toswap[[col1]][[rownum]]<-df_toswap[[col2]][[rownum]]
  df_toswap[[col2]][[rownum]]<-tmp
  return(df_toswap)
}
