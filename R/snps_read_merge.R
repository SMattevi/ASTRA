#' Read and merge vcf containing all and only (Shapeit4) phased variants
#'
#' @param variant_calling_out_path A character vector with vcf path.
#' @param phased_out_path A character vector with phased vcf path.
#' @param sample_name A character vector with sample identifier in vcf.
#'
#' @import biomaRt data.table forcats tidyr
#'
#' @return A dataframe.
#' @export
#'
read_shapeit<- function(variant_calling_out_path,phased_out_path,sample_name="SAMPLE1"){

  all_snps<- data.table::fread(variant_calling_out_path , skip = "CHROM")
  prephased_snps<- data.table::fread(phased_out_path , skip = "CHROM")
  suppressWarnings(all_snps<-all_snps%>%tidyr::separate(sample_name , into=c("Ref","Alt")))
  suppressWarnings(prephased_snps<-prephased_snps%>%tidyr::separate(sample_name , into=c("Ref","Alt")))
  all_snps<-dplyr::select(all_snps,c("#CHROM","POS","ID","REF","ALT"))
  prephased_snps$CHROM<-as.factor(prephased_snps$`#CHROM`)
  all_snps$CHROM<-as.factor(all_snps$`#CHROM`)
  prephased_snps<-prephased_snps[,c("CHROM","POS","ID","REF","ALT","Ref","Alt")]
  df<- merge(prephased_snps,all_snps,by=c("CHROM","POS","ID","REF","ALT"),all=TRUE)
  return(df[,c("CHROM","POS","ID","REF","ALT","Ref","Alt")])

}
#' Read and merge vcf containing all and only (Haptreex) phased variants
#'
#' @param variant_calling_out_path A character vector with vcf path.
#' @param phased_out_path A character vector with phased txt path.
#' @param sample_name A character vector with sample identifier in vcf.
#'
#' @import biomaRt data.table forcats tidyr
#'
#' @return A dataframe.
#' @export
#'
read_haptreex<- function(variant_calling_out_path,phased_out_path,sample_name="SAMPLE1"){

  all_snps<- data.table::fread(variant_calling_out_path , skip = "CHROM")
  prephased_snps<- data.table::fread(cmd=paste0("grep -v '*' ", phased_out_path,"| grep -v 'BLOCK'"), header=FALSE, sep="\t",col.names = c("row","Ref","Alt","CHROM","POS","NA"))
  suppressWarnings(all_snps<-all_snps%>%tidyr::separate(sample_name , into=c("Ref","Alt")))

  all_snps<-dplyr::select(all_snps,c("#CHROM","POS","ID","REF","ALT"))
  all_snps$CHROM<-as.factor(all_snps$`#CHROM`)
  prephased_snps$CHROM<-as.factor(prephased_snps$CHROM)
  prephased_snps<-prephased_snps[,c("CHROM","POS","Ref","Alt")]
  df<- merge(prephased_snps,all_snps,by=c("CHROM","POS"),all=TRUE)
  return(df[,c("CHROM","POS","ID","REF","ALT","Ref","Alt")])

}
