#' Read and merge vcf containing all and only (Shapeit4) phased variants
#'
#' @param df_initial Dataframe originated by conshap() or snps_read_merge()/_haptreex()
#' @param df_phased Result of monoallelic() function
#'
#' @return A dataframe.
#' @export
#'
prep_outcome<-function(df_initial,df_phased){
  phasing_out<-unique(df_phased[,c("contig","position","Ref_code")])
  phasing_out$Alt_code<-ifelse(phasing_out$Ref_code==0,1,0)

  phasing_out$Ref_code<-as.factor(phasing_out$Ref_code)
  phasing_out$Alt_code<-as.factor(phasing_out$Alt_code)
  phasing_out$contig<-as.factor(phasing_out$contig)
  df_initial$Ref<-as.factor(df_initial$Ref)
  df_initial$Alt<-as.factor(df_initial$Alt)
  df_initial$CHROM<-as.factor(df_initial$CHROM)

  df_final<-merge(phasing_out,df_initial[,c("CHROM","POS","ID","REF","ALT")],by.x=c("contig","position"),by.y=c("CHROM","POS"))

  df_final<-merge(df_final,df_initial[!is.na(df_initial$Ref),],by.x=c("contig","position","ID","REF","ALT","Ref_code","Alt_code"),
                  by.y=c("CHROM","POS","ID","REF","ALT","Ref","Alt"),all=TRUE)
  if(nrow(unique(df_final[duplicated(df_final[,c("contig","position")])]))>0){
    stop("Duplicated haplotypes individuated")
  }
  return(df_final)
}
