#' Read and hapotype phasing results from Shapeit, Haptreex,...
#'
#' @param df1 First dataframe with haplotype calling results
#' @param df2 Second dataframe with haplotype calling results
#' @param BMsel BM dataframe containing only chromosome under analysis
#' @param CHROM Chromosome column
#' @param POS Position column
#' @param Ref Reference haplotype column
#' @param REF Reference allele column
#' @param ALT Alternative allele column
#' @param DF1 Name of the tool used for the haplotype calling in df1
#' @param DF2 Name of the tool used for the haplotype calling in df2
#'
#' @return A dataframe.
#' @export
#'
conshap<- function(df1,df2,BMsel,CHROM="CHROM",POS="POS",Ref="Ref",REF="REF",ALT="ALT",DF1="DF1",DF2="DF2"){

  df1$Ref<-df1[[Ref]]
  df2$Ref<-df2[[Ref]]

  df_merged<-merge(df1,df2,by=c("CHROM","POS","ID","REF","ALT"))

  #Gene annotation and extraction of SNPs overlapping a gene
  df_merged$variantID<-paste0(df_merged[[CHROM]],":",df_merged[[POS]],":",df_merged[[REF]],":",df_merged[[ALT]])

  isnps <- with(df_merged, IRanges::IRanges(POS, width=1, names=variantID))
  igenes <- with(BMsel, IRanges::IRanges(start_position_l, end_position_l, names=ensembl_gene_id))
  olaps <- GenomicRanges::findOverlaps(isnps, igenes)

  overlapping<-cbind(df_merged[IRanges::from(olaps),], BMsel[IRanges::to(olaps),])

  #in Ref column will be saved final haplotype
  overlapping$Ref<-as.integer(NA)
  df_merged$Ref<-as.integer(NA)
  df_merged$Tool<-NA

  overlapping$Ref.x<-as.integer(overlapping$Ref.x)
  overlapping$Ref.y<-as.integer(overlapping$Ref.y)

  overlapping$Alt.x<-as.integer(overlapping$Alt.x)
  overlapping$Alt.y<-as.integer(overlapping$Alt.y)

  n=length(unique(overlapping$ensembl_gene_id))
  ii=1

  for(gene in unique(overlapping$ensembl_gene_id)){

    cat(paste0(round(ii / n * 100), '% completed'))

    #look at one gene at the time
    tmp_df<-overlapping[overlapping$ensembl_gene_id==gene,]

    l<-length(tmp_df$variantID)

    tmp_df$Ref.x<-as.integer(tmp_df$Ref.x)
    tmp_df$Ref.y<-as.integer(tmp_df$Ref.y)

    for(i in 1:l){
      #for each snp if
      #1. not already defined haplotype
      if(is.na(tmp_df$Ref[i])){
        #1a. df1 has haplotype definition
        if(!is.na(tmp_df$Ref.x[i])){
          #1a_bis. df2 has NOT haplotype
          if(is.na(tmp_df$Ref.y[i])){
            #then Ref is df1 haplotype
            tmp_df$Ref[i]<-tmp_df$Ref.x[i]
            df_merged$Tool[df_merged$variantID==tmp_df$variantID[i]]<-DF1
          }
          #1a_bis. df2 has haplotype = to df1
          else if(tmp_df$Ref.x[i]==tmp_df$Ref.y[i]){
            #then Ref is equal to both the haplotypes
            tmp_df$Ref[i]<-tmp_df$Ref.x[i]
            df_merged$Tool[df_merged$variantID==tmp_df$variantID[i]]<-paste0(DF1,"+",DF2)
          }
          #1a_bis. df2 has haplotype != from df1
          else if(tmp_df$Ref.x[i]!=tmp_df$Ref.y[i]){
            #assign the haplotype of the df with minimum number of NAs, swap remaining haplotypes for the other df
            if(sum(is.na(tmp_df$Ref.x))>=sum(is.na(tmp_df$Ref.y))){
              tmp_df$Ref.x[i:l]<-1-tmp_df$Ref.x[i:l]
              tmp_df$Ref[i]<-tmp_df$Ref.y[i]
              df_merged$Tool[df_merged$variantID==tmp_df$variantID[i]]<-DF2

              overlapping$Ref.x[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]<-1-overlapping$Ref.x[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]
              overlapping$Alt.x[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]<-1-overlapping$Alt.x[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]
            }
            else{
              tmp_df$Ref.y[i:l]<-1-tmp_df$Ref.y[i:l]
              tmp_df$Ref[i]<-tmp_df$Ref.x[i]
              df_merged$Tool[df_merged$variantID==tmp_df$variantID[i]]<-DF1

              overlapping$Ref.y[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]<-1-overlapping$Ref.y[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]
              overlapping$Alt.y[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]<-1-overlapping$Alt.y[which(overlapping$variantID==tmp_df$variantID[i])[1]:nrow(overlapping)]
            }
          }
        }
        #1b. if df1 has NOT haplotype definition
        else{
          #1b_bis. if df2 has haplotype definition use this as haplotype
          if(!is.na(tmp_df$Ref.y[i])){
            tmp_df$Ref[i]<-tmp_df$Ref.y[i]
            df_merged$Tool[df_merged$variantID==tmp_df$variantID[i]]<-DF2
          }
        }
        #2. assign results to overlapping and final dataframes
        overlapping[overlapping$variantID==tmp_df$variantID[i],Ref]<-tmp_df$Ref[i]
        df_merged[df_merged$variantID==tmp_df$variantID[i],Ref]<-tmp_df$Ref[i]
      }
    }
    if (ii == n) cat(': Done')
    else cat('\014')

    ii=ii+1

  }
  df_merged$Alt<-1-df_merged$Ref
  return(df_merged[,c("CHROM","POS","ID","REF","ALT","Ref","Alt","Tool")])
}
