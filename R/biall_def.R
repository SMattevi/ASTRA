#' Function for the definition of biallelic/monoallelic SNPs based on allelic counts
#'
#' @param df A dataframe with columns variantID, ensembl_gene_id, hgnc_symbol, refFrac.
#' @param refFracValue refFrac values where to cut biall/monoall default 0.15,0.85
#' @param biallObs Minimum number of SNPs supporting biallelic observation
#' @param monoallObs Minimum number of SNPs supporting monoallelic observation
#' @param bt biallelic threshold (min number of SNPs to consider a gene always biallelic)
#' @param br min number of reads on minor allele for biallelic definition
#'
#' @importFrom stats aggregate
#' @import data.table
#' @import ggplot2
#' @importFrom scales percent
#'
#' @return A df
#' @export
#'
biallmonoall<-function(df,refFracValue=c(0.15,0.85),biallObs=1,monoallObs=1,bt=2,br=2){
  #considered as biallelic when reference fraction (refCount/totCount) comprised between refFracValues
  bial<-stats::aggregate(variantID~ensembl_gene_id,df[df$refFrac>refFracValue[1]&df$refFrac<refFracValue[2]&df$refCount>=br&df$altCount>=br,],length)

  #annotate with hgnc_symbol when present
  genesymbol<-base::unique(df[,c("hgnc_symbol","ensembl_gene_id")])
  bial<-base::merge(genesymbol,bial,by=c("ensembl_gene_id"),all.y=T)

  #consider as bial only genes with at least 1 evidences
  data<-bial[bial$variantID>=biallObs,]

  base::colnames(data)[base::colnames(data) == 'variantID'] <- 'variantID.biall'


  ##considered as monoallelic when reference fraction (refCount/totCount) <= refFracValue[1] or >= refFracValue[2]
  mono<-stats::aggregate(variantID~ensembl_gene_id,df[df$refFrac<=refFracValue[1]|df$refFrac>=refFracValue[2],],length)
  mono<-base::merge(genesymbol,mono,by=c("ensembl_gene_id"),all.y=T)

  #consider as monoal only genes with at least monoallObs evidences
  data1<-mono[mono$variantID>=monoallObs,]

  base::colnames(data1)[base::colnames(data1) == 'variantID'] <- 'variantID.mono'

  data_merge<-base::merge(data,data1,by=c("ensembl_gene_id"),all=TRUE)
  data_merge[is.na(data_merge)] <- 0

  # if(dim(data_merge)[1]!=0){
  #   for(i in 1:nrow(data_merge)){
  #     if(data_merge$variantID.biall[i]>data_merge$variantID.mono[i]){
  #       data1<-data1[!(data1$ensembl_gene_id==data_merge$ensembl_gene_id[i]),]
  #       df<-df[!(df$ensembl_gene_id==data_merge$ensembl_gene_id[i]&(df$refFrac<=refFracValue[1]|df$refFrac>=refFracValue[2]))]
  #     }else if(data_merge$variantID.biall[i]<data_merge$variantID.mono[i]){
  #       data<-data[!(data$ensembl_gene_id==data_merge$ensembl_gene_id[i]),]
  #       df<-df[!(df$ensembl_gene_id==data_merge$ensembl_gene_id[i]&df$refFrac>refFracValue[1]&df$refFrac<refFracValue[2])]
  #     }else{
  #       data<-data[!(data$ensembl_gene_id==data_merge$ensembl_gene_id[i]),]
  #       data1<-data1[!(data1$ensembl_gene_id==data_merge$ensembl_gene_id[i]),]
  #       df<-df[!(df$ensembl_gene_id==data_merge$ensembl_gene_id[i])]
  #     }
  #   }
  # }

  if(dim(data_merge)[1]!=0){
    for(i in 1:nrow(data_merge)){
      if(data_merge$variantID.biall[i]>=data_merge$variantID.mono[i]|data_merge$variantID.biall[i]>=bt){
        data1<-data1[!(data1$ensembl_gene_id==data_merge$ensembl_gene_id[i]),]
        #to remove non concordant snps uncomment following line
        #df<-df[!(df$ensembl_gene_id==data_merge$ensembl_gene_id[i]&(df$refFrac<=refFracValue[1]|df$refFrac>=refFracValue[2]))]
      }else{
        data<-data[!(data$ensembl_gene_id==data_merge$ensembl_gene_id[i]),]
        #df<-df[!(df$ensembl_gene_id==data_merge$ensembl_gene_id[i]&df$refFrac>refFracValue[1]&df$refFrac<refFracValue[2])]
      }
    }
  }

  df$type=ifelse(df$ensembl_gene_id%in%data$ensembl_gene_id, "Biallelic" , ifelse(df$ensembl_gene_id%in%data1$ensembl_gene_id , "Monoallelic" , "None"))

  b <- length(unique(data$ensembl_gene_id))
  m <- length(unique(data1$ensembl_gene_id))
  r <- length(unique(genesymbol$ensembl_gene_id)) - b - m
  sum_df <- data.frame(Genes=c('Biallelic','Monoallelic','Remaning'),
                           number=c(b,m,r))
  sum_df$percentage <- sum_df$number/sum(sum_df$number)


  # creating a pie chart
  pie <- ggplot2::ggplot(sum_df, ggplot2::aes(x="", y=number, fill=Genes))+
    ggplot2::geom_bar(width = 1, stat = "identity") + ggplot2::coord_polar("y", start=0) +
    ggplot2::theme(axis.text=ggplot2::element_blank(),axis.ticks = ggplot2::element_blank(),
                   axis.title = ggplot2::element_blank()) +
    ggplot2::geom_text(ggplot2::aes(label = scales::percent(round(percentage,2))), position = ggplot2::position_stack(vjust = 0.4))

  print(pie)
  print(base::paste0("Number of biallelic genes: ",length(unique(data$ensembl_gene_id)),". Number of presumable monoallelic genes: ",length(unique(data1$ensembl_gene_id))))
  return(df)
}
