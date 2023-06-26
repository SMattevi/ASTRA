#' Apply manual phasing
#'
#' @param df A df containg annotated SNPs to phase (and phased)
#' @importFrom dplyr arrange
#' @importFrom stats sd
#'
#' @return A dataframe.
#'
#' @export
#'
manualP<- function(df){

  df$Ref_code<-df$Ref
  df$Ref_code<-ifelse(is.na(df$Ref_code),0,df$Ref_code)

  df$Ref<-ifelse(is.na(df$Ref),"NP","PHASED_WS")
  df_in<- manual_prep_phasing(df)

  df_in$Ref<-ifelse(df_in$Ref=="PHASED","NP","PHASED_WS")

  df_in<-df_in%>%arrange(-.data$position)

  df_out<- manual_prep_phasing(df_in)
  return(df_out)

}

manual_prep_phasing<- function(df) {

  mpp<-df

  for(i in unique(mpp$ensembl_gene_id)){

    #select one gene
    df_tmp<-mpp[mpp$ensembl_gene_id==i,]
    #if everything already phased go to the next gene
    if(all(df_tmp$Ref == "PHASED" | df_tmp$Ref == "PHASED_WS")){

      next

    }else if(length(unique(df_tmp$variantID)) == 1){

      snps_num <- length(unique(mpp$variantID[mpp$variantID==unique(df_tmp$variantID)]))
      if(snps_num == 1){

        mpp$Ref[mpp$ensembl_gene_id==i] <- "PHASED"
        df_tmp$Ref[df_tmp$ensembl_gene_id==i] <- "PHASED"
        next

      }
    }

    sd_tmp<-(stats::sd(df_tmp$refFrac))
    # #add to visited all SNPs already phased (by other phaser or previously in this step)
    visited<-c()
    visited<-append(visited,df_tmp[df_tmp$Ref=="PHASED_WS" | df_tmp$Ref=="PHASED",refFrac])

    for(j in unique(df_tmp$variantID[df_tmp$Ref=="NP"])){

      sd_tmp<-ifelse(is.na(stats::sd(visited)),sd_tmp,(stats::sd(c(visited,df_tmp$refFrac[df_tmp$variantID==j]))))

      copy_tmp<-df_tmp

      tmp<-df_tmp$refAllele[df_tmp$variantID==j]
      df_tmp$refAllele[df_tmp$variantID==j]<-df_tmp$altAllele[df_tmp$variantID==j]
      df_tmp$altAllele[df_tmp$variantID==j]<-tmp
      tmp<-df_tmp$refCount[df_tmp$variantID==j]
      df_tmp$refCount[df_tmp$variantID==j]<-df_tmp$altCount[df_tmp$variantID==j]
      df_tmp$altCount[df_tmp$variantID==j]<-tmp
      df_tmp$refFrac<-df_tmp$refCount/df_tmp$totalCount

      sd_v1<-(stats::sd(c(visited,df_tmp$refFrac[df_tmp$variantID==j])))
      sd_vis<-ifelse(is.na(sd_v1),(stats::sd(df_tmp$refFrac)),sd_v1)

      if(sd_tmp<sd_vis){ #(stats::sd(df_tmp$refFrac))){

        df_tmp<-copy_tmp
        mpp$Ref[mpp$variantID==j]<-"PHASED"
        df_tmp$Ref[df_tmp$variantID==j]<-"PHASED"

      }
      else{

        tmp<-mpp$refAllele[mpp$variantID==j]
        mpp$refAllele[mpp$variantID==j]<-mpp$altAllele[mpp$variantID==j]
        mpp$altAllele[mpp$variantID==j]<-tmp
        tmp<-mpp$refCount[mpp$variantID==j]
        mpp$refCount[mpp$variantID==j]<-mpp$altCount[mpp$variantID==j]
        mpp$altCount[mpp$variantID==j]<-tmp
        df_tmp$Ref[df_tmp$variantID==j]<-"PHASED"
        mpp$Ref[mpp$variantID==j]<-"PHASED"
        mpp$refFrac<-mpp$refCount/mpp$totalCount
        sd_tmp<-(stats::sd(df_tmp$refFrac))

        mpp$Ref_code[mpp$variantID==j]<-ifelse(mpp$Ref_code[mpp$variantID==j]==0,1,0)

      }
      visited<-append(visited,df_tmp[df_tmp$variantID==j,refFrac])
    }
  }

  return(mpp)
}
