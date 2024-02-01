#' Read and merge vcf containing all and only (Shapeit4) phased variants
#'
#' @param df A dataframe.
#' @param genes_to_plot A dataframe.
#' @param title Title for the plot
#'
#' @import ggplot2
#' @import ggpattern
#' @importFrom stats median
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#'
#' @return A ggplot
#' @export
#'
plot_ase<- function(df,genes_to_plot,title){
  medianFrac<-stats::median(df$refFrac)

  sel_genes<-melt(df[df$hgnc_symbol%in%genes_to_plot,], measure.vars = c("refCount", "altCount"))#&(df$refFrac<=(medianFrac-0.15)|df$refFrac>=(medianFrac+0.15)),], measure.vars = c("refCount", "altCount"))
  if (dim(sel_genes)[1]>0) {
    sel_genes$Ref<-ifelse(is.na(sel_genes$Ref),"NOT_PHASED",ifelse(sel_genes$Ref==0|sel_genes$Ref==1,"PHASED",sel_genes$Ref))
    sel_genes$variantID<-paste0(sel_genes$variantID,"_",sel_genes$hgnc_symbol)
    plot_res<- sel_genes %>%
      mutate(variantID = fct_reorder(.data$variantID, .data$position,.desc=TRUE)) %>% ggplot2::ggplot(ggplot2::aes(.data$value,.data$variantID,fill=.data$variable))+
      ggplot2::geom_bar(position = "fill", stat = "identity")+
      ggpattern::geom_bar_pattern(
        aes(pattern=.data$Ref),stat = "identity",position = "fill",
        color = "black",
        pattern_fill = "black",
        pattern_angle = 45,
        pattern_density = 0.1,
        pattern_spacing = 0.018,
        pattern_key_scale_factor = 0.6)+
      guides(pattern = guide_legend(override.aes = list(fill = "white"), order = 2),
             fill = guide_legend(override.aes = list(pattern = "none", order=1)))+
      ggplot2::geom_text(aes(x = 0.85, label = .data$tech),color="white")+
      ggplot2::geom_point(aes(x=1.2,size=.data$totalCount),color="black")+
      ggplot2::ggtitle(title)+
      ggplot2::scale_fill_discrete(name = "", labels = c("A1_counts", "A2_counts"))+
      ggplot2::geom_vline(xintercept = .85)+ggplot2::geom_vline(xintercept = .15)
    #  color = "black") #+#+
    #ggplot2::facet_grid(~cluster)
    return(plot_res)
  }
}

#' Plot functional annotation of ASE/ASA analysis results. Plotting function from https://r-graph-gallery.com/297-circular-barplot-with-groups.html
#'
#' @param df A dataframe.
#' @param ASC Allele Specific categorization in "biallelic/monoallelic/other" column identifier.
#' @param gene Gene-id column identifier.
#' @param annotation Functional annotation column to plot
#'
#' @import fmsb
#' @importFrom reshape2 dcast
#'
#' @return A plot
#' @export
#'
plot_annotation<- function(df,ASC,gene,annotation){
  data<-aggregate(df[[gene]] ~ as.factor(df[[ASC]]) + df[[annotation]], FUN = length)
  colnames(data)<-c("ASC","annotation","gene")
  #data = data %>% arrange(ASC, gene)
  data=reshape2::dcast(data,ASC~annotation,value.var = "gene")
  data[is.na(data)]=0
  rownames(data)<-data$ASC
  data<-data[,-1]

  data <- rbind(rep(max(data),ncol(data)) , rep(0,ncol(data)) , data)

  # Color vector
  colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
  colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

  # plot with default options:
  fmsb::radarchart( data  , axistype= 1,
              #custom polygon
              pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
              #custom the grid
              cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,max(data),max(data)/4), cglwd=0.8,
              #custom labels
              vlcex=0.8
  )

  # Add a legend
  legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)
}

#' Read and merge vcf containing all and only (Shapeit4) phased variants
#'
#' @param df A dataframe.
#' @param msigdb_species species.
#' @param msigdb_cat Category for msibdg, see msigdbr_collections()
#' @param msigdb_subcat msigdb subcategory
#' @param ASC Allele Specific categorization in "biallelic/monoallelic/other" column identifier in df.
#' @param gene Gene-id column identifier in df.
#'
#' @import msigdbr
#' @import tidyverse
#' @importFrom reshape2 dcast
#'
#' @return A ggplot
#' @export
#'
plot_msigdb<- function(df,msigdb_species="human",msigdb_cat="C5",msigdb_subcat="BP",ASC="type",gene="ensembl_gene_id"){
  msigdb_gene_sets = msigdbr(species = msigdb_species, category = msigdb_cat,subcategory = msigdb_subcat)

  df_ann<-merge(df,msigdb_gene_sets[,c("human_ensembl_gene","gs_name")],by.x=c(gene),by.y=c("human_ensembl_gene"))
  df_ann<-unique(df_ann[,c(..gene,..ASC,"gs_name")])

  gs<-merge(reshape2::dcast(df_ann,gs_name~get(ASC), fun=length, value.var = 'gs_name')%>%top_n(10,Biallelic)%>%dplyr::select(gs_name),
            reshape2::dcast(df_ann,gs_name~get(ASC), fun=length, value.var = 'gs_name')%>%top_n(10,Monoallelic)%>%dplyr::select(gs_name),all=TRUE)

  df_ann_plot<-df_ann[df_ann$gs_name%in%gs$gs_name,]
  plot_annotation(df=df_ann_plot,ASC,gene,annotation="gs_name")
}
