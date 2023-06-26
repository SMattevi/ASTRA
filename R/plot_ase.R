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
      ggplot2::geom_point(aes(x=1.2,size=.data$totalCount),color="gray")+
      ggplot2::ggtitle(title)+
      ggplot2::scale_fill_discrete(name = "", labels = c("A1_counts", "A2_counts"))+
      ggplot2::geom_vline(xintercept = .85)+ggplot2::geom_vline(xintercept = .15)
    #  color = "black") #+#+
    #ggplot2::facet_grid(~cluster)
    return(plot_res)
  }
}
