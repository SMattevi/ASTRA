#' Read ASE info and merge with vcf info. Check for mismatching alleles and prepare dataframe for manual phasing.
#'
#' @param df_var A df obtained from snps_read_and_merge()/conshap()
#' @param ASE_path A path of the ASE obtained from gex bulk analysis
#' @param BM A df for SNPs annotation
#' @param reads Number of minimum reads overlapping a SNPs
#' @param refAllele refAllele ASE
#' @param altAllele altAllele ASE
#' @param position position ASE
#' @param refCount refCount ASE
#' @param altCount altCount ASE
#' @param contig contig ASE
#' @param totalCount totalCount ASE
#' @param expand_region expand gene region by indicated bp
#' @param od overall depth, default ≥ 10
#'
#'
#' @importFrom data.table fread
#' @importFrom IRanges IRanges
#' @importFrom IRanges from
#' @importFrom IRanges to
#' @importFrom GenomicRanges findOverlaps
#' @importFrom stats reorder
#' @import ggplot2
#' @return A dataframe.
#'
#' @export
#'
#'
#'
addASE<- function(df_var,ASE_path,BM,reads=2,od=10,expand_region=0,refAllele="refAllele",altAllele="altAllele",position="position",refCount="refCount",altCount="altCount",contig="contig", totalCount="totalCount"){

  ASE_in<- data.table::fread(ASE_path)
  ASE_in[[contig]]<-as.factor(ASE_in[[contig]])
  snpsr<-base::merge(ASE_in,df_var,by.x=c(contig,position),by.y=c("CHROM","POS"))

  #check for eventually mismatching alleles
  snpsturned<-base::merge(ASE_in,df_var,by.x=c(contig,position,refAllele,altAllele),by.y=c("CHROM","POS","ALT","REF"))
  if(dim(snpsturned)[1]>0){
    print("mismatching phasing and ASE information!")
    stop()
  }
  #2. flip according to haplotype REF and ALT allele
  for(i in rownames(snpsr)){
    i<-as.integer(i)
    if(!is.na(snpsr$Ref[i])&snpsr$Ref[i]=="1"){
      snpsr<-flip_alleles(snpsr,refAllele,altAllele,i)
      snpsr<-flip_alleles(snpsr,refCount,altCount,i)
    }
  }

  #variantIDs creation as chr:pos:refAllele_altAllele
  snpsr$variantID<-paste0(snpsr[[contig]],":",snpsr[[position]],":",snpsr[[refAllele]],"_",snpsr[[altAllele]])

  #Gene annotation and extraction of SNPs overlapping a gene
  isnps <- with(snpsr, IRanges::IRanges(position, width=1, names=variantID))
  igenes <- with(BM, IRanges::IRanges(start_position-expand_region, end_position+expand_region, names=ensembl_gene_id))
  olaps <- GenomicRanges::findOverlaps(isnps, igenes)

  overlapping<-cbind(snpsr[IRanges::from(olaps),], BM[IRanges::to(olaps),])

  #3. keep only SNPs overlapping a gene with at least 2 reads supporting ref or alt count
  minmax<- dplyr::filter(overlapping, refCount>=reads & altCount>=reads, )
  odpass<- dplyr::filter(minmax, totalCount>=od, )

  #prepare df to be merged with other clusters
  df<-odpass
  df$refFrac<-df[[refCount]]/df[[totalCount]]
  df<-df[,c(..position,..contig,"variantID","Ref","hgnc_symbol","ensembl_gene_id","refFrac","totalCount",..refCount,..altCount,..refAllele,..altAllele)]

  snp_df<-data.frame(SNPS=as.factor(c("Processed (has ASEReadCount)",paste0("With min ref/alt counts ",reads),paste0("With overall depth >= ",od),"Overlapping at least 1 genes")),NUM=c(length(unique(snpsr$variantID)),length(unique(minmax$variantID)),length(unique(odpass$variantID)),length(unique(overlapping$variantID))))

  plot_snp<-ggplot2::ggplot(snp_df, ggplot2::aes(x = reorder(.data$SNPS,-.data$NUM), y = .data$NUM ,fill = .data$SNPS )) +
    ggplot2::geom_bar(width = 0.85, stat="identity")+ggplot2::coord_polar(theta = "y") +
    ggplot2::geom_text(aes(label=.data$NUM),position = position_stack(vjust = 0.5),color="white") +
    ggplot2::theme_void()
  print(plot_snp)

  return(df)
}
