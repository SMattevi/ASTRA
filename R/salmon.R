#' Read and merge vcf containing all and only (Shapeit4) phased variants
#'
#' @param salmon_path A character vector with quant.sf path.
#' @param txdb txdb for salmon import
#' @param counts Filtering value on salmon counts
#'
#' @import tximport
#'
#' @return A dataframe.
#' @export
#'
read_salmon<- function(salmon_path,txdb,counts){
  k <- keys(txdb, keytype = "TXID")
  tx2gene <- select(txdb, k, "GENEID", "TXID")

  txi.salmon <- tximport::tximport(salmon_path, type = "salmon", tx2gene = tx2gene, ignoreAfterBar=T)
  salmon_counts <- as.data.frame(txi.salmon$counts)
  colnames(salmon_counts) <- "counts"

  return(salmon_counts)
}
#' SEESAW analysis plotting
#'
#' @param names A character vector with the sample id
#' @param edb_id ensembl database id from AnnotationHub i.e. "AH113665"
#' @param tss_gap gap size for tss grouping
#' @param gene_id ensebl gene id to plot
#' @param files path to quant.sf file from the seesaw pipeline
#' @param cellline_id cellular line description i.e. "naive"
#'
#' @import AnnotationHub
#' @import fishpond
#' @import SummarizedExperiment
#' @import plyranges
#' @import Gviz
#'
#' @return A plot.
#' @export
#'
seesaw<-function(names, edb_id="AH113665",tss_gap=50,gene_id,files, cellline_id){
  cellline<-c(cellline_id)

  coldata <- data.frame(names,files,cellline)

  #query annotation for Homo Sapiens, download and save in edb
  ah <- AnnotationHub()
  edb <- ah[[edb_id]]

  #extract transcripts
  txps <- transcripts(edb)
  #group transcripts by TSS (i.e. 50 bp gap)
  txps <- makeTx2Tss(edb, maxgap=tss_gap) %>%
    select(tx_id, gene_id, group_id, tss)

  #Effectively import allelic counts obtained after g2gtools and salmon. Change a1 and a2!
  gse <- importAllelicCounts(
    coldata, a1="L", a2="R",
    format="wide", tx2gene=txps, importer=read.delim
  )

  # filtering out lowly expressed features (assay(gse) = expression for each allele, rowsums = number of allele in which it has to be verified):
  keep <- rowSums(assay(gse) >= 10) >= 2
  gse <- gse[keep,]

  #The counts assay provides the read counts for features x samples,
  #as Salmon is used for quantification these are estimated counts that are non-negative but fractional,
  #not integer. The abundance assay gives estimates in Transcripts Per Million (TPM).
  #The length assay gives the average effective transcript length for a given transcript or a given gene.
  #Finally the infRep’s are the bootstrap samples for the counts.
  assayNames(gse)

  #plot isoform specific expression
  plotAllelicGene(gse, gene=gene_id,
                  db=edb,cov="cellline",qvalue = FALSE,log2FC = FALSE)
}
