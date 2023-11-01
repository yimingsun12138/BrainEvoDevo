#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: modern human fixed SNC genome distribution                      ##
## Data: 2022.10.22                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#sleep
ii <- 1
while(1){
  i <- matrix(data = rnorm(n = 2000*2000),nrow = 2000,ncol = 2000)
  i <- i %*% i
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#general setting
setwd('/content/data/sunym/project/Brain/')
.libPaths('/content/data/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(Seurat)
library(ggplot2)
library(dplyr)
library(scibet)
library(Matrix)
library(tidyverse)
library(cowplot)
library(viridis)
library(ComplexHeatmap)
library(parallel)
library(ggsignif)
library(RColorBrewer)
library(ggsci)
library(scales)
library(patchwork)
library(ggpointdensity)
library(latex2exp)
library(ArchR)
library(scales)
library(circlize)
library(ggpubr)
library(ggtext)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(DESeq2)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# global distribution -----------------------------------------------------
#load data
SNC_GRanges <- readRDS(file = './res/step_67_fig_221021/SNC_GRanges.rds')

#make Txdb object
gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
gene_anno <- rtracklayer::as.data.frame(gene_anno)
gene_anno$seqnames <- as.character(gene_anno$seqnames)
gene_anno$seqnames <- unlist(base::lapply(X = gene_anno$seqnames,FUN = function(x){
  if(x %in% c(as.character(1:22),'X','Y','MT')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
}))
gene_anno <- makeGRangesFromDataFrame(df = gene_anno,keep.extra.columns = TRUE)
gene_anno <- GenomicFeatures::makeTxDbFromGRanges(gene_anno)

#plot
peakAnno <- annotatePeak(SNC_GRanges,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")

pdf(file = './res/step_68_fig_221022/human_fixed_SNC_genome_distribution.pdf',width = 6,height = 4)
plotAnnoPie(peakAnno)
dev.off()

# non-extent SNC genome distribution --------------------------------------
#load data
SNC_GRanges <- readRDS(file = './res/step_67_fig_221021/SNC_GRanges.rds')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#make Txdb object
gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
gene_anno <- rtracklayer::as.data.frame(gene_anno)
gene_anno$seqnames <- as.character(gene_anno$seqnames)
gene_anno$seqnames <- unlist(base::lapply(X = gene_anno$seqnames,FUN = function(x){
  if(x %in% c(as.character(1:22),'X','Y','MT')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
}))
gene_anno <- makeGRangesFromDataFrame(df = gene_anno,keep.extra.columns = TRUE)
gene_anno <- GenomicFeatures::makeTxDbFromGRanges(gene_anno)

## active SNC in each cell type --------------------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

SNC_list <- base::lapply(X = cell_type_list,FUN = function(x){
  return(SNC_GRanges[SNC_GRanges@elementMetadata[,x] == 'YES'])
})
cell_type_list <- gsub(pattern = '.',replacement = '-',x = cell_type_list,fixed = TRUE)
names(SNC_list) <- cell_type_list

#plot
peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_68_fig_221022/non_extent_cell_type_active_SNC_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

## SNC with different active intersect size --------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- as.data.frame(SNC_GRanges@elementMetadata[,cell_type_list])
temp <- rowSums(temp == 'YES')
SNC_GRanges$intersect <- temp

SNC_list <- base::lapply(X = 0:13,FUN = function(x){
  return(SNC_GRanges[SNC_GRanges$intersect == x])
})
names(SNC_list) <- 0:13

#plot
peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_68_fig_221022/non_extent_SNC_with_intersect_size_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

# extent SNC genome distribution --------------------------------------
#load data
SNC_GRanges <- readRDS(file = './res/step_67_fig_221021/extent_plot/SNC_GRanges.rds')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#make Txdb object
gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
gene_anno <- rtracklayer::as.data.frame(gene_anno)
gene_anno$seqnames <- as.character(gene_anno$seqnames)
gene_anno$seqnames <- unlist(base::lapply(X = gene_anno$seqnames,FUN = function(x){
  if(x %in% c(as.character(1:22),'X','Y','MT')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
}))
gene_anno <- makeGRangesFromDataFrame(df = gene_anno,keep.extra.columns = TRUE)
gene_anno <- GenomicFeatures::makeTxDbFromGRanges(gene_anno)

## active SNC in each cell type --------------------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

SNC_list <- base::lapply(X = cell_type_list,FUN = function(x){
  return(SNC_GRanges[SNC_GRanges@elementMetadata[,x] == 'YES'])
})
cell_type_list <- gsub(pattern = '.',replacement = '-',x = cell_type_list,fixed = TRUE)
names(SNC_list) <- cell_type_list

#plot
peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_68_fig_221022/extent_cell_type_active_SNC_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

## SNC with different active intersect size --------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- as.data.frame(SNC_GRanges@elementMetadata[,cell_type_list])
temp <- rowSums(temp == 'YES')
SNC_GRanges$intersect <- temp

SNC_list <- base::lapply(X = 0:13,FUN = function(x){
  return(SNC_GRanges[SNC_GRanges$intersect == x])
})
names(SNC_list) <- 0:13

#plot
peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_68_fig_221022/extent_SNC_with_intersect_size_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()