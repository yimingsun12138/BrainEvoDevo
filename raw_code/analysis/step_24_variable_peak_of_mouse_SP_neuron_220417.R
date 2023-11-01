#####################################################################################
## Project: mouse fetal Brain scATAC_seq                                           ##
## Script Purpose: variable peak of mouse SP neuron                                ##
## Data: 2022.04.17                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#general setting
setwd('/data/User/sunym/project/Brain/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(parallel)
library(ArchR)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(harmony)
library(BSgenome.Mmusculus.UCSC.mm10)
library(rtracklayer)
library(org.Mm.eg.db)
library(clusterProfiler)
library(scibet)
library(viridis)
library(networkD3)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')
p1 <- plotEmbedding(ArchRProj = mouse_ATAC_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = mouse_ATAC_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+plot_layout(ncol = 2)

# calculate variable peak -------------------------------------------------
getAvailableMatrices(mouse_ATAC_ArchR)
markersPeaks <- getMarkerFeatures(
  ArchRProj = mouse_ATAC_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment","log10(nFrags)"),
  testMethod = "wilcoxon"
)
markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE)
markerList

Ex_3_marker <- markerList$`Ex-3`

# public peakset ----------------------------------------------------------

## liuyt_H3K27ac_peakset ---------------------------------------------------
# load public peakset
liuyt_H3K27ac_peakset <- read.table(file = './data/public/Multiscale_3D_Genome_Rewiring_during_Mouse_Neural_Development/liuyt/mouse.Enhancer.GSE96107_CN_H3K27ac.IDR0.05.filt.narrowPeak.subtract.TSSup2kb.bed',header = FALSE)
liuyt_H3K27ac_peakset <- liuyt_H3K27ac_peakset[,1:3]
colnames(liuyt_H3K27ac_peakset) <- c('chrom','start','end')
temp <- paste(liuyt_H3K27ac_peakset$chrom,liuyt_H3K27ac_peakset$start,liuyt_H3K27ac_peakset$end,sep = '-')
temp <- c(!duplicated(temp))
liuyt_H3K27ac_peakset <- liuyt_H3K27ac_peakset[temp,]
temp <- paste(liuyt_H3K27ac_peakset$chrom,liuyt_H3K27ac_peakset$start,liuyt_H3K27ac_peakset$end,sep = '-')
rownames(liuyt_H3K27ac_peakset) <- temp
liuyt_H3K27ac_peakset <- as(liuyt_H3K27ac_peakset,'GRanges')

#overlap with public peakset
temp <- GenomicRanges::countOverlaps(query = Ex_3_marker,subject = liuyt_H3K27ac_peakset)
Ex_3_marker <- Ex_3_marker[temp > 0]
temp <- paste(Ex_3_marker@seqnames,as.character(Ex_3_marker@ranges),sep = '-')
names(Ex_3_marker) <- temp

mouse_peakset <- getPeakSet(mouse_ATAC_ArchR)
temp <- paste(mouse_peakset@seqnames,as.character(mouse_peakset@ranges),sep = '-')
names(mouse_peakset) <- temp
Ex_3_marker$nearestGene <- mouse_peakset[names(Ex_3_marker)]$nearestGene

#track plot
my_trackplot <- function(Granges_obj,ext = 250000){
  temp <- liuyt_H3K27ac_peakset[countOverlaps(query = liuyt_H3K27ac_peakset,subject = Granges_obj) >= 1]
  Granges_obj <- rtracklayer::as.data.frame(Granges_obj)
  Granges_obj$start <- Granges_obj$start-ext
  Granges_obj$end <- Granges_obj$end+ext
  if(Granges_obj$start <= 0){
    Granges_obj$start <- 0
  }
  Granges_obj <- Granges_obj[,1:3]
  colnames(Granges_obj) <- c('chrom','start','end')
  Granges_obj <- as(Granges_obj,'GRanges')
  p <- plotBrowserTrack(
    ArchRProj = mouse_ATAC_ArchR, 
    groupBy = "cell_type", 
    region = Granges_obj, 
    loops = getPeak2GeneLinks(mouse_ATAC_ArchR,resolution = 1000),
    useGroups = c('Ex-3','Ex-1','IP','RG'),
    features = temp
  )
  return(p)
}

my_save_trackplot <- function(i,Granges_obj,ext = 250000){
  p <- my_trackplot(Granges_obj = Granges_obj,ext = ext)
  char <- paste(Granges_obj@seqnames,as.character(Granges_obj@ranges),sep = '-')
  path <- paste('./res/step_24_fig_220417',paste0(i,'_',char,'_trackplot.pdf'),sep = '/')
  pdf(file = path,width = 14,height = 7)
  grid::grid.newpage()
  grid::grid.draw(p)
  dev.off()
  Granges_obj <- append(Granges_obj,liuyt_H3K27ac_peakset[countOverlaps(query = liuyt_H3K27ac_peakset,subject = Granges_obj) >= 1])
  Granges_obj <- rtracklayer::as.data.frame(Granges_obj)
  path <- paste('./res/step_24_fig_220417',paste0(i,'_',char,'_peak.csv'),sep = '/')
  write.csv(Granges_obj,file = path)
}

for(i in 1:100){
  my_save_trackplot(i = i,Ex_3_marker[i])
}

## liuyt_H3K4me1_peakset ---------------------------------------------------
Ex_3_marker <- markerList$`Ex-3`
# load public peakset
liuyt_H3K4me1_peakset <- read.table(file = './data/public/Multiscale_3D_Genome_Rewiring_during_Mouse_Neural_Development/liuyt/mouse.Enhancer.GSE96107_CN_H3K4me1.peaks.braodPeak.subtract.TSSup2kb.bed',header = FALSE)
liuyt_H3K4me1_peakset <- liuyt_H3K4me1_peakset[,1:3]
colnames(liuyt_H3K4me1_peakset) <- c('chrom','start','end')
temp <- paste(liuyt_H3K4me1_peakset$chrom,liuyt_H3K4me1_peakset$start,liuyt_H3K4me1_peakset$end,sep = '-')
temp <- c(!duplicated(temp))
liuyt_H3K4me1_peakset <- liuyt_H3K4me1_peakset[temp,]
temp <- paste(liuyt_H3K4me1_peakset$chrom,liuyt_H3K4me1_peakset$start,liuyt_H3K4me1_peakset$end,sep = '-')
rownames(liuyt_H3K4me1_peakset) <- temp
liuyt_H3K4me1_peakset <- as(liuyt_H3K4me1_peakset,'GRanges')

#overlap with public peakset
temp <- GenomicRanges::countOverlaps(query = Ex_3_marker,subject = liuyt_H3K4me1_peakset)
Ex_3_marker <- Ex_3_marker[temp > 0]
temp <- paste(Ex_3_marker@seqnames,as.character(Ex_3_marker@ranges),sep = '-')
names(Ex_3_marker) <- temp

mouse_peakset <- getPeakSet(mouse_ATAC_ArchR)
temp <- paste(mouse_peakset@seqnames,as.character(mouse_peakset@ranges),sep = '-')
names(mouse_peakset) <- temp
Ex_3_marker$nearestGene <- mouse_peakset[names(Ex_3_marker)]$nearestGene

#track plot
my_trackplot <- function(Granges_obj,ext = 250000){
  temp <- liuyt_H3K4me1_peakset[countOverlaps(query = liuyt_H3K4me1_peakset,subject = Granges_obj) >= 1]
  Granges_obj <- rtracklayer::as.data.frame(Granges_obj)
  Granges_obj$start <- Granges_obj$start-ext
  Granges_obj$end <- Granges_obj$end+ext
  if(Granges_obj$start <= 0){
    Granges_obj$start <- 0
  }
  Granges_obj <- Granges_obj[,1:3]
  colnames(Granges_obj) <- c('chrom','start','end')
  Granges_obj <- as(Granges_obj,'GRanges')
  p <- plotBrowserTrack(
    ArchRProj = mouse_ATAC_ArchR, 
    groupBy = "cell_type", 
    region = Granges_obj, 
    loops = getPeak2GeneLinks(mouse_ATAC_ArchR,resolution = 1000),
    useGroups = c('Ex-3','Ex-1','IP','RG'),
    features = temp
  )
  return(p)
}

my_save_trackplot <- function(i,Granges_obj,ext = 250000){
  p <- my_trackplot(Granges_obj = Granges_obj,ext = ext)
  char <- paste(Granges_obj@seqnames,as.character(Granges_obj@ranges),sep = '-')
  path <- paste('./res/step_24_fig_220418',paste0(i,'_',char,'_trackplot.pdf'),sep = '/')
  pdf(file = path,width = 14,height = 7)
  grid::grid.newpage()
  grid::grid.draw(p)
  dev.off()
  Granges_obj <- append(Granges_obj,liuyt_H3K4me1_peakset[countOverlaps(query = liuyt_H3K4me1_peakset,subject = Granges_obj) >= 1])
  Granges_obj <- rtracklayer::as.data.frame(Granges_obj)
  path <- paste('./res/step_24_fig_220418',paste0(i,'_',char,'_peak.csv'),sep = '/')
  write.csv(Granges_obj,file = path)
}

for(i in 1:100){
  my_save_trackplot(i = i,Ex_3_marker[i])
}
