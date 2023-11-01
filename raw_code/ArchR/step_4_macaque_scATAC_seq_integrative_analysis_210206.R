#####################################################################################
## Project: macaque fetal Brain scATAC_seq                                         ##
## Script Purpose: macaque scATAC_seq integrative analysis                         ##
## Data: 2022.02.05                                                                ##
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
setwd('/data/User/sunym/project/Brain/ArchR/')
.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.3/',
            '/data/User/sunym/software/R_lib/R_4.1.3/'))
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
library(harmony,lib.loc = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.3/')
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(rtracklayer)
library(org.Mmu.eg.db)
library(clusterProfiler)
library(circlize)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

#load data
macaque_ATAC_ArchR <- readRDS(file = './processed_data/macaque_ATAC_ArchR_220109/Save-ArchR-Project.rds')


# peak co-accessibility ----------------------------------------------------
macaque_ATAC_ArchR <- addCoAccessibility(
  ArchRProj = macaque_ATAC_ArchR,
  reducedDims = "IterativeLSI"
)

my_send_sms('co-accessibility calculate done!')

cA <- getCoAccessibility(
  ArchRProj = macaque_ATAC_ArchR,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA

markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(macaque_ATAC_ArchR,resolution = 2000),
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG')
)

#RG
pdf(file = './res/step_4_fig_220206/macaque_ATAC_RG_VIM_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$VIM)
dev.off()

#RG-3
pdf(file = './res/step_4_fig_220206/macaque_ATAC_RG_3_EGFR_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$EGFR)
dev.off()

#OPC
pdf(file = './res/step_4_fig_220206/macaque_ATAC_OPC_OLIG2_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$OLIG2)
dev.off()

#IP
pdf(file = './res/step_4_fig_220206/macaque_ATAC_IP_EOMES_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$EOMES)
dev.off()

#Ex
pdf(file = './res/step_4_fig_220206/macaque_ATAC_Ex_NEUROD6_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$NEUROD6)
dev.off()

pdf(file = './res/step_4_fig_220206/macaque_ATAC_Ex_FEZF2_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$FEZF2)
dev.off()

#In
pdf(file = './res/step_4_fig_220206/macaque_ATAC_In_DLX5_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$DLX5)
dev.off()

pdf(file = './res/step_4_fig_220206/macaque_ATAC_In_LHX6_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$LHX6)
dev.off()

pdf(file = './res/step_4_fig_220206/macaque_ATAC_In_SP8_trackplot_coaccessibility.pdf',width = 12,height = 6)
grid::grid.newpage()
grid::grid.draw(p$SP8)
dev.off()

# peak and gene express correlation ---------------------------------------
macaque_ATAC_ArchR <- addPeak2GeneLinks(
  ArchRProj = macaque_ATAC_ArchR,
  reducedDims = "IterativeLSI",
  useMatrix = 'GeneIntegrationMatrix'
)

my_send_sms('peak and gene correlation done!')

p2g <- getPeak2GeneLinks(
  ArchRProj = macaque_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g

markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_ATAC_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_ATAC_ArchR),
  useGroups = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG')
)

grid::grid.newpage()
grid::grid.draw(p$EOMES)

grid::grid.newpage()
grid::grid.draw(p$SOX9)

p <- plotPeak2GeneHeatmap(ArchRProj = macaque_ATAC_ArchR, groupBy = "cell_type", k = 10)

pdf(file = './res/step_4_fig_220206/macaque_ATAC_RNA_EP_link_heatmap.pdf',width = 8,height = 10)
ComplexHeatmap::draw(p, heatmap_legend_side = "bot", annotation_legend_side = "bot", width = unit(2,'inches'), heigh = unit(8,'inches'))
dev.off()

# save data ---------------------------------------------------------------
saveArchRProject(ArchRProj = macaque_ATAC_ArchR, outputDirectory = "./processed_data/macaque_ATAC_ArchR_220109/", load = FALSE, overwrite = TRUE)
