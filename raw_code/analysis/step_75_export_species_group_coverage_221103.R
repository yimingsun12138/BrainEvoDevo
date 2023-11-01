#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: export species group coverage                                   ##
## Data: 2022.11.03                                                                ##
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
library(org.Mmu.eg.db)
library(DESeq2)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# using reads in peaks ----------------------------------------------------

## human -------------------------------------------------------------------
#load data
ArchR_project <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

#create temp
temp <- getArrowFiles(ArchRProj = ArchR_project)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/Greenleaf_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = ArchR_project),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = ArchR_project)
)
temp <- temp[rownames(ArchR_project@cellColData)]
temp$cell_type <- ArchR_project$cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = getPeakSet(ArchRProj = ArchR_project))
temp <- addPeakMatrix(ArchRProj = temp)

#export coverage
temp$ReadsInPromoter <- temp$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = temp,groupBy = 'cell_type',normMethod = 'ReadsInPromoter',tileSize = 25,maxCells = NULL,verbose = TRUE)
ii <- gsub(pattern = 'ReadsInPromoter',replacement = 'ReadsInPeaks',fixed = TRUE,x = coverage_file)
for (i in 1:length(coverage_file)) {
  char <- paste('mv',coverage_file[i],ii[i],sep = ' ')
  system(char)
}


## macaque -----------------------------------------------------------------
#load data
ArchR_project <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

#create temp
temp <- getArrowFiles(ArchRProj = ArchR_project)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = ArchR_project),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = ArchR_project)
)
temp <- temp[rownames(ArchR_project@cellColData)]
temp$cell_type <- ArchR_project$cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = getPeakSet(ArchRProj = ArchR_project))
temp <- addPeakMatrix(ArchRProj = temp)

#export coverage
temp$ReadsInPromoter <- temp$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = temp,groupBy = 'cell_type',normMethod = 'ReadsInPromoter',tileSize = 25,maxCells = NULL,verbose = TRUE)
ii <- gsub(pattern = 'ReadsInPromoter',replacement = 'ReadsInPeaks',fixed = TRUE,x = coverage_file)
for (i in 1:length(coverage_file)) {
  char <- paste('mv',coverage_file[i],ii[i],sep = ' ')
  system(char)
}

## mouse -------------------------------------------------------------------
#load data
ArchR_project <- loadArchRProject(path = './processed_data/221008_summary/mouse_multiome_ArchR_221009/')

#create temp
temp <- getArrowFiles(ArchRProj = ArchR_project)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/mouse_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = ArchR_project),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = ArchR_project)
)
temp <- temp[rownames(ArchR_project@cellColData)]
temp$cell_type <- ArchR_project$Gex_macaque_cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = getPeakSet(ArchRProj = ArchR_project))
temp <- addPeakMatrix(ArchRProj = temp)

#export coverage
temp$ReadsInPromoter <- temp$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = temp,groupBy = 'cell_type',normMethod = 'ReadsInPromoter',tileSize = 25,maxCells = NULL,verbose = TRUE)
ii <- gsub(pattern = 'ReadsInPromoter',replacement = 'ReadsInPeaks',fixed = TRUE,x = coverage_file)
for (i in 1:length(coverage_file)) {
  char <- paste('mv',coverage_file[i],ii[i],sep = ' ')
  system(char)
}

# use reads in TSS --------------------------------------------------------

## human -------------------------------------------------------------------
setwd('/home/sunym')
#load data
ArchR_project <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

#create temp
temp <- getArrowFiles(ArchRProj = ArchR_project)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/Greenleaf_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = ArchR_project),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = ArchR_project)
)
temp <- temp[rownames(ArchR_project@cellColData)]
temp$cell_type <- ArchR_project$cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = getPeakSet(ArchRProj = ArchR_project))
temp <- addPeakMatrix(ArchRProj = temp)

#export coverage
coverage_file <- getGroupBW(ArchRProj = temp,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)

## macaque -----------------------------------------------------------------
setwd('/home/sunym')
#load data
ArchR_project <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')

#create temp
temp <- getArrowFiles(ArchRProj = ArchR_project)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = ArchR_project),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = ArchR_project)
)
temp <- temp[rownames(ArchR_project@cellColData)]
temp$cell_type <- ArchR_project$cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = getPeakSet(ArchRProj = ArchR_project))
temp <- addPeakMatrix(ArchRProj = temp)

#export coverage
coverage_file <- getGroupBW(ArchRProj = temp,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)

## mouse -------------------------------------------------------------------
setwd('/home/sunym')
#load data
ArchR_project <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/mouse_multiome_ArchR_221009/')

#create temp
temp <- getArrowFiles(ArchRProj = ArchR_project)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/mouse_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = ArchR_project),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = ArchR_project)
)
temp <- temp[rownames(ArchR_project@cellColData)]
temp$cell_type <- ArchR_project$Gex_macaque_cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = getPeakSet(ArchRProj = ArchR_project))
temp <- addPeakMatrix(ArchRProj = temp)

#export coverage
coverage_file <- getGroupBW(ArchRProj = temp,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 25,maxCells = NULL,verbose = TRUE)
