#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: motif enrichment analysis on DAP and HAR                        ##
## Data: 2023.02.14                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.2.2/')
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
library(BSgenome.Mmusculus.UCSC.mm10)
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(universalmotif)
library(topGO)
library(future.apply)
library(transPlotR)
library(aplot)
library(Biostrings)
library(riverplot)
library(CellChat)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

# function to export bed file ---------------------------------------------
temp_fun <- function(DAP_list,group,method = 'wilcox'){
  #convert DAP list to bed file
  
  cell_type <- DAP_list
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #file file
  DAP_list <- list.files(path = './res/step_74_fig_221104')
  DAP_list <- DAP_list[grep(pattern = method,x = DAP_list,fixed = TRUE)]
  DAP_list <- DAP_list[grep(pattern = cell_type_dot,x = DAP_list,fixed = TRUE)]
  DAP_list <- paste('./res/step_74_fig_221104',DAP_list,sep = '/')
  DAP_list <- readRDS(file = DAP_list)
  
  #filter DAP list
  DAP_list <- names(DAP_list)[DAP_list == group]
  DAP_list <- base::do.call(what = rbind,args = base::lapply(X = DAP_list,FUN = function(x){
    x <- strsplit(x = x,split = '-')
    return(x[[1]])
  }))
  colnames(DAP_list) <- c('chrom','start','end')
  DAP_list <- as.data.frame(DAP_list)
  DAP_list <- as(DAP_list,'GRanges')
  
  #export
  char <- paste0('./res/step_103_fig_230214/',cell_type_dot,'_',group,'_DAP_by_',method,'.bed')
  export.bed(object = DAP_list,con = char)
  
  #return
  print(paste(cell_type,'done!',sep = ' '))
  return(NULL)
}

# export wilcox based DAP -------------------------------------------------

#cell type list
cell_type_list <- list.files(path = './res/step_74_fig_221104')
cell_type_list <- gsub(pattern = '_DAP_list.*.rds$',replacement = '',x = cell_type_list,fixed = FALSE)
cell_type_list <- unique(gsub(pattern = '.',replacement = '-',x = cell_type_list,fixed = TRUE))

group_list <- c('human_specific','human_macaque_conserved','species_conserved')

for (i in cell_type_list) {
  for (j in group_list) {
    temp_fun(DAP_list = i,group = j,method = 'wilcox')
  }
}
