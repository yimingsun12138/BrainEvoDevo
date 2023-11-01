#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: export DAP analysis                                             ##
## Data: 2023.08.24                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.1/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(dplyr)
library(Seurat)
library(ggplot2)
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
library(OpenAI4R)
library(paletteer)
library(ggpattern)
library(ggrepel)
library(openxlsx)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0314')

# get cell type list ------------------------------------------------------
cell_type_list <- list.files('./res/step_118_fig_230422/Wilcox')
cell_type_list <- gsub(pattern = '_df_wilcox.rds',replacement = '',x = cell_type_list,fixed = TRUE)
names(cell_type_list) <- c('End','ExN','ExM','ExUp','ExDp','InCGE','InMGE','IP','Mic','Per','RG-neu')
cell_type_list <- cell_type_list[c("RG-neu","IP","ExN","ExM","ExUp","ExDp")]

# wilcox ------------------------------------------------------------------
work_sheet <- createWorkbook()

for (i in names(cell_type_list)) {
  cell_type <- cell_type_list[i]
  temp <- paste0('./res/step_118_fig_230422/Wilcox/',cell_type,'_df_wilcox.rds')
  temp <- readRDS(temp)
  temp <- temp[,c("peak","compare","log2FC","pval","fdr","group")]
  temp <- temp[which(temp$group != 'others'),]
  
  addWorksheet(wb = work_sheet,sheetName = i)
  writeDataTable(wb = work_sheet,sheet = i,x = temp)
}

saveWorkbook(wb = work_sheet,file = './res/step_130_fig_230824/DAP_analysis_by_Wilcox.xlsx',overwrite = TRUE)

# DESeq2 ------------------------------------------------------------------
work_sheet <- createWorkbook()

for (i in names(cell_type_list)) {
  cell_type <- cell_type_list[i]
  temp <- paste0('./res/step_118_fig_230422/DESeq2_random/',cell_type,'_df_DESeq2.rds')
  temp <- readRDS(temp)
  temp <- temp[,c("peak","compare","log2FC","pvalue","fdr","group")]
  colnames(temp) <- c("peak","compare","log2FC","pval","fdr","group")
  temp <- temp[which(temp$group != 'others'),]
  
  addWorksheet(wb = work_sheet,sheetName = i)
  writeDataTable(wb = work_sheet,sheet = i,x = temp)
}

saveWorkbook(wb = work_sheet,file = './res/step_130_fig_230824/DAP_analysis_by_DESeq2.xlsx',overwrite = TRUE)
