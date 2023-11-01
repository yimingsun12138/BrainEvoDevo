#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: zooHAR scATAC profile                                           ##
## Data: 2023.05.29                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.0/')
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

# process zooHAR ----------------------------------------------------------
zooHAR_file <- read.csv(file = './res/step_120_fig_230610/zooHARs.csv')

temp <- zooHAR_file[,c("chrom","start","end")]
colnames(temp) <- c('chrom','start','end')
zooHAR <- as(temp,'GRanges')

mcols(zooHAR) <- zooHAR_file[,c("simple_name","name","translated_category","Genoa_score","Genoa_score_1400bp")]

#save zooHAR
saveRDS(object = zooHAR,file = './res/step_120_fig_230610/zooHARs.rds')

# overlap with human cell type peaks --------------------------------------------

#load ArchR object
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#cell type list
cell_type_lict <- unique(x = Greenleaf_ATAC_ArchR$cell_type)
cell_type_dot <- sub(pattern = '-',replacement = '.',x = cell_type_lict,fixed = TRUE)

#load peak_list
for (i in cell_type_dot) {
  cell_type <- i
  peak_file <- list.files(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls')
  peak_file <- peak_file[grep(pattern = cell_type,x = peak_file,fixed = TRUE)]
  peak_file <- paste('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls',peak_file,sep = '/')
  peak_file <- readRDS(file = peak_file)
  
  temp_char <- paste(cell_type,'peak',sep = '_')
  assign(x = temp_char,value = peak_file)
}

#overlap peak numbers
peak_number <- base::unlist(base::lapply(X = cell_type_dot,FUN = function(x){
  #get peak file
  cell_type <- x
  temp_char <- paste(cell_type,'peak',sep = '_')
  peak_file <- get(x = temp_char)
  
  #count overlaps
  temp <- sum(countOverlaps(query = peak_file,subject = zooHAR) > 0)
  return(temp)
}))

#ggplot
plot_table <- data.frame(cell_type = cell_type_lict,peak_number = peak_number)
plot_table[which(plot_table$cell_type == 'RG-1'),"cell_type"] <- 'RG-neu'
plot_table[which(plot_table$cell_type == 'RG-2'),"cell_type"] <- 'RG-glia'
plot_table[which(plot_table$cell_type == 'Ex-1'),"cell_type"] <- 'ExN'
plot_table[which(plot_table$cell_type == 'Ex-2'),"cell_type"] <- 'ExM'
plot_table[which(plot_table$cell_type == 'Ex-3'),"cell_type"] <- 'ExUp'
plot_table[which(plot_table$cell_type == 'Ex-4'),"cell_type"] <- 'ExDp'

table(plot_table$cell_type %in% names(color_param$unify_cell_type))
plot_table$cell_type <- factor(plot_table$cell_type,levels = names(color_param$unify_cell_type)[names(color_param$unify_cell_type) %in% plot_table$cell_type])

