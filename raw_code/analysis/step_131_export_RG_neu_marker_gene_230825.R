#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: export RG-neu marker gene                                       ##
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

# load marker -------------------------------------------------------------
Human_specific_signature <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
Primate_conserved_signature <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
Species_conserved_signature <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')

gene_list <- c(rep('Human specific',times = length(Human_specific_signature)),
               rep('Primate specific',times = length(Primate_conserved_signature)),
               rep('Species conserved',times = length(Species_conserved_signature)))
names(gene_list) <- c(Human_specific_signature,Primate_conserved_signature,Species_conserved_signature)

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

meta_data <- readRDS(file = './res/step_96_fig_230104/mouse_multiome_annotation.rds')
mouse_multiome_Seurat$ReAnno_celltype <- meta_data[colnames(mouse_multiome_Seurat),"ReAnno_celltype"]

# human -------------------------------------------------------------------
species_label <- 'human'
temp_Seurat <- Greenleaf_RNA_Seurat
temp_Seurat$cell_type <- Greenleaf_RNA_Seurat$ReAnno_celltype

cell_type_list <- unique(Greenleaf_RNA_Seurat$ReAnno_celltype)
cell_type_list <- c('RG-1','IP','Ex-1','Ex-2','Ex-3','Ex-4')
names(cell_type_list) <- c('RG-neu','IP','ExN','ExM','ExUp','ExDp')

stat_table <- base::do.call(what = rbind,args = base::lapply(X = names(gene_list),FUN = function(x){
  exp_pct <- base::lapply(X = cell_type_list,FUN = function(y){
    all_cell <- sum(temp_Seurat$cell_type == y)
    expressed_cell <- sum(temp_Seurat@assays$RNA@counts[x,temp_Seurat$cell_type == y] > 0)
    return(expressed_cell/all_cell)
  })
  exp_pct <- unlist(exp_pct) * 100
  names(exp_pct) <- names(cell_type_list)
  
  expression_level <- base::lapply(X = cell_type_list,FUN = function(y){
    return(mean(exp(temp_Seurat@assays$RNA@data[x,temp_Seurat$cell_type == y])))
  })
  expression_level <- unlist(expression_level)
  expression_level <- round(expression_level,digits = 2)
  names(expression_level) <- names(cell_type_list)
  
  temp <- c(x,gene_list[x],
            exp_pct['RG-neu'],max(exp_pct[names(exp_pct) != 'RG-neu']),
            expression_level['RG-neu'],max(expression_level[names(expression_level) != 'RG-neu']))
  names(temp) <- c('gene','group','target Pct %','max non-target Pct %','target Exp','max non-target Exp')
  return(temp)
}))
