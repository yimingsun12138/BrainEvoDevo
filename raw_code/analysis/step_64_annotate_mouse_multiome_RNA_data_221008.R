#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: annotate mouse multiome RNA data                                ##
## Data: 2022.10.08                                                                ##
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
library(Rmisc)
library(Seurat)
library(SingleCellExperiment)
library(sctransform)
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
library(patchwork)
library(networkD3)
library(htmlwidgets)
library(circlize)
library(harmony)
library(ArchR)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# load data ---------------------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/mouse_multiome_Seurat_220922.rds')
mouse_multiome_Seurat_human_symbol <- readRDS(file = './processed_data/220802_summary/mouse_multiome_Seurat_human_symbol_220923.rds')
cell_list <- readRDS(file = './ArchR/res/step_15_fig_221003/cell_list.rds')

# filter and save data ----------------------------------------------------
mouse_multiome_Seurat <- mouse_multiome_Seurat[,cell_list]
mouse_multiome_Seurat_human_symbol <- mouse_multiome_Seurat_human_symbol[,cell_list]

#save data
saveRDS(object = mouse_multiome_Seurat,file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
saveRDS(object = mouse_multiome_Seurat_human_symbol,file = './processed_data/221008_summary/mouse_multiome_Seurat_human_symbol_221009.rds')

# display data ------------------------------------------------------------
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
temp_col <- col_param$celltype.mouse

p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cols = temp_col) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'cell_type',label = TRUE,repel = TRUE,cols = temp_col) + theme(aspect.ratio = 1)
p4 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+p4+plot_layout(ncol = 2)
