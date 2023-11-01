#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: create mouse multiome human symbol Seurat object                ##
## Data: 2022.09.23                                                                ##
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

# convert to human symbol -------------------------------------------------
mouse_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human_anno <- mouse_to_human_anno[,c(2,1,4,3)]
express_matrix <- mouse_multiome_Seurat@assays$RNA@counts
express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = express_matrix,anno = mouse_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 200*(1024^3),workers = 6)

# create seurat object ----------------------------------------------------
meta_data <- mouse_multiome_Seurat@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_RNA','nFeature_RNA','seurat_clusters'))]
temp <- colnames(meta_data)[!(grepl(pattern = 'RNA_snn_res',x = colnames(meta_data),fixed = TRUE))]
meta_data <- meta_data[colnames(express_matrix),temp]
mouse_multiome_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'mouse',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

# process -----------------------------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 18,resolution = 1,group.by = 'cell_type',label = TRUE)

#save data
saveRDS(object = mouse_multiome_Seurat,file = './processed_data/220802_summary/mouse_multiome_Seurat_human_symbol_220923.rds')
