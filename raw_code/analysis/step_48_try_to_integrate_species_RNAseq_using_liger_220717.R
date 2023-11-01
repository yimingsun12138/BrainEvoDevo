#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: try to integrate species RNAseq using liger                     ##
## Data: 2022.07.17                                                                ##
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
library(clustree)
library(ggvenn)
library(harmony)
library(rliger)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# try to integrate using basic dataset ------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_multiome_Seurat_220712.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/220711_summary/Greenleaf_RNA_Seurat_220714.rds')
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_meta_data.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]
mouse_multiome_Seurat$cell_type <- meta_data$cell_type

#convert gene id
temp <- macaque_multiome_Seurat@assays$RNA@counts
macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
macaque_to_human_anno <- macaque_to_human_anno[,c(2,1,4,3)]
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)
macaque_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

temp <- mouse_multiome_Seurat@assays$RNA@counts
mouse_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human_anno <- mouse_to_human_anno[,c(2,1,4,3)]
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 10*(1024^3),workers = 6)
mouse_multiome_Seurat[['converted']] <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

#create liger object
Brain_RNA_liger <- createLiger(raw.data = list(macaque_multiome = macaque_multiome_Seurat@assays$converted@counts,
                                               Greenleaf_RNA = Greenleaf_RNA_Seurat@assays$converted@counts,
                                               mouse_multiome = mouse_multiome_Seurat@assays$converted@counts))

Brain_RNA_liger <- rliger::normalize(object = Brain_RNA_liger,verbose = TRUE)
Brain_RNA_liger <- rliger::selectGenes(object = Brain_RNA_liger,do.plot = TRUE)
Brain_RNA_liger <- rliger::scaleNotCenter(object = Brain_RNA_liger,verbose = TRUE)
gc()

#iNMF
Brain_RNA_liger <- rliger::optimizeALS(object = Brain_RNA_liger,k = 30)
Brain_RNA_liger <- quantile_norm(Brain_RNA_liger)
Brain_RNA_liger <- runUMAP(object = Brain_RNA_liger,distance = 'cosine',n_neighbors = 50,min_dist = 0.6)
all.plots <- plotByDatasetAndCluster(Brain_RNA_liger,axis.labels = c('UMAP 1', 'UMAP 2'),return.plots = TRUE)
all.plots[[1]] + all.plots[[2]]

#converte to seurat
Brain_RNA_Seurat <- ligerToSeurat(object = Brain_RNA_liger,nms = NULL)
Brain_RNA_Seurat$raw_cell_type <- NA
Brain_RNA_Seurat@meta.data[colnames(macaque_multiome_Seurat),"raw_cell_type"] <- macaque_multiome_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"raw_cell_type"] <- Greenleaf_RNA_Seurat$cell_type
Brain_RNA_Seurat@meta.data[colnames(mouse_multiome_Seurat),"raw_cell_type"] <- mouse_multiome_Seurat$cell_type
Brain_RNA_Seurat$species <- NA
Brain_RNA_Seurat@meta.data[colnames(macaque_multiome_Seurat),"species"] <- 'macaque'
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"species"] <- 'human'
Brain_RNA_Seurat@meta.data[colnames(mouse_multiome_Seurat),"species"] <- 'mouse'

DimPlot(Brain_RNA_Seurat,group.by = 'raw_cell_type',split.by = 'species',label = TRUE,repel = TRUE)

#at some points, it works!

# try different liger integration parameters ------------------------------

## select variable features without mouse --------------------------------------
Brain_RNA_liger <- rliger::selectGenes(object = Brain_RNA_liger,do.plot = FALSE,datasets.use = 1:2)