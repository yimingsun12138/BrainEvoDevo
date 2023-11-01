#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque RNA data convert gene symbol to human symbol            ##
## Data: 2022.07.18                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# convert macaque multiome Seurat -----------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')

#convert gene id
macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
macaque_to_human_anno <- macaque_to_human_anno[,c(2,1,4,3)]

temp <- macaque_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^4),workers = 6)

meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech","donor","cell_type","sub_cell_type")]

#re-create macaque multiome data
macaque_multiome_Seurat_human_symbol <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#process data
ndims <- 33
macaque_multiome_Seurat_human_symbol <- my_process_seurat(object = macaque_multiome_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat_human_symbol <- FindNeighbors(object = macaque_multiome_Seurat_human_symbol,reduction = 'pca',dims = 1:ndims,assay = 'RNA',verbose = TRUE)
macaque_multiome_Seurat_human_symbol <- FindClusters(object = macaque_multiome_Seurat_human_symbol,resolution = seq(from = 0.3,to = 2.0,by = 0.1),verbose = TRUE)
macaque_multiome_Seurat_human_symbol <- RunUMAP(object = macaque_multiome_Seurat_human_symbol,dims = 1:ndims,reduction = 'pca',assay = 'RNA',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_multiome_Seurat_human_symbol,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = macaque_multiome_Seurat_human_symbol,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = macaque_multiome_Seurat_human_symbol,file = './processed_data/220718_summary/macaque_multiome_Seurat_human_symbol_220718.rds')

# convert macaque RNA Seurat ----------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220718_summary/macaque_RNA_Seurat_220718.rds')

#convert gene name
macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
macaque_to_human_anno <- macaque_to_human_anno[,c(2,1,4,3)]

temp <- macaque_RNA_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 4)

meta_data <- macaque_RNA_Seurat@meta.data
meta_data <- meta_data[,c("batch","donor","tech","cell_type","sub_cell_type")]

#re-create macaque RNA data
macaque_RNA_Seurat_human_symbol <- CreateSeuratObject(counts = temp,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#process data
ndims <- 31
macaque_RNA_Seurat_human_symbol <- my_process_seurat(object = macaque_RNA_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
macaque_RNA_Seurat_human_symbol <- FindNeighbors(object = macaque_RNA_Seurat_human_symbol,reduction = 'pca',dims = 1:ndims,assay = 'RNA',verbose = TRUE)
macaque_RNA_Seurat_human_symbol <- FindClusters(object = macaque_RNA_Seurat_human_symbol,resolution = seq(from = 0.3,to = 2.0,by = 0.1),verbose = TRUE)
macaque_RNA_Seurat_human_symbol <- RunUMAP(object = macaque_RNA_Seurat_human_symbol,dims = 1:ndims,reduction = 'pca',assay = 'RNA',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_RNA_Seurat_human_symbol,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = macaque_RNA_Seurat_human_symbol,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = macaque_RNA_Seurat_human_symbol,file = './processed_data/220718_summary/macaque_RNA_Seurat_human_symbol_220719.rds')
