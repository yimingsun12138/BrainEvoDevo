#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: seurat CCA integrate species Brain RNA data                     ##
## Data: 2022.07.19                                                                ##
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

# Greenleaf data process --------------------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/220711_summary/Greenleaf_RNA_Seurat_220714.rds')
temp <- Greenleaf_RNA_Seurat@assays$converted@counts
meta_data <- Greenleaf_RNA_Seurat@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_originalexp','nFeature_originalexp','nCount_converted','nFeature_converted','macaque_cell_type','macaque_cell_type_score','macaque_sub_cell_type','macaque_sub_cell_type_score'))]

#re-create seurat object
ndims <- 45
Greenleaf_RNA_Seurat <- CreateSeuratObject(counts = temp,project = 'human',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
Greenleaf_RNA_Seurat <- my_process_seurat(object = Greenleaf_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
Greenleaf_RNA_Seurat <- my_process_seurat(object = Greenleaf_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = ndims,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5),group.by = 'cell_type',label = TRUE)
p1 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'RNA_snn_res.0.3',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#save temp data
saveRDS(object = Greenleaf_RNA_Seurat,file = '/data/User/sunym/temp/Greenleaf_RNA_Seurat_temp.rds')

# mouse data process ------------------------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_meta_data.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]
mouse_multiome_Seurat$cell_type <- meta_data$cell_type

#convert gene name
mouse_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human_anno <- mouse_to_human_anno[,c(2,1,4,3)]
temp <- mouse_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)

meta_data <- mouse_multiome_Seurat@meta.data
meta_data <- meta_data[,c("percent.mt","donor","Age","cell_type")]

#re-create mouse multiome Seurat
mouse_multiome_Seurat_human_symbol <- CreateSeuratObject(counts = temp,project = 'mouse',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
mouse_multiome_Seurat_human_symbol <- my_process_seurat(object = mouse_multiome_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat_human_symbol <- my_process_seurat(object = mouse_multiome_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 16,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5),group.by = 'cell_type',label = TRUE)
p1 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'RNA_snn_res.0.3',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#save temp data
saveRDS(object = mouse_multiome_Seurat_human_symbol,file = '/data/User/sunym/temp/mouse_multiome_Seurat_temp.rds')

# run seurat cca ----------------------------------------------------------
#load data
macaque_multiome_Seurat_human_symbol <- readRDS(file = './processed_data/220718_summary/macaque_multiome_Seurat_human_symbol_220718.rds')
Greenleaf_RNA_Seurat <- readRDS(file = '/data/User/sunym/temp/Greenleaf_RNA_Seurat_temp.rds')
mouse_multiome_Seurat_human_symbol <- readRDS(file = '/data/User/sunym/temp/mouse_multiome_Seurat_temp.rds')

#select integration features
gene_list <- SelectIntegrationFeatures(object.list = list(Greenleaf_RNA_Seurat,mouse_multiome_Seurat_human_symbol,macaque_multiome_Seurat_human_symbol),verbose = TRUE)

#find integration anchor
anchors <- FindIntegrationAnchors(object.list = list(Greenleaf_RNA_Seurat,mouse_multiome_Seurat_human_symbol,macaque_multiome_Seurat_human_symbol),
                                  anchor.features = gene_list,reduction = 'cca',dims = 1:50,verbose = TRUE)

#integrate data
Brain_RNA_Seurat <- IntegrateData(anchorset = anchors,dims = 1:50,verbose = TRUE)

#process
Brain_RNA_Seurat$species <- NA
Brain_RNA_Seurat@meta.data[colnames(macaque_multiome_Seurat_human_symbol),"species"] <- 'macaque'
Brain_RNA_Seurat@meta.data[colnames(Greenleaf_RNA_Seurat),"species"] <- 'human'
Brain_RNA_Seurat@meta.data[colnames(mouse_multiome_Seurat_human_symbol),"species"] <- 'mouse'

DefaultAssay(Brain_RNA_Seurat) <- 'integrated'
Brain_RNA_Seurat <- ScaleData(object = Brain_RNA_Seurat,features = rownames(Brain_RNA_Seurat@assays$integrated@data),assay = 'integrated',vars.to.regress = NULL,verbose = TRUE)
Brain_RNA_Seurat <- RunPCA(object = Brain_RNA_Seurat,assay = 'integrated',features = rownames(Brain_RNA_Seurat@assays$integrated@data),npcs = 50,verbose = TRUE)
Brain_RNA_Seurat <- FindNeighbors(object = Brain_RNA_Seurat,reduction = 'pca',dims = 1:29,assay = 'integrated',verbose = TRUE)
Brain_RNA_Seurat <- FindClusters(object = Brain_RNA_Seurat,resolution = c(0.3,0.5,0.7,0.9,1.1,1.3,1.5),verbose = TRUE)
Brain_RNA_Seurat <- RunUMAP(object = Brain_RNA_Seurat,dims = 1:29,reduction = 'pca',assay = 'integrated',verbose = TRUE)

#plot
DimPlot(object = Brain_RNA_Seurat,group.by = 'species',label = FALSE)
DimPlot(object = Brain_RNA_Seurat,group.by = 'cell_type',split.by = 'species',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))

#save data
saveRDS(object = Brain_RNA_Seurat,file = './processed_data/220718_summary/Greenleaf_RNA_Seurat_AND_mouse_multiome_Seurat_AND_macaque_multiome_Seurat_CCA_integrated_Seurat_220721.rds')
