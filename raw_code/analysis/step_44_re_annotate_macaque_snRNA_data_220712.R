#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: re_annotate macaque snRNA data                                  ##
## Data: 2022.07.12                                                                ##
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

# load data ---------------------------------------------------------------
macaque_RNA_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_RNA_Seurat <- macaque_RNA_Seurat[,macaque_RNA_Seurat$tech == 'scRNA']
macaque_multiome_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_multiome_Seurat_220712.rds')

meta_data <- macaque_RNA_Seurat@meta.data
meta_data <- meta_data[,c("batch","donor","tech","cell_type")]
macaque_RNA_Seurat <- macaque_RNA_Seurat@assays$RNA@counts
macaque_RNA_Seurat <- CreateSeuratObject(counts = macaque_RNA_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

# label transfer of multiome data -----------------------------------------
macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = VariableFeatures(macaque_multiome_Seurat),vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
temp <- projectMatrix_SeuratUMAP(X_scaled = macaque_RNA_Seurat@assays$RNA@scale.data,object = macaque_multiome_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
macaque_RNA_Seurat@reductions$projected_PCA <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
macaque_RNA_Seurat@reductions$projected_UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')
DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE)

anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = macaque_RNA_Seurat,ref_reduction = 'pca',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'RNA',l2.norm = TRUE,dims = 1:31,features = VariableFeatures(macaque_multiome_Seurat),verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$cell_type,l2.norm = TRUE,dims = 1:31,eps = TRUE,verbose = TRUE)
macaque_RNA_Seurat <- AddMetaData(object = macaque_RNA_Seurat,metadata = predictions)
macaque_RNA_Seurat$predicted_cell_type <- macaque_RNA_Seurat$predicted.id
p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'predicted_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$sub_cell_type,l2.norm = TRUE,dims = 1:31,eps = TRUE,verbose = TRUE)
macaque_RNA_Seurat <- AddMetaData(object = macaque_RNA_Seurat,metadata = predictions)
macaque_RNA_Seurat$predicted_sub_cell_type <- macaque_RNA_Seurat$predicted.id
p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'predicted_sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

scibet::Confusion_heatmap(ori = macaque_RNA_Seurat$cell_type,prd = macaque_RNA_Seurat$predicted_cell_type)
scibet::Confusion_heatmap(ori = macaque_RNA_Seurat$cell_type,prd = macaque_RNA_Seurat$predicted_sub_cell_type)

#It performs well

# generate the final macaque RNA Seurat -----------------------------------
meta_data <- macaque_RNA_Seurat@meta.data
macaque_RNA_Seurat <- macaque_RNA_Seurat@assays$RNA@counts
macaque_RNA_Seurat <- CreateSeuratObject(counts = macaque_RNA_Seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
macaque_RNA_Seurat$batch <- meta_data[colnames(macaque_RNA_Seurat),"batch"]
macaque_RNA_Seurat$donor <- meta_data[colnames(macaque_RNA_Seurat),"donor"]
macaque_RNA_Seurat$tech <- meta_data[colnames(macaque_RNA_Seurat),"tech"]
macaque_RNA_Seurat$cell_type <- meta_data[colnames(macaque_RNA_Seurat),"predicted_cell_type"]
macaque_RNA_Seurat$sub_cell_type <- meta_data[colnames(macaque_RNA_Seurat),"predicted_sub_cell_type"]

#little modify
macaque_RNA_Seurat$cell_type <- macaque_RNA_Seurat$sub_cell_type
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Cyc-S','Cyc-G2M'),"cell_type"] <- 'Cycling'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Ex-3_1','Ex-3_2'),"cell_type"] <- 'Ex-3'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Ex-4_1','Ex-4_2'),"cell_type"] <- 'Ex-4'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('RG-2_1','RG-2_2'),"cell_type"] <- 'RG-2'
scibet::Confusion_heatmap(ori = macaque_RNA_Seurat$cell_type,prd = macaque_RNA_Seurat$sub_cell_type)

macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
macaque_RNA_Seurat <- FindNeighbors(object = macaque_RNA_Seurat,dims = 1:27,assay = 'RNA')
macaque_RNA_Seurat <- FindClusters(object = macaque_RNA_Seurat,resolution = 0.5)
macaque_RNA_Seurat <- RunUMAP(object = macaque_RNA_Seurat,dims = 1:27,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#re-subcluster Cycling
Idents(macaque_RNA_Seurat) <- 'cell_type'
macaque_RNA_Seurat <- FindSubCluster(object = macaque_RNA_Seurat,cluster = 'Cycling',graph.name = 'RNA_snn',resolution = 0.1)
p1 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'sub.cluster',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)
FeaturePlot(object = macaque_RNA_Seurat,features = c('CLSPN','AURKA'))

macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub.cluster == 'Cycling_0',"sub_cell_type"] <- 'Cyc-S'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub.cluster == 'Cycling_1',"sub_cell_type"] <- 'Cyc-G2M'
p1 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

macaque_RNA_Seurat$sub.cluster <- macaque_RNA_Seurat$sub_cell_type

#save data
saveRDS(object = macaque_RNA_Seurat,file = './processed_data/220711_summary/macaque_RNA_Seurat_220712.rds')
