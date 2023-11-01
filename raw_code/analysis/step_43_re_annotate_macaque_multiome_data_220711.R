#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: re_annotate macaque multiome data                               ##
## Data: 2022.07.11                                                                ##
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

# load macaque multiome data ----------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']

#recreate macaque multiome data
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech",'donor',"cell_type")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
gc()

# process of macaque multiome data ----------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,
                                             nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)

ndims <- 31
macaque_multiome_Seurat <- RunUMAP(object = macaque_multiome_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
macaque_multiome_Seurat <- FindNeighbors(object = macaque_multiome_Seurat,reduction = 'pca',dims = 1:ndims)
macaque_multiome_Seurat <- FindClusters(object = macaque_multiome_Seurat,resolution = c(0.7))
DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))

#show donor contribution
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'donor',label = FALSE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',split.by = 'donor',ncol = 2,label = FALSE) + theme(aspect.ratio = 1)

# annotation --------------------------------------------------------------
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_Seurat_reannotate_round_1_meta_data.rds')

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#sub-cluster on cluster 7
Idents(macaque_multiome_Seurat) <- 'RNA_snn_res.0.7'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = '7',graph.name = 'RNA_snn',resolution = 0.15)

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'sub.cluster',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#annotate
macaque_multiome_Seurat$cell_type <- NA
#OPC
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('13'),"cell_type"] <- 'OPC'
#RG-2
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('12'),"cell_type"] <- 'RG-2'
#RG-1
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('14'),"cell_type"] <- 'RG-1'
#Cycling
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('11'),"cell_type"] <- 'Cycling'
#IP
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('7_1'),"cell_type"] <- 'IP'
#Ex-1
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('7_0'),"cell_type"] <- 'Ex-1'
#Ex-2
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('5'),"cell_type"] <- 'Ex-2'
#Ex-3
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('0','3','6'),"cell_type"] <- 'Ex-3'
#Ex-4
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('10','15'),"cell_type"] <- 'Ex-4'
#InCGE
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('2','4'),"cell_type"] <- 'InCGE'
#InMGE
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('1','8','9'),"cell_type"] <- 'InMGE'
#Mic
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('17'),"cell_type"] <- 'Mic'
#End
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('16'),"cell_type"] <- 'End'
#Per
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('18'),"cell_type"] <- 'Per'
#VLMC
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster %in% c('19'),"cell_type"] <- 'VLMC'

table(macaque_multiome_Seurat$cell_type)
DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#sub-cluster
macaque_multiome_Seurat$sub_cell_type <- macaque_multiome_Seurat$cell_type
#Cycling
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'Cycling',graph.name = 'RNA_snn',resolution = 0.1)
FeaturePlot(object = macaque_multiome_Seurat,features = c('TOP2A','MKI67','CLSPN','AURKA'),pt.size = 0.1)
DimPlot(object = macaque_multiome_Seurat,group.by = 'sub.cluster',label = TRUE,repel = TRUE)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'Cycling_2',"cell_type"] <- 'IP'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'Cycling_2',"sub_cell_type"] <- 'IP'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'Cycling_0',"sub_cell_type"] <- 'Cyc-G2M'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'Cycling_1',"sub_cell_type"] <- 'Cyc-S'

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#RG-2
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'RG-2',graph.name = 'RNA_snn',resolution = 0.08)
DimPlot(object = macaque_multiome_Seurat,group.by = 'sub.cluster',label = TRUE,repel = TRUE)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'RG-2_0',"sub_cell_type"] <- 'RG-2_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'RG-2_1',"sub_cell_type"] <- 'RG-2_2'

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#Ex-3
p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('0','3'),"sub_cell_type"] <- 'Ex-3_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('6'),"sub_cell_type"] <- 'Ex-3_2'

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#Ex-4
p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('10'),"sub_cell_type"] <- 'Ex-4_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('15'),"sub_cell_type"] <- 'Ex-4_2'

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

p1 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = macaque_multiome_Seurat@reductions$umap@cell.embeddings,meta_data = meta_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

macaque_multiome_Seurat$sub.cluster <- macaque_multiome_Seurat$sub_cell_type

p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = macaque_multiome_Seurat,file = './processed_data/220711_summary/macaque_multiome_Seurat_220712.rds')