#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: filter Cycling cells with In fate in macaque multiome RNA       ##
## Data: 2022.08.01                                                                ##
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
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/220718_summary/macaque_multiome_ArchR_220720/')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')

p1 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),group.by = 'Clusters',label = TRUE,repel = TRUE)
p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#filter C7 cycling
C7_cell_list <- rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters == 'C7']
#saveRDS(object = C7_cell_list,file = './res/step_52_fig_220802/In_fate_Cycling_list.rds')

# re-create macaque multiome data -----------------------------------------
macaque_multiome_old <- macaque_multiome_Seurat
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(colnames(macaque_multiome_Seurat) %in% C7_cell_list)]
DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech","donor","cell_type","sub_cell_type")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

# re-process macaque multiome data ----------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
ndims <- 31
macaque_multiome_Seurat <- FindNeighbors(object = macaque_multiome_Seurat,reduction = 'pca',dims = 1:ndims)
macaque_multiome_Seurat <- FindClusters(object = macaque_multiome_Seurat,resolution = 0.7)
macaque_multiome_Seurat <- RunUMAP(object = macaque_multiome_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#refine annotation

#Ex-2
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('6'),"cell_type"] <- 'Ex-2'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('6'),"sub_cell_type"] <- 'Ex-2'

#Ex-3
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('0','2','7'),"cell_type"] <- 'Ex-3'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('0','2'),"sub_cell_type"] <- 'Ex-3_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('7'),"sub_cell_type"] <- 'Ex-3_2'

#Ex-4
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('9','14'),"cell_type"] <- 'Ex-4'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('9'),"sub_cell_type"] <- 'Ex-4_1'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('14'),"sub_cell_type"] <- 'Ex-4_2'

#InMGE
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('1','8','10'),"cell_type"] <- 'InMGE'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('1','8','10'),"sub_cell_type"] <- 'InMGE'

#InCGE
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('3','4'),"cell_type"] <- 'InCGE'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('3','4'),"sub_cell_type"] <- 'InCGE'

#Mic
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('17'),"cell_type"] <- 'Mic'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('17'),"sub_cell_type"] <- 'Mic'

#End
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('15'),"cell_type"] <- 'End'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('15'),"sub_cell_type"] <- 'End'

#Per
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('18'),"cell_type"] <- 'Per'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('18'),"sub_cell_type"] <- 'Per'

#VLMC
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('19'),"cell_type"] <- 'VLMC'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('19'),"sub_cell_type"] <- 'VLMC'

#OPC
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('12'),"cell_type"] <- 'OPC'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('12'),"sub_cell_type"] <- 'OPC'

#Cycling
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type %in% c('Cycling'),"cell_type"] <- 'Cycling'
macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type %in% c('Cycling'),"sub_cell_type"] <- 'Cycling'

table(macaque_multiome_Seurat$cell_type,macaque_multiome_Seurat$sub_cell_type)
p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)
p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#cells in cluster 5
cell_list <- colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('5') & !(macaque_multiome_Seurat$cell_type %in% c('IP','Ex-1'))]
table(macaque_multiome_Seurat@meta.data[cell_list,"cell_type"])
reference_seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$RNA_snn_res.0.7 %in% c('5') & !(colnames(macaque_multiome_Seurat) %in% cell_list)]
query_seurat <- macaque_multiome_Seurat[,cell_list]
anchors <- my_FindTransferAnchors(reference = reference_seurat,query = query_seurat,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'RNA',query_assay = 'RNA',l2.norm = FALSE,dims = 1:31,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = reference_seurat$cell_type,l2.norm = FALSE,dims = 1:31,verbose = TRUE)
query_seurat <- AddMetaData(object = query_seurat,metadata = predictions)
p1 <- DimPlot(object = reference_seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = query_seurat,reduction = 'umap',group.by = 'predicted.id',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

macaque_multiome_Seurat@meta.data[colnames(query_seurat),"cell_type"] <- query_seurat$predicted.id
macaque_multiome_Seurat@meta.data[colnames(query_seurat),"sub_cell_type"] <- query_seurat$predicted.id

p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)
p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = macaque_multiome_Seurat,file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')

#check data
table(macaque_multiome_Seurat$cell_type,macaque_multiome_Seurat$sub_cell_type)
table(macaque_multiome_Seurat$RNA_snn_res.0.7,macaque_multiome_Seurat$cell_type)

# convert to human symbol -------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')

#convert to human symbol
macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
macaque_to_human_anno <- macaque_to_human_anno[,c(2,1,4,3)]
express_matrix <- macaque_multiome_Seurat@assays$RNA@counts
express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = express_matrix,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)

#create macaque multiome seurat human symbol
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("batch","tech","donor","cell_type","sub_cell_type")]
macaque_multiome_Seurat_human_symbol <- CreateSeuratObject(counts = express_matrix,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#process macaque multiome seurat human symbol
macaque_multiome_Seurat_human_symbol <- my_process_seurat(object = macaque_multiome_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
ndims <- 32
macaque_multiome_Seurat_human_symbol <- FindNeighbors(object = macaque_multiome_Seurat_human_symbol,reduction = 'pca',dims = 1:ndims)
macaque_multiome_Seurat_human_symbol <- FindClusters(object = macaque_multiome_Seurat_human_symbol,resolution = 0.7)
macaque_multiome_Seurat_human_symbol <- RunUMAP(object = macaque_multiome_Seurat_human_symbol,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_multiome_Seurat_human_symbol,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat_human_symbol,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

table(macaque_multiome_Seurat_human_symbol$seurat_clusters,macaque_multiome_Seurat_human_symbol$cell_type)

#save data
saveRDS(object = macaque_multiome_Seurat_human_symbol,file = './processed_data/220802_summary/macaque_multiome_Seurat_human_symbol_220802.rds')
