#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: filter Cycling cells with In fate in macaque snRNA data         ##
## Data: 2022.08.03                                                                ##
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

# filter In lineage Cycling -----------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220718_summary/macaque_RNA_Seurat_220718.rds')

p1 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'RNA_snn_res.0.5',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#further cluster
macaque_RNA_Seurat <- FindClusters(object = macaque_RNA_Seurat,resolution = 2.5)
p1 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'RNA_snn_res.2.5',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

table(macaque_RNA_Seurat$RNA_snn_res.2.5,macaque_RNA_Seurat$cell_type)

DimPlot(object = macaque_RNA_Seurat,group.by = 'RNA_snn_res.2.5',label = TRUE,repel = TRUE)
FeaturePlot(object = macaque_RNA_Seurat,features = c('GAD1','GAD2','DLX2','DLX5'),order = TRUE)
VlnPlot(object = macaque_RNA_Seurat,features = c('GAD1','GAD2','DLX2','DLX5'),pt.size = 0,assay = 'RNA',group.by = 'RNA_snn_res.2.5',slot = 'data')
#cluster 22 can be filted

cell_list <- colnames(macaque_RNA_Seurat)[!(macaque_RNA_Seurat$RNA_snn_res.2.5 %in% c('22'))]
#saveRDS(cell_list,file = './res/step_53_fig_220803/macaque_RNA_Seurat_cell_list.rds')

macaque_RNA_Seurat <- macaque_RNA_Seurat[,cell_list]
DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
DimPlot(object = macaque_RNA_Seurat,group.by = 'RNA_snn_res.2.5',label = TRUE,repel = TRUE)
FeaturePlot(object = macaque_RNA_Seurat,features = c('GAD1','GAD2','DLX2','DLX5'),order = TRUE)
#seems ok

# label transfer from macaque multiome Seurat -----------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
macaque_RNA_Seurat <- readRDS(file = './processed_data/220718_summary/macaque_RNA_Seurat_220718.rds')
cell_list <- readRDS(file = './res/step_53_fig_220803/macaque_RNA_Seurat_cell_list.rds')

macaque_RNA_Seurat <- macaque_RNA_Seurat[,cell_list]

#recreate macaque RNA Seurat
meta_data <- macaque_RNA_Seurat@meta.data
meta_data <- meta_data[,c("batch","donor","tech","cell_type","sub_cell_type")]
macaque_RNA_Seurat <- macaque_RNA_Seurat@assays$RNA@counts
macaque_RNA_Seurat <- CreateSeuratObject(counts = macaque_RNA_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#re-process macaque RNA Seurat
macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = VariableFeatures(macaque_multiome_Seurat),vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
temp <- projectMatrix_SeuratUMAP(X_scaled = macaque_RNA_Seurat@assays$RNA@scale.data,object = macaque_multiome_Seurat,assayUsed = 'RNA',missing_gene = FALSE)
macaque_RNA_Seurat@reductions$projected_PCA <- CreateDimReducObject(embeddings = temp$pcaCoord_proj,assay = 'RNA')
macaque_RNA_Seurat@reductions$projected_UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')
DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE)

#label transfer
anchors <- my_FindTransferAnchors(reference = macaque_multiome_Seurat,query = macaque_RNA_Seurat,ref_reduction = 'pca',query_reduction = 'projected_PCA',ref_assay = 'RNA',query_assay = 'RNA',l2.norm = FALSE,dims = 1:31,features = VariableFeatures(macaque_multiome_Seurat),verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$cell_type,l2.norm = FALSE,dims = 1:31,eps = TRUE,verbose = TRUE)
macaque_RNA_Seurat <- AddMetaData(object = macaque_RNA_Seurat,metadata = predictions)
macaque_RNA_Seurat$predicted_cell_type <- macaque_RNA_Seurat$predicted.id
p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'predicted_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

predictions <- TransferData(anchorset = anchors,refdata = macaque_multiome_Seurat$sub_cell_type,l2.norm = FALSE,dims = 1:31,eps = TRUE,verbose = TRUE)
macaque_RNA_Seurat <- AddMetaData(object = macaque_RNA_Seurat,metadata = predictions)
macaque_RNA_Seurat$predicted_sub_cell_type <- macaque_RNA_Seurat$predicted.id
p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'projected_UMAP',group.by = 'predicted_sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

table(macaque_RNA_Seurat$predicted_cell_type,macaque_RNA_Seurat$predicted_sub_cell_type)
scibet::Confusion_heatmap(ori = macaque_RNA_Seurat$cell_type,prd = macaque_RNA_Seurat$predicted_cell_type)
scibet::Confusion_heatmap(ori = macaque_RNA_Seurat$cell_type,prd = macaque_RNA_Seurat$predicted_sub_cell_type)

#seems perform well

#save meta_data
meta_data <- macaque_RNA_Seurat@meta.data
saveRDS(object = meta_data,file = './res/step_53_fig_220803/macaque_RNA_Seurat_predicted_meta_data.rds')

# generate final macaque RNA Seurat ---------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220718_summary/macaque_RNA_Seurat_220718.rds')
meta_data <- readRDS(file = './res/step_53_fig_220803/macaque_RNA_Seurat_predicted_meta_data.rds')
macaque_RNA_Seurat <- macaque_RNA_Seurat[,rownames(meta_data)]

#re-create macaque RNA Seurat
macaque_RNA_Seurat <- macaque_RNA_Seurat@assays$RNA@counts
macaque_RNA_Seurat <- CreateSeuratObject(counts = macaque_RNA_Seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
macaque_RNA_Seurat$batch <- meta_data[colnames(macaque_RNA_Seurat),"batch"]
macaque_RNA_Seurat$donor <- meta_data[colnames(macaque_RNA_Seurat),"donor"]
macaque_RNA_Seurat$tech <- meta_data[colnames(macaque_RNA_Seurat),"tech"]
macaque_RNA_Seurat$cell_type <- meta_data[colnames(macaque_RNA_Seurat),"predicted_cell_type"]
macaque_RNA_Seurat$sub_cell_type <- meta_data[colnames(macaque_RNA_Seurat),"predicted_sub_cell_type"]

#re-process macaque RNA Seurat
macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
ndims <- 28
macaque_RNA_Seurat <- FindNeighbors(object = macaque_RNA_Seurat,reduction = 'pca',dims = 1:ndims)
macaque_RNA_Seurat <- FindClusters(object = macaque_RNA_Seurat,resolution = 0.7)
macaque_RNA_Seurat <- RunUMAP(object = macaque_RNA_Seurat,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#refine annotation
#cluster level
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.0.7 %in% c('3','4'),"sub_cell_type"] <- 'InCGE'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.0.7 %in% c('2','6','12','21'),"sub_cell_type"] <- 'InMGE'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.0.7 %in% c('15'),"sub_cell_type"] <- 'OPC'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.0.7 %in% c('19'),"sub_cell_type"] <- 'Mic'

#cell type level
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Cycling'),"cell_type"] <- 'Cycling'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('End'),"cell_type"] <- 'End'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Ex-1'),"cell_type"] <- 'Ex-1'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Ex-2'),"cell_type"] <- 'Ex-2'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Ex-3_1','Ex-3_2'),"cell_type"] <- 'Ex-3'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Ex-4_1','Ex-4_2'),"cell_type"] <- 'Ex-4'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('InCGE'),"cell_type"] <- 'InCGE'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('InMGE'),"cell_type"] <- 'InMGE'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('IP'),"cell_type"] <- 'IP'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Mic'),"cell_type"] <- 'Mic'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('OPC'),"cell_type"] <- 'OPC'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('Per'),"cell_type"] <- 'Per'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('RG-1'),"cell_type"] <- 'RG-1'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('RG-2_1','RG-2_2'),"cell_type"] <- 'RG-2'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$sub_cell_type %in% c('VLMC'),"cell_type"] <- 'VLMC'

#check data
table(macaque_RNA_Seurat$sub_cell_type,macaque_RNA_Seurat$cell_type)
table(macaque_RNA_Seurat$RNA_snn_res.0.7,macaque_RNA_Seurat$cell_type)
p1 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = macaque_RNA_Seurat,file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')

# generate macaque RNA Seurat human symbol --------------------------------
# load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')

#convert to human symbol
macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
macaque_to_human_anno <- macaque_to_human_anno[,c(2,1,4,3)]
express_matrix <- macaque_RNA_Seurat@assays$RNA@counts
express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = express_matrix,anno = macaque_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)

#create macaque RNA Seurat human symbol
meta_data <- macaque_RNA_Seurat@meta.data
meta_data <- meta_data[,c("batch","donor","tech","cell_type","sub_cell_type")]
macaque_RNA_Seurat_human_symbol <- CreateSeuratObject(counts = express_matrix,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#process of macaque RNA Seurat human symbol
macaque_RNA_Seurat_human_symbol <- my_process_seurat(object = macaque_RNA_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
ndims <- 27
macaque_RNA_Seurat_human_symbol <- FindNeighbors(object = macaque_RNA_Seurat_human_symbol,reduction = 'pca',dims = 1:ndims)
macaque_RNA_Seurat_human_symbol <- FindClusters(object = macaque_RNA_Seurat_human_symbol,resolution = 0.7)
macaque_RNA_Seurat_human_symbol <- RunUMAP(object = macaque_RNA_Seurat_human_symbol,dims = 1:ndims,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
p1 <- DimPlot(object = macaque_RNA_Seurat_human_symbol,reduction = 'umap',group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_RNA_Seurat_human_symbol,reduction = 'umap',group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion2,ArchRPalettes$bear)))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = macaque_RNA_Seurat_human_symbol,file = './processed_data/220802_summary/macaque_RNA_Seurat_human_symbol_220803.rds')
