#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: generate processed mouse multiome RNA data                      ##
## Data: 2022.07.22                                                                ##
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
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_meta_data.rds')
Brain_RNA_Seurat <- readRDS(file = './processed_data/220718_summary/Greenleaf_RNA_Seurat_AND_mouse_multiome_Seurat_AND_macaque_multiome_Seurat_CCA_integrated_Seurat_220721.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/220718_summary/macaque_multiome_Seurat_220718.rds')

mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]
mouse_multiome_Seurat$cell_type <- meta_data$cell_type

# label transfer ----------------------------------------------------------
macaque_data <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'macaque']
mouse_data <- Brain_RNA_Seurat[,Brain_RNA_Seurat$species == 'mouse']

#add annotation
macaque_data$cell_type <- macaque_multiome_Seurat@meta.data[colnames(macaque_data),"cell_type"]
macaque_data$sub_cell_type <- macaque_multiome_Seurat@meta.data[colnames(macaque_data),"sub_cell_type"]
macaque_data <- macaque_data[,!(macaque_data$cell_type %in% c('Mic','OPC'))]

#label transfer
#cell type
anchors <- my_FindTransferAnchors(reference = macaque_data,query = mouse_data,ref_reduction = 'pca',query_reduction = 'pca',ref_assay = 'integrated',query_assay = 'integrated',l2.norm = FALSE,dims = 1:29,verbose = TRUE)
predictions <- TransferData(anchorset = anchors,refdata = macaque_data$cell_type,l2.norm = FALSE,dims = 1:29,verbose = TRUE)
mouse_data <- AddMetaData(object = mouse_data,metadata = predictions)
DimPlot(mouse_data,group.by = 'predicted.id',label = TRUE,repel = TRUE)

mouse_multiome_Seurat$macaque_cell_type <- mouse_data@meta.data[colnames(mouse_multiome_Seurat),"predicted.id"]
mouse_multiome_Seurat$macaque_cell_type_score <- mouse_data@meta.data[colnames(mouse_multiome_Seurat),"prediction.score.max"]

#sub cell type
predictions <- TransferData(anchorset = anchors,refdata = macaque_data$sub_cell_type,l2.norm = FALSE,dims = 1:29,verbose = TRUE)
mouse_data <- AddMetaData(object = mouse_data,metadata = predictions)
DimPlot(mouse_data,group.by = 'predicted.id',label = TRUE,repel = TRUE)

mouse_multiome_Seurat$macaque_sub_cell_type <- mouse_data@meta.data[colnames(mouse_multiome_Seurat),"predicted.id"]
mouse_multiome_Seurat$macaque_sub_cell_type_score <- mouse_data@meta.data[colnames(mouse_multiome_Seurat),"prediction.score.max"]

# re-process mouse multiome data ------------------------------------------
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 27,resolution = 0.5,group.by = 'macaque_cell_type',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1+p2+plot_layout(ncol = 2)

table(mouse_multiome_Seurat$macaque_cell_type,mouse_multiome_Seurat$seurat_clusters)

#marker express
FeaturePlot(object = mouse_multiome_Seurat,features = c('Slc17a6','Dpy19l1','Nrp1'))
FeaturePlot(object = mouse_multiome_Seurat,features = c('Kcnk2','Hmcn1','Prss12','Ccbe1'))
FeaturePlot(object = mouse_multiome_Seurat,features = c('Rxfp1','Pdzd2','Trpm3','Tle4','Hs3st4'))
VlnPlot(object = mouse_multiome_Seurat,features = c('Rxfp1','Pdzd2','Trpm3','Tle4','Hs3st4'),group.by = 'seurat_clusters')

# find specific marker of Ex-1 and Ex-2 -----------------------------------
#Ex-1
Ex_1_marker <- my_seurat_find_specific_marker(seu.obj = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4')],assay = 'RNA',ident.1 = 'Ex-1',group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.5,pct.1 = 0.3,overlap_ratio = 0.85,workers = 6,future.globals.maxSize = 100*(1024^3))
VlnPlot(object = macaque_multiome_Seurat,features = Ex_1_marker[85:92],group.by = 'cell_type',pt.size = 0)
#PALMD HS3ST1 SLIT2 NDST3 DPY19L1 NRP1 EPHA4 RAB12 SLA
FeaturePlot(object = macaque_multiome_Seurat,features = c('PALMD','HS3ST1','SLIT2','NDST3','DPY19L1','NRP1','EPHA4','RAB12','SLA'))
#PALMD HS3ST1 SLIT2 DPY19L1 NRP1 RAB12 checked
VlnPlot(object = mouse_multiome_Seurat,features = c('Palmd','Hs3st1','Slit2','Dpy19l1','Nrp1','Rab12'),group.by = 'seurat_clusters')

#Ex-2
Ex_2_marker <- my_seurat_find_specific_marker(seu.obj = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4')],assay = 'RNA',ident.1 = 'Ex-2',group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.5,pct.1 = 0.3,overlap_ratio = 0.85,workers = 6,future.globals.maxSize = 100*(1024^3))
VlnPlot(object = macaque_multiome_Seurat,features = Ex_2_marker[97:98],group.by = 'cell_type',pt.size = 0)
#KIF26B HMCN1 KCNK2 MB21D2 BCL6 ITPR1 MYRIP STK32B BMPR1B PRSS12 GALNT10 RNF182 MCUR1 SVEP1 SNTG2 NUP93 CCBE1
FeaturePlot(object = macaque_multiome_Seurat,features = c('KIF26B','HMCN1','KCNK2','MB21D2','BCL6','ITPR1','MYRIP','STK32B','BMPR1B','PRSS12','GALNT10','RNF182','MCUR1','SVEP1','SNTG2','NUP93','CCBE1'))
#KIF26B HMCN1 KCNK2 BCL6 ITPR1 STK32B BMPR1B PRSS12 GALNT10 RNF182 CCBE1 checked
VlnPlot(object = mouse_multiome_Seurat,features = c('Kif26b','Hmcn1','Kcnk2','Bcl6','Itpr1','Stk32b','Bmpr1b','Prss12','Galnt10','Rnf182','Ccbe1'),group.by = 'seurat_clusters')
FeaturePlot(object = mouse_multiome_Seurat,features = c('Kif26b','Hmcn1','Kcnk2','Bcl6','Itpr1','Stk32b','Bmpr1b','Prss12','Galnt10','Rnf182','Ccbe1'))
FeaturePlot(object = mouse_multiome_Seurat,features = c('Ppp1r17','Eomes'))

# re-process mouse multiome data ------------------------------------------
meta_data <- mouse_multiome_Seurat@meta.data
meta_data <- meta_data[,c("percent.mt","donor","Age","macaque_cell_type","macaque_cell_type_score","macaque_sub_cell_type","macaque_sub_cell_type_score")]
mouse_multiome_Seurat <- mouse_multiome_Seurat@assays$RNA@counts
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

#re-process
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 27,resolution = 0.5,group.by = 'macaque_cell_type',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1+p2+plot_layout(ncol = 2)

#re-annotate
mouse_multiome_Seurat$cell_type <- mouse_multiome_Seurat$macaque_cell_type
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('2'),"cell_type"] <- 'InMGE'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('0'),"cell_type"] <- 'Ex-4'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('4'),"cell_type"] <- 'Ex-2'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('1','9','5'),"cell_type"] <- 'Ex-1'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('3'),"cell_type"] <- 'IP'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('10'),"cell_type"] <- 'End'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$RNA_snn_res.0.5 %in% c('11'),"cell_type"] <- 'Per'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$macaque_cell_type %in% c('RG-1','RG-2'),"cell_type"] <- 'RG'

p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1+p2+plot_layout(ncol = 2)

# save data ---------------------------------------------------------------
saveRDS(object = mouse_multiome_Seurat,file = './processed_data/220718_summary/mouse_multiome_Seurat_220723.rds')

# generate mouse multiome data human symbol -------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/220718_summary/mouse_multiome_Seurat_220723.rds')

#symbol convert
mouse_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_to_human_anno <- mouse_to_human_anno[,c(2,1,4,3)]
temp <- mouse_multiome_Seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human_anno,filter_anno = TRUE,future.globals.maxSize = 100*(1024^3),workers = 6)

meta_data <- mouse_multiome_Seurat@meta.data
meta_data <- meta_data[,c("percent.mt","donor","Age","macaque_cell_type","macaque_cell_type_score","macaque_sub_cell_type","macaque_sub_cell_type_score","cell_type")]
mouse_multiome_Seurat_human_symbol <- CreateSeuratObject(counts = temp,project = 'mouse',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)

# process of mouse multiom data human symbol ------------------------------
mouse_multiome_Seurat_human_symbol <- my_process_seurat(object = mouse_multiome_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_multiome_Seurat_human_symbol <- my_process_seurat(object = mouse_multiome_Seurat_human_symbol,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 18,resolution = 0.5,group.by = 'cell_type',label = TRUE)

p1 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = mouse_multiome_Seurat_human_symbol,file = './processed_data/220718_summary/mouse_multiome_Seurat_human_symbol_220723.rds')

# re-generate mouse multiome Seurat ---------------------------------------
#filter some cells do not pass the QC in ATAC processing

## mouse multiome seurat ---------------------------------------------------
mouse_multiome_Seurat <- readRDS(file = './processed_data/220718_summary/mouse_multiome_Seurat_220723.rds')
cell_list <- readRDS(file = './ArchR/res/step_12_fig_220723/cell_list.rds')

#filted cells
temp <- colnames(mouse_multiome_Seurat)[!(colnames(mouse_multiome_Seurat) %in% cell_list)]
table(mouse_multiome_Seurat@meta.data[temp,"cell_type"])
#the filted cells seems on enrichment in certain cell types.

#filter
mouse_multiome_Seurat <- mouse_multiome_Seurat[,cell_list]
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = mouse_multiome_Seurat,file = './processed_data/220718_summary/mouse_multiome_Seurat_220723.rds')

## mouse multiome Seurat human symbol --------------------------------------
#load data
mouse_multiome_Seurat_human_symbol <- readRDS(file = './processed_data/220718_summary/mouse_multiome_Seurat_human_symbol_220723.rds')
mouse_multiome_Seurat_human_symbol <- mouse_multiome_Seurat_human_symbol[,cell_list]

p1 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p2 <- DimPlot(object = mouse_multiome_Seurat_human_symbol,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(ArchRPalettes$stallion))
p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = mouse_multiome_Seurat_human_symbol,file = './processed_data/220718_summary/mouse_multiome_Seurat_human_symbol_220723.rds')
