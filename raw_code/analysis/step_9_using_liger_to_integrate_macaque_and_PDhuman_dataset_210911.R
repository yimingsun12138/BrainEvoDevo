#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: using liger to integrate macaque and PDhuman dataset            ##
## Data: 2021.09.11                                                                ##
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
library(rliger)
library(factoextra)
library(circlize)
library(ggsci)
library(SeuratWrappers)
library(magrittr)
library(harmony)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

# load data ---------------------------------------------------------------
# macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
# DefaultAssay(macaque_RNA_seurat) <- 'RNA'
# #modify the annotation
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('oRG','vRG'),"cell_type"] <- 'RG'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('30','9','1','4','22','2','12'),'cell_type'] <- 'Ex-1'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('14','8'),'cell_type'] <- 'Ex-2'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('16','15','24','6','18'),'cell_type'] <- 'Ex-3'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('19','20'),'cell_type'] <- 'Ex-4'
# 
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type == 'oRG',"sub_cell_type"] <- 'RG-1'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type == 'vRG-1',"sub_cell_type"] <- 'RG-2'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type == 'vRG-2',"sub_cell_type"] <- 'RG-3'
# 
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-1',"sub_cell_type"] <- 'Ex-1'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-2',"sub_cell_type"] <- 'Ex-2'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-3',"sub_cell_type"] <- 'Ex-3'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-4',"sub_cell_type"] <- 'Ex-4'
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
DefaultAssay(PDhuman_RNA_seurat) <- 'RNA'
PDhuman_RNA_seurat <- RenameCells(object = PDhuman_RNA_seurat,new.names = paste('human',colnames(PDhuman_RNA_seurat),sep = '_'))


# filter more low quality data --------------------------------------------
# DimPlot(macaque_RNA_seurat,label = TRUE,repel = TRUE,group.by = 'seurat_clusters')
# FeaturePlot(macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'))
# VlnPlot(macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'),group.by = 'seurat_clusters',pt.size = 0)
# macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('12','6'))]
# DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
# cell_list <- colnames(macaque_RNA_seurat)
# macaque_RNA_seurat$old_cell_type <- macaque_RNA_seurat$sub_cell_type
# macaque_RNA_seurat <- macaque_RNA_seurat[,cell_list]
# 
# #round 1
# meta_data <- macaque_RNA_seurat@meta.data[,c('percent.mt','species','donor','Phase','old_cell_type')]
# macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
# macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',
#                                          meta.data = meta_data,min.cells = 0,min.features = 0)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2500,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = c(0.5,1,1.8),group.by = 'old_cell_type',label = TRUE)
# DimPlot(object = macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
# FeaturePlot(macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'))
# VlnPlot(macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'),pt.size = 0,group.by = 'seurat_clusters')
# macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('20','25'))]
# 
# #round 2
# meta_data <- macaque_RNA_seurat@meta.data[,c('percent.mt','species','donor','Phase','old_cell_type')]
# macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
# macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',
#                                          meta.data = meta_data,min.cells = 0,min.features = 0)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 33,resolution = c(0.5,1,1.5),group.by = 'old_cell_type',label = TRUE)
# DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
# FeaturePlot(macaque_RNA_seurat,features = c('nCount_RNA','nFeature_RNA'))
# DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# macaque_RNA_seurat@meta.data[,"cell_type"] <- as.character(macaque_RNA_seurat$seurat_clusters)
# 
# #re annotation
# #round 1
# DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
# DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA',
#                              col.max = 2.5, col.min = -2.5, scale = TRUE,
#                              features = list(End=c('CLDN5','PECAM1'),
#                                              Per=c('PDGFRB'),
#                                              Mic=c('CX3CR1'),
#                                              Ast=c('AQP4','APOE'),
#                                              RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
#                                              OPC=c('SOX10','OLIG2','EGFR'),
#                                              Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
#                                              IP=c('EOMES','PPP1R17'),
#                                              Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
#                                              In=c('DLX5','GAD2','GAD1','DLX2'),
#                                              MGE=c('LHX6','SST'),
#                                              CGE=c('SP8','NR2F2'),
#                                              PSB=c('MEIS2','ETV1')),
#                              group.by = 'seurat_clusters', cols = c('#2CA02CFF','white','#D62728FF'),
#                              return_data_plot = TRUE)
# 
# my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
#            cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) +
#   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
#         panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
#         axis.line = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 14),
#         strip.background = element_rect(fill = 'grey',colour = 'black'),
#         legend.position = 'bottom',
#         panel.grid = element_line(color="grey",size = 0.1)) +
#   xlab('') + ylab('')
# 
# #cluster 27 Mic
# #cluster 23 31 End
# #cluster 29 Per
# #cluster 30 Astrocyte
# #cluster 20 OPC
# #cluster 12 Cyc
# #cluster 17 8 IP
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('27'),"cell_type"] <- 'Mic'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('23','31'),"cell_type"] <- 'End'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('29'),"cell_type"] <- 'Per'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('30'),"cell_type"] <- 'Astrocyte'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('20'),"cell_type"] <- 'OPC'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('12'),"cell_type"] <- 'Cyc'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('8','17'),"cell_type"] <- 'IP'
# 
# 
# #round2
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA',
#                              col.max = 2.5, col.min = -2.5, scale = TRUE,
#                              features = list(End=c('CLDN5','PECAM1'),
#                                              Per=c('PDGFRB'),
#                                              Mic=c('CX3CR1'),
#                                              Ast=c('AQP4','APOE'),
#                                              RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
#                                              OPC=c('SOX10','OLIG2','EGFR'),
#                                              Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
#                                              IP=c('EOMES','PPP1R17'),
#                                              Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
#                                              In=c('DLX5','GAD2','GAD1','DLX2'),
#                                              MGE=c('LHX6','SST'),
#                                              CGE=c('SP8','NR2F2'),
#                                              PSB=c('MEIS2','ETV1')),
#                              group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
#                              return_data_plot = TRUE)
# 
# my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
#            cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) +
#   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
#         panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
#         axis.line = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 14),
#         strip.background = element_rect(fill = 'grey',colour = 'black'),
#         legend.position = 'bottom',
#         panel.grid = element_line(color="grey",size = 0.1)) +
#   xlab('') + ylab('')
# 
# #cluster 3 28 25 2 10 InMGE
# #cluster 6 4 18 InCGE
# #cluster 26 InPSB
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('3','28','25','2','10'),"cell_type"] <- 'InMGE'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('6','4','18'),"cell_type"] <- 'InCGE'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('26'),"cell_type"] <- 'InPSB'
# 
# #round 3
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA',
#                              col.max = 2.5, col.min = -2.5, scale = TRUE,
#                              features = list(End=c('CLDN5','PECAM1'),
#                                              Per=c('PDGFRB'),
#                                              Mic=c('CX3CR1'),
#                                              Ast=c('AQP4','APOE'),
#                                              RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
#                                              OPC=c('SOX10','OLIG2','EGFR'),
#                                              Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
#                                              IP=c('EOMES','PPP1R17'),
#                                              Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
#                                              In=c('DLX5','GAD2','GAD1','DLX2'),
#                                              MGE=c('LHX6','SST'),
#                                              CGE=c('SP8','NR2F2'),
#                                              PSB=c('MEIS2','ETV1')),
#                              group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
#                              return_data_plot = TRUE)
# 
# my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
#            cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) +
#   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
#         panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
#         axis.line = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 14),
#         strip.background = element_rect(fill = 'grey',colour = 'black'),
#         legend.position = 'bottom',
#         panel.grid = element_line(color="grey",size = 0.1)) +
#   xlab('') + ylab('')
# 
# #cluster 0 1 5 Ex-1
# #cluster 7 Ex-2
# #cluster 9 11 14 Ex-3
# #cluster 21 24 15 Ex-4
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('0','1','5'),"cell_type"] <- 'Ex-1'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('7'),"cell_type"] <- 'Ex-2'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('9','11','14'),"cell_type"] <- 'Ex-3'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('21','24','15'),"cell_type"] <- 'Ex-4'
# 
# #round 4
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA',
#                              col.max = 2.5, col.min = -2.5, scale = TRUE,
#                              features = list(End=c('CLDN5','PECAM1'),
#                                              Per=c('PDGFRB'),
#                                              Mic=c('CX3CR1'),
#                                              Ast=c('AQP4','APOE'),
#                                              RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
#                                              OPC=c('SOX10','OLIG2','EGFR'),
#                                              Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
#                                              IP=c('EOMES','PPP1R17'),
#                                              Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
#                                              In=c('DLX5','GAD2','GAD1','DLX2'),
#                                              MGE=c('LHX6','SST'),
#                                              CGE=c('SP8','NR2F2'),
#                                              PSB=c('MEIS2','ETV1')),
#                              group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
#                              return_data_plot = TRUE)
# 
# my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
#            cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) +
#   theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
#         panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
#         axis.line = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text = element_text(size = 12),
#         strip.text = element_text(size = 14),
#         strip.background = element_rect(fill = 'grey',colour = 'black'),
#         legend.position = 'bottom',
#         panel.grid = element_line(color="grey",size = 0.1)) +
#   xlab('') + ylab('')
# 
# #cluster 16 22 13 RG
# #cluster 19 Ex-3
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('16','22','13'),"cell_type"] <- 'RG'
# macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters %in% c('19'),"cell_type"] <- 'Ex-3'
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# 
# 
# #do sub cluster
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# Idents(macaque_RNA_seurat) <- 'cell_type'
# macaque_RNA_seurat@meta.data[,'sub_cell_type'] <- as.character(macaque_RNA_seurat$cell_type)
# #Cyc
# temp <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type == 'Cyc']
# temp <- RunPCA(object = temp,assay = 'RNA',npcs = 50)
# ElbowPlot(temp,ndims = 50)
# temp <- FindNeighbors(object = temp,dims = 1:5)
# temp <- FindClusters(object = temp,resolution = 0.3)
# p1 <- DimPlot(temp,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
# p2 <- DimPlot(temp,group.by = 'Phase',label = TRUE,repel = TRUE)
# p1+p2
# 
# macaque_RNA_seurat@meta.data[colnames(temp[,temp$seurat_clusters %in% c('2','3')]),"sub_cell_type"] <- 'Cyc-G2M'
# macaque_RNA_seurat@meta.data[colnames(temp[,temp$seurat_clusters %in% c('0','1','4')]),"sub_cell_type"] <- 'Cyc-S'
# DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
# #RG
# temp <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type == 'RG']
# temp <- RunPCA(object = temp,assay = 'RNA',npcs = 50)
# ElbowPlot(temp,ndims = 50)
# temp <- FindNeighbors(object = temp,dims = 1:4)
# temp <- FindClusters(object = temp,resolution = 0.6)
# p1 <- DimPlot(temp,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
# p2 <- DimPlot(temp,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# p1+p2
# 
# macaque_RNA_seurat@meta.data[colnames(temp[,temp$seurat_clusters %in% c('5','7')]),"sub_cell_type"] <- 'RG-3'
# macaque_RNA_seurat@meta.data[colnames(temp[,temp$seurat_clusters %in% c('0','6')]),"sub_cell_type"] <- 'RG-2'
# macaque_RNA_seurat@meta.data[colnames(temp[,temp$seurat_clusters %in% c('1','2','3','4','8')]),"sub_cell_type"] <- 'RG-1'
# DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
# 
# #plot
# color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
# color_paramater <- color_paramater[c(1:22,25,26),]
# color_paramater[which(color_paramater$id == 'oRG'),'id'] <- 'RG'
# color_paramater[which(color_paramater$id == 'vRG'),'id'] <- 'Ex-2'
# color_paramater[which(color_paramater$id == 'Ex'),'id'] <- 'Ex-1'
# color_paramater[which(color_paramater$id == 'Ex-U'),'id'] <- 'Ex-3'
# color_paramater[which(color_paramater$id == 'Ex-SP'),'id'] <- 'Ex-4'
# 
# temp_col <- data.frame(id = c('RG-1','RG-2','RG-3'),col = colorRampPalette(c("#00A087FF","#3C5488FF"))(5)[2:4])
# color_paramater <- rbind(color_paramater,temp_col)
# write.csv(color_paramater,file = './data/parameter/color_paramater.csv')
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id
# DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cols = temp_col)
# DefaultAssay(macaque_RNA_seurat) <- 'RNA'
# saveRDS(macaque_RNA_seurat,file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210913.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$sub_cell_type %in% c('InPSB','Astrocyte'))]
# convert gene name -------------------------------------------------------
dim(macaque_RNA_seurat)
dim(PDhuman_RNA_seurat)

macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')
table(rownames(PDhuman_RNA_seurat) %in% c(human_to_human_anno[,1],human_to_human_anno[,2]))

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
dim(temp)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp

temp <- PDhuman_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno,filter_anno = FALSE)
dim(temp)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
PDhuman_RNA_seurat[['converted']] <- temp

# integrate using liger ---------------------------------------------------
# #create liger object
# Brain_RNA <- createLiger(list(human=PDhuman_RNA_seurat@assays$converted@counts,macaque=macaque_RNA_seurat@assays$converted@counts))
# Brain_RNA <- rliger::normalize(Brain_RNA)
# Brain_RNA <- rliger::selectGenes(Brain_RNA)
# Brain_RNA <- rliger::scaleNotCenter(Brain_RNA)
# 
# Brain_RNA <- optimizeALS(Brain_RNA,lambda = 20,k = 30)
# Brain_RNA <- quantile_norm(Brain_RNA,knn_k = 20,quantiles = 50,min_cells = 20,
#                            do.center = FALSE,max_sample = 1000,refine.knn = TRUE,eps = 0.9)
# Brain_RNA <- louvainCluster(Brain_RNA,resolution = 0.4)
# Brain_RNA <- runUMAP(Brain_RNA,distance = 'cosine',min_dist = 0.1)
# all.plots <- plotByDatasetAndCluster(Brain_RNA, axis.labels = c('UMAP 1','UMAP 2'), return.plots = T)

# try harmony -------------------------------------------------------------
macaque_RNA_seurat$dataset <- 'macaque'
PDhuman_RNA_seurat$dataset <- 'human'
#find variable gene
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',nfeatures = 6000,npcs = 50,preprocess = TRUE)
PDhuman_RNA_seurat <- my_process_seurat(object = PDhuman_RNA_seurat,assay = 'converted',nfeatures = 6000,npcs = 50,preprocess = TRUE)
gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),VariableFeatures(PDhuman_RNA_seurat))

#integration
Brain_RNA_seurat <- my_harmony_integration(named_seurat_list = list(macaque=macaque_RNA_seurat,human=PDhuman_RNA_seurat),
                                           assay = 'converted',variable_feature = gene_list,
                                           var_to_regress_list = list(macaque=c('donor','nCount_converted'),human=c('Donor','Library','nCount_converted')),
                                           npcs = 50,reference_loading = 'macaque',integration_var = 'dataset',harmony_input_dim = 40,max.iter.harmony = 50,
                                           reference_dataset = NULL,UMAP_dim = 40,resolution = 1,kmeans_init_nstart = 10,kmeans_init_iter_max = 200,sigma = 0.3)

Brain_RNA_seurat@meta.data[,'cell_type'] <- NA
Brain_RNA_seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$sub_cell_type
Brain_RNA_seurat@meta.data[colnames(PDhuman_RNA_seurat),"cell_type"] <- PDhuman_RNA_seurat$Cluster

Brain_RNA_seurat <- my_process_seurat(object = Brain_RNA_seurat,assay = 'integration',preprocess = FALSE,dim_to_use = 40,resolution = 1,group.by = 'dataset',label = FALSE)
DimPlot(Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset')
DimPlot(Brain_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,split.by = 'dataset')

pdf(file = './res/step_9_fig_210911/macaque_PDhuman_harmony_integration.pdf',width = 12,height = 6)
DimPlot(Brain_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        axis.line = element_blank(),
        legend.position = 'none') + 
  labs(title = '')
dev.off()
#label transfer
predicted_table <- my_knn_label_transfer(data = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'macaque'],
                                         query = Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human'],
                                         reference_var = 'cell_type',query_var = 'predicted_label',
                                         reduction = 'pca',knn_k = 500,knn_n = 500,n_core = 4,replace = TRUE,SNN = TRUE)
saveRDS(predicted_table,file = './res/step_9_fig_210911/macaque_PDhuman_SNN_k_500_n_500_predicted_table.rds')
histogram(predicted_table$predicted_score)
table(predicted_table$predicted_score > 0.5)

#recreate PDhuman_RNA_seurat
temp <- PDhuman_RNA_seurat@assays$RNA@counts
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
meta_data <- PDhuman_RNA_seurat@meta.data
meta_data <- meta_data[,c('Subcluster','Donor','Layer','Gestation_week','Index','Library','Percentage_mitochondrial','Phase')]
PDhuman_RNA_seurat <- Brain_RNA_seurat[,Brain_RNA_seurat$dataset == 'human']

meta_data <- meta_data[colnames(PDhuman_RNA_seurat),]
predicted_table <- predicted_table[colnames(PDhuman_RNA_seurat),]

PDhuman_RNA_seurat <- AddMetaData(object = PDhuman_RNA_seurat,metadata = meta_data)
PDhuman_RNA_seurat <- AddMetaData(PDhuman_RNA_seurat,metadata = predicted_table[,c(1,2)])
PDhuman_RNA_seurat[['original']] <- temp
PDhuman_RNA_seurat <- my_process_seurat(object = PDhuman_RNA_seurat,assay = 'original',reduction.name = 'PCA',nfeatures = 3000,vars.to.regress = c('nCount_original','Donor','Library'),npcs = 50,preprocess = TRUE)
DimPlot(PDhuman_RNA_seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE,reduction = 'umap',cols = temp_col)
FeaturePlot(PDhuman_RNA_seurat,features = c('EGFR','predicted_score'),reduction = 'umap',slot = 'data')
