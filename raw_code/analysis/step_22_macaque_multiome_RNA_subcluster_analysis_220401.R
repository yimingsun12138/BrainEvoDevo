#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque multiome RNA subcluster analysis                        ##
## Data: 2022.03.15                                                                ##
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
library(parallel)
library(ArchR)
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(ComplexHeatmap)
library(dendextend)
library(Matrix)
library(matrixStats)
library(patchwork)
library(cowplot)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']
DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

meta_data <- readRDS(file = './res/step_22_fig_220330/macaque_multiome_Seurat_meta_data.rds')
macaque_multiome_Seurat$sub_cell_type <- meta_data[colnames(macaque_multiome_Seurat),"sub_cell_type"]
DimPlot(macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

macaque_multiome_Seurat <- NormalizeData(object = macaque_multiome_Seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)

# find marker -------------------------------------------------------------
#RG
RG_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('RG-1','RG-2','RG-3'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(RG_marker[RG_marker$pct.1 > 0.6 & RG_marker$pct.2 < 0.2,])
RG_marker <- RG_marker[rownames(RG_marker) %in% gene_list,]
write.csv(RG_marker,file = './res/step_22_fig_220401/RG_marker.csv')
pdf(file = './res/step_22_fig_220401/RG_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(RG_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#RG-1
RG_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(RG_1_marker[RG_1_marker$pct.1 > 0.6 & RG_1_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/RG_1_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-1',ident.2 = c('RG-2','RG-3'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
RG_1_marker <- RG_1_marker[rownames(RG_1_marker) %in% gene_list,]
RG_1_marker$subcluster_marker <- 'FALSE'
RG_1_marker[rownames(RG_1_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(RG_1_marker,file = './res/step_22_fig_220401/RG_1_marker.csv')
pdf(file = './res/step_22_fig_220401/RG_1_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(RG_1_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#RG-2
RG_2_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-2',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(RG_2_marker[RG_2_marker$pct.1 > 0.6 & RG_2_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/RG_2_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-2',ident.2 = c('RG-1','RG-3'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
RG_2_marker <- RG_2_marker[rownames(RG_2_marker) %in% gene_list,]
RG_2_marker$subcluster_marker <- 'FALSE'
RG_2_marker[rownames(RG_2_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(RG_2_marker,file = './res/step_22_fig_220401/RG_2_marker.csv')
pdf(file = './res/step_22_fig_220401/RG_2_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(RG_2_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#RG-3
RG_3_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-3',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(RG_3_marker[RG_3_marker$pct.1 > 0.6 & RG_3_marker$pct.2 < 0.2,])
temp <- 'EGFR'
gene_list <- dplyr::union(gene_list,temp)
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-3',ident.2 = c('RG-1','RG-2'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
RG_3_marker <- RG_3_marker[rownames(RG_3_marker) %in% gene_list,]
RG_3_marker$subcluster_marker <- 'FALSE'
RG_3_marker[rownames(RG_3_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(RG_3_marker,file = './res/step_22_fig_220401/RG_3_marker.csv')
pdf(file = './res/step_22_fig_220401/RG_3_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(RG_3_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex
Ex_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('Ex-1','Ex-2','Ex-3'),group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_marker[Ex_marker$pct.1 > 0.6 & Ex_marker$pct.2 < 0.2,])
Ex_marker <- Ex_marker[rownames(Ex_marker) %in% gene_list,]
write.csv(Ex_marker,file = './res/step_22_fig_220401/Ex_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-1_0
Ex_1_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_1_0_marker[Ex_1_0_marker$pct.1 > 0.6 & Ex_1_0_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_1_0_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1_0',ident.2 = c('Ex-1_1','Ex-2_0','Ex-2_1','Ex-3_0','Ex-3_1'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_1_0_marker <- Ex_1_0_marker[rownames(Ex_1_0_marker) %in% gene_list,]
Ex_1_0_marker$subcluster_marker <- 'FALSE'
Ex_1_0_marker[rownames(Ex_1_0_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_1_0_marker,file = './res/step_22_fig_220401/Ex_1_0_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_1_0_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_1_0_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-1_1
Ex_1_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_1_1_marker[Ex_1_1_marker$pct.1 > 0.6 & Ex_1_1_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_1_1_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1_1',ident.2 = c('Ex-1_0','Ex-2_0','Ex-2_1','Ex-3_0','Ex-3_1'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_1_1_marker <- Ex_1_1_marker[rownames(Ex_1_1_marker) %in% gene_list,]
Ex_1_1_marker$subcluster_marker <- 'FALSE'
Ex_1_1_marker[rownames(Ex_1_1_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_1_1_marker,file = './res/step_22_fig_220401/Ex_1_1_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_1_1_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_1_1_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-2_0
Ex_2_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_2_0_marker[Ex_2_0_marker$pct.1 > 0.6 & Ex_2_0_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_2_0_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2_0',ident.2 = c('Ex-1_1','Ex-1_0','Ex-2_1','Ex-3_0','Ex-3_1'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_2_0_marker <- Ex_2_0_marker[rownames(Ex_2_0_marker) %in% gene_list,]
Ex_2_0_marker$subcluster_marker <- 'FALSE'
Ex_2_0_marker[rownames(Ex_2_0_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_2_0_marker,file = './res/step_22_fig_220401/Ex_2_0_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_2_0_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_2_0_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-2_1
Ex_2_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_2_1_marker[Ex_2_1_marker$pct.1 > 0.6 & Ex_2_1_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_2_1_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2_1',ident.2 = c('Ex-1_1','Ex-2_0','Ex-1_0','Ex-3_0','Ex-3_1'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_2_1_marker <- Ex_2_1_marker[rownames(Ex_2_1_marker) %in% gene_list,]
Ex_2_1_marker$subcluster_marker <- 'FALSE'
Ex_2_1_marker[rownames(Ex_2_1_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_2_1_marker,file = './res/step_22_fig_220401/Ex_2_1_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_2_1_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_2_1_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-3_0
Ex_3_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_3_0_marker[Ex_3_0_marker$pct.1 > 0.6 & Ex_3_0_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_3_0_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3_0',ident.2 = c('Ex-1_1','Ex-2_0','Ex-2_1','Ex-1_0','Ex-3_1'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_3_0_marker <- Ex_3_0_marker[rownames(Ex_3_0_marker) %in% gene_list,]
Ex_3_0_marker$subcluster_marker <- 'FALSE'
Ex_3_0_marker[rownames(Ex_3_0_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_3_0_marker,file = './res/step_22_fig_220401/Ex_3_0_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_3_0_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_3_0_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-3_1
Ex_3_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_3_1_marker[Ex_3_1_marker$pct.1 > 0.6 & Ex_3_1_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_3_1_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3_1',ident.2 = c('Ex-1_1','Ex-2_0','Ex-2_1','Ex-3_0','Ex-1_0'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_3_1_marker <- Ex_3_1_marker[rownames(Ex_3_1_marker) %in% gene_list,]
Ex_3_1_marker$subcluster_marker <- 'FALSE'
Ex_3_1_marker[rownames(Ex_3_1_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_3_1_marker,file = './res/step_22_fig_220401/Ex_3_1_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_3_1_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_3_1_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#IP
IP_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('IP'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(IP_marker[IP_marker$pct.1 > 0.6 & IP_marker$pct.2 < 0.2,])
IP_marker <- IP_marker[rownames(IP_marker) %in% gene_list,]
write.csv(IP_marker,file = './res/step_22_fig_220401/IP_marker.csv')
pdf(file = './res/step_22_fig_220401/IP_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(IP_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-1
Ex_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_1_marker[Ex_1_marker$pct.1 > 0.6 & Ex_1_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_1_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1',ident.2 = c('Ex-2','Ex-3'),group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_1_marker <- Ex_1_marker[rownames(Ex_1_marker) %in% gene_list,]
Ex_1_marker$subcluster_marker <- 'FALSE'
Ex_1_marker[rownames(Ex_1_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_1_marker,file = './res/step_22_fig_220401/Ex_1_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_1_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_1_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-2
Ex_2_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_2_marker[Ex_2_marker$pct.1 > 0.6 & Ex_2_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_2_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2',ident.2 = c('Ex-1','Ex-3'),group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_2_marker <- Ex_2_marker[rownames(Ex_2_marker) %in% gene_list,]
Ex_2_marker$subcluster_marker <- 'FALSE'
Ex_2_marker[rownames(Ex_2_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_2_marker,file = './res/step_22_fig_220401/Ex_2_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_2_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_2_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#Ex-3
Ex_3_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(Ex_3_marker[Ex_3_marker$pct.1 > 0.6 & Ex_3_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/Ex_3_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3',ident.2 = c('Ex-1','Ex-2'),group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
Ex_3_marker <- Ex_3_marker[rownames(Ex_3_marker) %in% gene_list,]
Ex_3_marker$subcluster_marker <- 'FALSE'
Ex_3_marker[rownames(Ex_3_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(Ex_3_marker,file = './res/step_22_fig_220401/Ex_3_marker.csv')
pdf(file = './res/step_22_fig_220401/Ex_3_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_3_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#InCGE
InCGE_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(InCGE_marker[InCGE_marker$pct.1 > 0.6 & InCGE_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/InCGE_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE',ident.2 = c('InMGE'),group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
InCGE_marker <- InCGE_marker[rownames(InCGE_marker) %in% gene_list,]
InCGE_marker$subcluster_marker <- 'FALSE'
InCGE_marker[rownames(InCGE_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(InCGE_marker,file = './res/step_22_fig_220401/InCGE_marker.csv')
pdf(file = './res/step_22_fig_220401/InCGE_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InCGE_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#InMGE
InMGE_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(InMGE_marker[InMGE_marker$pct.1 > 0.6 & InMGE_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/InMGE_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE',ident.2 = c('InCGE'),group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
InMGE_marker <- InMGE_marker[rownames(InMGE_marker) %in% gene_list,]
InMGE_marker$subcluster_marker <- 'FALSE'
InMGE_marker[rownames(InMGE_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(InMGE_marker,file = './res/step_22_fig_220401/InMGE_marker.csv')
pdf(file = './res/step_22_fig_220401/InMGE_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#InCGE_0
InCGE_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(InCGE_0_marker[InCGE_0_marker$pct.1 > 0.6 & InCGE_0_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/InCGE_0_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE_0',ident.2 = c('InCGE_1','InMGE_0','InMGE_1','InMGE_2'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
InCGE_0_marker <- InCGE_0_marker[rownames(InCGE_0_marker) %in% gene_list,]
InCGE_0_marker$subcluster_marker <- 'FALSE'
InCGE_0_marker[rownames(InCGE_0_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(InCGE_0_marker,file = './res/step_22_fig_220401/InCGE_0_marker.csv')
pdf(file = './res/step_22_fig_220401/InCGE_0_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InCGE_0_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#InCGE_1
InCGE_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(InCGE_1_marker[InCGE_1_marker$pct.1 > 0.6 & InCGE_1_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/InCGE_1_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE_1',ident.2 = c('InCGE_0','InMGE_0','InMGE_1','InMGE_2'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
InCGE_1_marker <- InCGE_1_marker[rownames(InCGE_1_marker) %in% gene_list,]
InCGE_1_marker$subcluster_marker <- 'FALSE'
InCGE_1_marker[rownames(InCGE_1_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(InCGE_1_marker,file = './res/step_22_fig_220401/InCGE_1_marker.csv')
pdf(file = './res/step_22_fig_220401/InCGE_1_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InCGE_1_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#InMGE_0
InMGE_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(InMGE_0_marker[InMGE_0_marker$pct.1 > 0.6 & InMGE_0_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/InMGE_0_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_0',ident.2 = c('InCGE_1','InCGE_0','InMGE_1','InMGE_2'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
InMGE_0_marker <- InMGE_0_marker[rownames(InMGE_0_marker) %in% gene_list,]
InMGE_0_marker$subcluster_marker <- 'FALSE'
InMGE_0_marker[rownames(InMGE_0_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(InMGE_0_marker,file = './res/step_22_fig_220401/InMGE_0_marker.csv')
pdf(file = './res/step_22_fig_220401/InMGE_0_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_0_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#InMGE_1
InMGE_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(InMGE_1_marker[InMGE_1_marker$pct.1 > 0.6 & InMGE_1_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/InMGE_1_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_1',ident.2 = c('InCGE_0','InCGE_1','InMGE_0','InMGE_2'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
InMGE_1_marker <- InMGE_1_marker[rownames(InMGE_1_marker) %in% gene_list,]
InMGE_1_marker$subcluster_marker <- 'FALSE'
InMGE_1_marker[rownames(InMGE_1_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(InMGE_1_marker,file = './res/step_22_fig_220401/InMGE_1_marker.csv')
pdf(file = './res/step_22_fig_220401/InMGE_1_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_1_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()

#InMGE_2
InMGE_2_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_2',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
gene_list <- rownames(InMGE_2_marker[InMGE_2_marker$pct.1 > 0.6 & InMGE_2_marker$pct.2 < 0.2,])
temp <- read.csv(file = './res/step_22_fig_220330/InMGE_2_marker.csv')
gene_list <- dplyr::union(gene_list,temp[,1])
temp <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_2',ident.2 = c('InCGE_0','InCGE_1','InMGE_0','InMGE_1'),group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.2,]
gene_list <- dplyr::union(gene_list,rownames(temp))
InMGE_2_marker <- InMGE_2_marker[rownames(InMGE_2_marker) %in% gene_list,]
InMGE_2_marker$subcluster_marker <- 'FALSE'
InMGE_2_marker[rownames(InMGE_2_marker) %in% rownames(temp),"subcluster_marker"] <- 'TRUE'
write.csv(InMGE_2_marker,file = './res/step_22_fig_220401/InMGE_2_marker.csv')
pdf(file = './res/step_22_fig_220401/InMGE_2_marker_vlnplot.pdf',width = 24,height = round(length(gene_list)/3)*4)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_2_marker),assay = 'RNA',group.by = 'sub_cell_type',pt.size = 0,ncol = 3,slot = 'data')
dev.off()