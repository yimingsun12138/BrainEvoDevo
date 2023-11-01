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
DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']

macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

macaque_multiome_Seurat$sub_cell_type <- macaque_multiome_Seurat$cell_type
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = 1,group.by = 'cell_type',label = TRUE)
# try subcluster at Ex-3 --------------------------------------------------
#do subcluster
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'Ex-3',resolution = 0.1,graph.name = 'RNA_snn')
DimPlot(macaque_multiome_Seurat,group.by = 'sub.cluster',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

macaque_multiome_ArchR$sub.cluster <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub.cluster"]
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "sub.cluster", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'Ex-3',"sub_cell_type"] <- macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'Ex-3',"sub.cluster"]

#find marker
Ex_3_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_3_0_marker <- Ex_3_0_marker[Ex_3_0_marker$pct.1 > 0.6 & Ex_3_0_marker$pct.2 < 0.4,]
temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'Ex-3'],ident.1 = 'Ex-3_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
Ex_3_0_marker <- Ex_3_0_marker[rownames(Ex_3_0_marker) %in% rownames(temp),]
write.csv(Ex_3_0_marker,file = './res/step_22_fig_220330/Ex_3_0_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_3_0_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_3_0_marker)[1:7],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

Ex_3_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_3_1_marker <- Ex_3_1_marker[Ex_3_1_marker$pct.1 > 0.6 & Ex_3_1_marker$pct.2 < 0.4,]
temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'Ex-3'],ident.1 = 'Ex-3_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
Ex_3_1_marker <- Ex_3_1_marker[rownames(Ex_3_1_marker) %in% rownames(temp),]
write.csv(Ex_3_1_marker,file = './res/step_22_fig_220330/Ex_3_1_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_3_1_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_3_1_marker)[1:12],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

Ex_3_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-3',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_3_marker <- Ex_3_marker[Ex_3_marker$pct.1 > 0.6 & Ex_3_marker$pct.2 < 0.2,]
write.csv(Ex_3_marker,file = './res/step_22_fig_220330/Ex_3_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_3_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_3_marker)[1:12],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

pdf(file = './res/step_22_fig_220330/Ex_3_marker_dotplot.pdf',width = 6,height = 6)
my_dotplot(object = macaque_multiome_Seurat,assay = 'RNA',
           features = list(Ex_3=c('GRIK3','FAM81A','KIAA1217','TRPM3'),
                           Ex_3_0=c('EPHB1','GLRA2','EPHA6'),
                           Ex_3_1=c('CLSTN2','POU6F2','CDH18','CDH9')),
           cols = c('#3B4992FF','white','#EE0000FF'),col.max = 2.5,col.min = -2.5,scale = TRUE,
           group.by = 'sub.cluster',return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

# try subcluster for Ex-2 --------------------------------------------------
#do subcluster
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'Ex-2',resolution = 0.05,graph.name = 'RNA_snn')
DimPlot(macaque_multiome_Seurat,group.by = 'sub.cluster',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

macaque_multiome_ArchR$sub.cluster <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub.cluster"]
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "sub.cluster", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'Ex-2',"sub_cell_type"] <- macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'Ex-2',"sub.cluster"]

#find marker
Ex_2_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_2_0_marker <- Ex_2_0_marker[Ex_2_0_marker$pct.1 > 0.6 & Ex_2_0_marker$pct.2 < 0.4,]
temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'Ex-2'],ident.1 = 'Ex-2_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
Ex_2_0_marker <- Ex_2_0_marker[rownames(Ex_2_0_marker) %in% rownames(temp),]
write.csv(Ex_2_0_marker,file = './res/step_22_fig_220330/Ex_2_0_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_2_0_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_2_0_marker)[1:6],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

Ex_2_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_2_1_marker <- Ex_2_1_marker[Ex_2_1_marker$pct.1 > 0.6 & Ex_2_1_marker$pct.2 < 0.4,]
temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'Ex-2'],ident.1 = 'Ex-2_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
Ex_2_1_marker <- Ex_2_1_marker[rownames(Ex_2_1_marker) %in% rownames(temp),]
write.csv(Ex_2_1_marker,file = './res/step_22_fig_220330/Ex_2_1_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_2_1_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_2_1_marker)[1:12],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

Ex_2_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-2',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_2_marker <- Ex_2_marker[Ex_2_marker$pct.1 > 0.6 & Ex_2_marker$pct.2 < 0.2,]
write.csv(Ex_2_marker,file = './res/step_22_fig_220330/Ex_2_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_2_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_2_marker)[13:22],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

pdf(file = './res/step_22_fig_220330/Ex_2_marker_dotplot.pdf',width = 6,height = 6)
my_dotplot(object = macaque_multiome_Seurat,assay = 'RNA',
           features = list(Ex_2=c('CHRM3','NWD2','SNTG1'),
                           Ex_2_0=c('VWC2L','SDK1'),
                           Ex_2_1=c('ARAP2','GALNT17','POU6F2','CDH12')),
           cols = c('#3B4992FF','white','#EE0000FF'),col.max = 2.5,col.min = -2.5,scale = TRUE,
           group.by = 'sub.cluster',return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

# try subcluster for Ex-1 -------------------------------------------------
#do subcluster
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'Ex-1',resolution = 0.1,graph.name = 'RNA_snn')
DimPlot(macaque_multiome_Seurat,group.by = 'sub.cluster',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

macaque_multiome_ArchR$sub.cluster <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub.cluster"]
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "sub.cluster", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'Ex-1',"sub_cell_type"] <- macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'Ex-1',"sub.cluster"]

#find marker
Ex_1_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_1_0_marker <- Ex_1_0_marker[Ex_1_0_marker$pct.1 > 0.6 & Ex_1_0_marker$pct.2 < 0.4,]
temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'Ex-1'],ident.1 = 'Ex-1_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
Ex_1_0_marker <- Ex_1_0_marker[rownames(Ex_1_0_marker) %in% rownames(temp),]
write.csv(Ex_1_0_marker,file = './res/step_22_fig_220330/Ex_1_0_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_1_0_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_1_0_marker)[1:12],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

Ex_1_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_1_1_marker <- Ex_1_1_marker[Ex_1_1_marker$pct.1 > 0.6 & Ex_1_1_marker$pct.2 < 0.4,]
# temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'Ex-1'],ident.1 = 'Ex-1_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
# temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
# Ex_1_1_marker <- Ex_1_1_marker[rownames(Ex_1_1_marker) %in% rownames(temp),]
write.csv(Ex_1_1_marker,file = './res/step_22_fig_220330/Ex_1_1_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_1_1_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_1_1_marker)[1:12],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

Ex_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'Ex-1',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
Ex_1_marker <- Ex_1_marker[Ex_1_marker$pct.1 > 0.6 & Ex_1_marker$pct.2 < 0.4,]
write.csv(Ex_1_marker,file = './res/step_22_fig_220330/Ex_1_marker.csv')
pdf(file = './res/step_22_fig_220330/Ex_1_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(Ex_1_marker)[1:14],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

pdf(file = './res/step_22_fig_220330/Ex_1_marker_dotplot.pdf',width = 6,height = 6)
my_dotplot(object = macaque_multiome_Seurat,assay = 'RNA',
           features = list(Ex_1=c('CUX2','VCAN','SLC22A23','HS6ST3'),
                           Ex_1_0=c('KIF26B','ITPR1','CCBE1'),
                           Ex_1_1=c('THSD7A','SEMA3C','PGAP1')),
           cols = c('#3B4992FF','white','#EE0000FF'),col.max = 2.5,col.min = -2.5,scale = TRUE,
           group.by = 'sub.cluster',return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

# try subcluster at InCGE --------------------------------------------------
#do subcluster
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'InCGE',resolution = 0.1,graph.name = 'RNA_snn')
DimPlot(macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'InCGE'],group.by = 'sub.cluster',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
# macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'InCGE_2',"sub.cluster"] <- 'InCGE_0'

macaque_multiome_ArchR$sub.cluster <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub.cluster"]
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "sub.cluster", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'InCGE',"sub_cell_type"] <- macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'InCGE',"sub.cluster"]

#find marker
InCGE_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
InCGE_0_marker <- InCGE_0_marker[InCGE_0_marker$pct.1 > 0.6 & InCGE_0_marker$pct.2 < 0.4,]
# temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'InCGE'],ident.1 = 'InCGE_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
# temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
# InCGE_0_marker <- InCGE_0_marker[rownames(InCGE_0_marker) %in% rownames(temp),]
write.csv(InCGE_0_marker,file = './res/step_22_fig_220330/InCGE_0_marker.csv')
pdf(file = './res/step_22_fig_220330/InCGE_0_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InCGE_0_marker)[1:8],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

InCGE_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
InCGE_1_marker <- InCGE_1_marker[InCGE_1_marker$pct.1 > 0.6 & InCGE_1_marker$pct.2 < 0.4,]
temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'InCGE'],ident.1 = 'InCGE_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
InCGE_1_marker <- InCGE_1_marker[rownames(InCGE_1_marker) %in% rownames(temp),]
write.csv(InCGE_1_marker,file = './res/step_22_fig_220330/InCGE_1_marker.csv')
pdf(file = './res/step_22_fig_220330/InCGE_1_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InCGE_1_marker)[1:8],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

InCGE_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InCGE',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
InCGE_marker <- InCGE_marker[InCGE_marker$pct.1 > 0.6 & InCGE_marker$pct.2 < 0.3,]
write.csv(InCGE_marker,file = './res/step_22_fig_220330/InCGE_marker.csv')
pdf(file = './res/step_22_fig_220330/InCGE_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InCGE_marker),assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

pdf(file = './res/step_22_fig_220330/InCGE_marker_dotplot.pdf',width = 6,height = 6)
my_dotplot(object = macaque_multiome_Seurat,assay = 'RNA',
           features = list(InCGE=c('THRB','ADARB2'),
                           InCGE_0=c('PDZRN3','CHD7'),
                           InCGE_1=c('SYNPR','NR3C2','GALNTL6','CNR1')),
           cols = c('#3B4992FF','white','#EE0000FF'),col.max = 2.5,col.min = -2.5,scale = TRUE,
           group.by = 'sub.cluster',return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

# try subcluster at InMGE --------------------------------------------------
#do subcluster
Idents(macaque_multiome_Seurat) <- 'cell_type'
macaque_multiome_Seurat <- FindSubCluster(object = macaque_multiome_Seurat,cluster = 'InMGE',resolution = 0.2,graph.name = 'RNA_snn')
DimPlot(macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'InMGE'],group.by = 'sub.cluster',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
# macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$sub.cluster == 'InMGE_3',"sub.cluster"] <- 'InMGE_2'

macaque_multiome_ArchR$sub.cluster <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub.cluster"]
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "sub.cluster", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'InMGE',"sub_cell_type"] <- macaque_multiome_Seurat@meta.data[macaque_multiome_Seurat$cell_type == 'InMGE',"sub.cluster"]

#find marker
InMGE_0_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
InMGE_0_marker <- InMGE_0_marker[InMGE_0_marker$pct.1 > 0.6 & InMGE_0_marker$pct.2 < 0.4,]
# temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'InMGE'],ident.1 = 'InMGE_0',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
# temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
# InMGE_0_marker <- InMGE_0_marker[rownames(InMGE_0_marker) %in% rownames(temp),]
write.csv(InMGE_0_marker,file = './res/step_22_fig_220330/InMGE_0_marker.csv')
pdf(file = './res/step_22_fig_220330/InMGE_0_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_0_marker),assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

InMGE_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
InMGE_1_marker <- InMGE_1_marker[InMGE_1_marker$pct.1 > 0.6 & InMGE_1_marker$pct.2 < 0.4,]
# temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'InMGE'],ident.1 = 'InMGE_1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
# temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
# InMGE_1_marker <- InMGE_1_marker[rownames(InMGE_1_marker) %in% rownames(temp),]
write.csv(InMGE_1_marker,file = './res/step_22_fig_220330/InMGE_1_marker.csv')
pdf(file = './res/step_22_fig_220330/InMGE_1_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_1_marker),assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

InMGE_2_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE_2',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
InMGE_2_marker <- InMGE_2_marker[InMGE_2_marker$pct.1 > 0.6 & InMGE_2_marker$pct.2 < 0.4,]
temp <- FindMarkers(object = macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type == 'InMGE'],ident.1 = 'InMGE_2',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
temp <- temp[temp$pct.1 > 0.6 & temp$pct.2 < 0.4,]
InMGE_2_marker <- InMGE_2_marker[rownames(InMGE_2_marker) %in% rownames(temp),]
write.csv(InMGE_2_marker,file = './res/step_22_fig_220330/InMGE_2_marker.csv')
pdf(file = './res/step_22_fig_220330/InMGE_2_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_2_marker)[13:24],assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

InMGE_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'InMGE',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
InMGE_marker <- InMGE_marker[InMGE_marker$pct.1 > 0.6 & InMGE_marker$pct.2 < 0.2,]
write.csv(InMGE_marker,file = './res/step_22_fig_220330/InMGE_marker.csv')
pdf(file = './res/step_22_fig_220330/InMGE_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(InMGE_marker),assay = 'RNA',group.by = 'sub_cell_type',slot = 'data',pt.size = 0)
dev.off()

pdf(file = './res/step_22_fig_220330/InMGE_marker_dotplot.pdf',width = 6,height = 6)
my_dotplot(object = macaque_multiome_Seurat,assay = 'RNA',
           features = list(InMGE=c('NXPH1'),
                           InMGE_0=c('DACH2','EML6'),
                           InMGE_1=c('RBMS3','PAM'),
                           InMGE_2=c('SATB1','PTCHD4','ELAVL2','NETO1')),
           cols = c('#3B4992FF','white','#EE0000FF'),col.max = 2.5,col.min = -2.5,scale = TRUE,
           group.by = 'sub.cluster',return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

#save data
saveRDS(object = macaque_multiome_Seurat@meta.data,file = './res/step_22_fig_220330/macaque_multiome_Seurat_meta_data.rds')

# final subcluster --------------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']
meta_data <- readRDS(file = './res/step_22_fig_220330/macaque_multiome_Seurat_meta_data.rds')
macaque_multiome_Seurat$sub_cell_type <- meta_data[colnames(macaque_multiome_Seurat),"sub_cell_type"]
DimPlot(macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_multiome_ArchR$sub_cell_type <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub_cell_type"]

#dimplot
RNA_UMAP <- macaque_multiome_Seurat@reductions$harmonyUMAP@cell.embeddings
ATAC_UMAP <- macaque_multiome_ArchR@embeddings$UMAP@listData$df

p1 <- my_dimplot(embedding = RNA_UMAP,meta_data = macaque_multiome_Seurat@meta.data,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = c(as.character(ArchRPalettes$stallion),'#7DD06F','#844081')) + 
  theme(aspect.ratio = 1)
p2 <- my_dimplot(embedding = ATAC_UMAP,meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + 
  scale_color_manual(values = c(as.character(ArchRPalettes$stallion),'#7DD06F','#844081')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_ArchR_dimplot.pdf',width = 16,height = 7)
p1+p2+plot_layout(ncol = 2)
dev.off()

# find marker -------------------------------------------------------------
macaque_multiome_Seurat <- NormalizeData(object = macaque_multiome_Seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)

#RG-1
RG_1_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
RG_1_marker <- RG_1_marker[RG_1_marker$pct.1 > 0.6 & RG_1_marker$pct.2 < 0.1,]
write.csv(RG_1_marker,file = './res/step_22_fig_220330/RG_1_marker.csv')

pdf(file = './res/step_22_fig_220330/RG_1_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(RG_1_marker)[1:12],pt.size = 0,group.by = 'sub_cell_type')
dev.off()

#RG-2
RG_2_marker <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = 'RG-2',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
RG_2_marker <- RG_2_marker[RG_2_marker$pct.1 > 0.6 & RG_2_marker$pct.2 < 0.1,]
write.csv(RG_2_marker,file = './res/step_22_fig_220330/RG_2_marker.csv')

pdf(file = './res/step_22_fig_220330/RG_2_marker_vlnplot.pdf',width = 14,height = 10)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(RG_2_marker)[13:14],pt.size = 0,group.by = 'sub_cell_type')
dev.off()

#dotplot
macaque_multiome_Seurat <- ScaleData(object = macaque_multiome_Seurat,features = rownames(macaque_multiome_Seurat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE)
dotplot_matrix <- my_dotplot(object = macaque_multiome_Seurat,assay = 'RNA',
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('FOXJ1'),
                                             RG=c('SOX9','NOTCH2','FGFR2','AQP4','EDNRB','EGFR'),
                                             OPC=c('SOX10','OLIG2'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6'),
                                             Ex_1=c('CUX2','HS6ST3','KIF26B','THSD7A'),
                                             Ex_2=c('NWD2','SNTG1','VWC2L','CDH12'),
                                             Ex_3=c('GRIK3','TRPM3','GLRA2','CLSTN2'),
                                             In=c('GAD1','GAD2'),
                                             InMGE=c('NXPH1','DACH2','RBMS3','PTCHD4'),
                                             InCGE=c('THRB','PDZRN3','NR3C2','CNR1')),
                             cols = c('#3B4992FF','white','#EE0000FF'),col.max = 2.5,col.min = -2.5,scale = TRUE,
                             group.by = 'sub_cell_type',return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE_1','InCGE_0','InMGE_0','InMGE_1','InMGE_2','Ex-3_1','Ex-3_0','Ex-2_1','Ex-2_0','Ex-1_0','Ex-1_1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))

pdf(file = './res/step_22_fig_220330/macaque_multiome_de_novel_marker_dotplot.pdf',width = 16,height = 7)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

# layer analysis ----------------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$tech == 'multiome']
meta_data <- readRDS(file = './res/step_22_fig_220330/macaque_multiome_Seurat_meta_data.rds')
macaque_multiome_Seurat$sub_cell_type <- meta_data[colnames(macaque_multiome_Seurat),"sub_cell_type"]
DimPlot(macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_multiome_ArchR$sub_cell_type <- macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"sub_cell_type"]

macaque_multiome_Seurat <- NormalizeData(object = macaque_multiome_Seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)

#classic marker
SP_marker <- c('CDH18','HAS3','HS3ST4','NXPH3','SLC35F3','ZFPM2')
SP_marker %in% rownames(macaque_multiome_Seurat)

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_SP_classic_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = SP_marker,object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.title = 'SP classic marker',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

CPo_marker <- c('DACT1','GALNT3','MAN1C1','MYCN','RHOU','THG1L','TMOD1')
CPo_marker %in% rownames(macaque_multiome_Seurat)

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_CPo_classic_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = CPo_marker,object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.title = 'CPo classic marker',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

CPi_marker <- c('LMO4','CD1D','CYTIP','FOXP1','HS3ST2','PCDH19','PTPRK','SAMD5','TRAF5')
CPi_marker %in% rownames(macaque_multiome_Seurat)

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_CPi_classic_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = CPi_marker,object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.title = 'CPi classic marker',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

IZ_marker <- c('ADM','COL20A1','NEU4','PDGFRA','PPP1R3G','SP5')
IZ_marker %in% rownames(macaque_multiome_Seurat)

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_IZ_classic_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = IZ_marker,object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.title = 'IZ classic marker',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

OSVZ_marker <- c('ACSBG1','BTBD17','LYN','PLAC9','PLAGL1','PREX2','TNC','EOMES')
OSVZ_marker %in% rownames(macaque_multiome_Seurat)

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_OSVZ_classic_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = OSVZ_marker,object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.title = 'OSVZ classic marker',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

ISVZ_marker <- c('PAX6','CDK6','GAS1','JAG1','PTH2R','SETD7')
ISVZ_marker %in% rownames(macaque_multiome_Seurat)

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_ISVZ_classic_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = ISVZ_marker,object = macaque_multiome_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.title = 'ISVZ classic marker',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

#bulk marker
marker_list <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')

#layer marker express
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

#CPo
CPo_marker <- marker_list[['CPo']]
CPo_marker <- unique(CPo_marker)
table(CPo_marker %in% macaque_annotation$gene_id)
CPo_marker <- CPo_marker[CPo_marker %in% macaque_annotation$gene_id]
CPo_marker <- macaque_annotation$gene_name[base::match(CPo_marker,table = macaque_annotation$gene_id)]
table(CPo_marker %in% rownames(macaque_multiome_Seurat))
CPo_marker <- CPo_marker[CPo_marker %in% rownames(macaque_multiome_Seurat)]
table(duplicated(CPo_marker))
CPo_marker <- unique(CPo_marker)
CPo_plot <- geneSetAveragePlot(genes = CPo_marker,object = macaque_multiome_Seurat,
                              object.class = 'seurat',plot.type = 'average',embedding = 'harmonyUMAP',assay = 'RNA',
                              aspectratio = 1,trim = c(0,0.99),color.palette = c("lightgrey",'red'),
                              print = FALSE,scaled = FALSE,plot.title = 'CPo bulk marker average plot')

#CPi
CPi_marker <- marker_list[['CPi']]
CPi_marker <- unique(CPi_marker)
table(CPi_marker %in% macaque_annotation$gene_id)
CPi_marker <- CPi_marker[CPi_marker %in% macaque_annotation$gene_id]
CPi_marker <- macaque_annotation$gene_name[base::match(CPi_marker,table = macaque_annotation$gene_id)]
table(CPi_marker %in% rownames(macaque_multiome_Seurat))
CPi_marker <- CPi_marker[CPi_marker %in% rownames(macaque_multiome_Seurat)]
table(duplicated(CPi_marker))
CPi_marker <- unique(CPi_marker)
CPi_plot <- geneSetAveragePlot(genes = CPi_marker,object = macaque_multiome_Seurat,
                              object.class = 'seurat',plot.type = 'average',embedding = 'harmonyUMAP',assay = 'RNA',
                              aspectratio = 1,trim = c(0,0.99),color.palette = c("lightgrey",'red'),
                              print = FALSE,scaled = FALSE,plot.title = 'CPi bulk marker average plot')

#SP
SP_marker <- marker_list[['SP']]
SP_marker <- unique(SP_marker)
table(SP_marker %in% macaque_annotation$gene_id)
SP_marker <- SP_marker[SP_marker %in% macaque_annotation$gene_id]
SP_marker <- macaque_annotation$gene_name[base::match(SP_marker,table = macaque_annotation$gene_id)]
table(SP_marker %in% rownames(macaque_multiome_Seurat))
SP_marker <- SP_marker[SP_marker %in% rownames(macaque_multiome_Seurat)]
table(duplicated(SP_marker))
SP_marker <- unique(SP_marker)
SP_plot <- geneSetAveragePlot(genes = SP_marker,object = macaque_multiome_Seurat,
                              object.class = 'seurat',plot.type = 'average',embedding = 'harmonyUMAP',assay = 'RNA',
                              aspectratio = 1,trim = c(0,0.99),color.palette = c("lightgrey",'red'),
                              print = FALSE,scaled = FALSE,plot.title = 'SP bulk marker average plot')

#IZ
IZ_marker <- marker_list[['IZ']]
IZ_marker <- unique(IZ_marker)
table(IZ_marker %in% macaque_annotation$gene_id)
IZ_marker <- IZ_marker[IZ_marker %in% macaque_annotation$gene_id]
IZ_marker <- macaque_annotation$gene_name[base::match(IZ_marker,table = macaque_annotation$gene_id)]
table(IZ_marker %in% rownames(macaque_multiome_Seurat))
IZ_marker <- IZ_marker[IZ_marker %in% rownames(macaque_multiome_Seurat)]
table(duplicated(IZ_marker))
IZ_marker <- unique(IZ_marker)
IZ_plot <- geneSetAveragePlot(genes = IZ_marker,object = macaque_multiome_Seurat,
                              object.class = 'seurat',plot.type = 'average',embedding = 'harmonyUMAP',assay = 'RNA',
                              aspectratio = 1,trim = c(0,0.99),color.palette = c("lightgrey",'red'),
                              print = FALSE,scaled = FALSE,plot.title = 'IZ bulk marker average plot')

#OSVZ
OSVZ_marker <- marker_list[['OSVZ']]
OSVZ_marker <- unique(OSVZ_marker)
table(OSVZ_marker %in% macaque_annotation$gene_id)
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% macaque_annotation$gene_id]
OSVZ_marker <- macaque_annotation$gene_name[base::match(OSVZ_marker,table = macaque_annotation$gene_id)]
table(OSVZ_marker %in% rownames(macaque_multiome_Seurat))
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% rownames(macaque_multiome_Seurat)]
table(duplicated(OSVZ_marker))
OSVZ_marker <- unique(OSVZ_marker)
OSVZ_plot <- geneSetAveragePlot(genes = OSVZ_marker,object = macaque_multiome_Seurat,
                              object.class = 'seurat',plot.type = 'average',embedding = 'harmonyUMAP',assay = 'RNA',
                              aspectratio = 1,trim = c(0,0.99),color.palette = c("lightgrey",'red'),
                              print = FALSE,scaled = FALSE,plot.title = 'OSVZ bulk marker average plot')

#ISVZ
ISVZ_marker <- marker_list[['ISVZ']]
ISVZ_marker <- unique(ISVZ_marker)
table(ISVZ_marker %in% macaque_annotation$gene_id)
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% macaque_annotation$gene_id]
ISVZ_marker <- macaque_annotation$gene_name[base::match(ISVZ_marker,table = macaque_annotation$gene_id)]
table(ISVZ_marker %in% rownames(macaque_multiome_Seurat))
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% rownames(macaque_multiome_Seurat)]
table(duplicated(ISVZ_marker))
ISVZ_marker <- unique(ISVZ_marker)
ISVZ_plot <- geneSetAveragePlot(genes = ISVZ_marker,object = macaque_multiome_Seurat,
                              object.class = 'seurat',plot.type = 'average',embedding = 'harmonyUMAP',assay = 'RNA',
                              aspectratio = 1,trim = c(0,0.99),color.palette = c("lightgrey",'red'),
                              print = FALSE,scaled = FALSE,plot.title = 'ISVZ bulk marker average plot')

pdf(file = './res/step_22_fig_220330/macaque_multiome_Seurat_bulk_layer_marker_.pdf',width = 16,height = 8)
CPo_plot+CPi_plot+SP_plot+IZ_plot+OSVZ_plot+ISVZ_plot + plot_layout(ncol = 3)
dev.off()