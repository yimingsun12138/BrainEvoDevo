#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: semi_fig_220714                                                 ##
## Data: 2022.07.14                                                                ##
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

# macaque multiome dimplot ------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_multiome_Seurat_220712.rds')

p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))

pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_dimplot.pdf',width = 24,height = 8)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

#IP marker
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('EOMES')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('PPP1R17')) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_IP_marker_featureplot.pdf',width = 16,height = 12)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()

#Ex-1 marker
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('NRP1')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('SLC17A6')) + theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('DPY19L1')) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_Ex_2_marker_featureplot.pdf',width = 18,height = 10)
p1+p2+p3+p4+p5+plot_layout(ncol = 3)
dev.off()

#Ex-2 marker
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('HMCN1')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('CCBE1')) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_Ex_2_marker_featureplot.pdf',width = 16,height = 12)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()

#Ex-3_1 marker
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('KCNH5')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('VWC2L')) + theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('RORB')) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_Ex_3_1_marker_featureplot.pdf',width = 18,height = 10)
p1+p2+p3+p4+p5+plot_layout(ncol = 3)
dev.off()

#Ex-3_2 marker
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('ADAMTS3')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('IL1RAPL2')) + theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('POSTN')) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_Ex_3_2_marker_featureplot.pdf',width = 18,height = 10)
p1+p2+p3+p4+p5+plot_layout(ncol = 3)
dev.off()

#Ex-4_1 marker
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('FOXP2')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('RGS6')) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_Ex_4_1_marker_featureplot.pdf',width = 16,height = 12)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()

#Ex-4_2 marker
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('TMEM178A')) + theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('SLIT3')) + theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('CLSTN2')) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_Ex_4_2_marker_featureplot.pdf',width = 18,height = 10)
p1+p2+p3+p4+p5+plot_layout(ncol = 3)
dev.off()

# macaque ATAC clustering -------------------------------------------------
macaque_multiome_ArchR <- loadArchRProject(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/')
meta_data <- readRDS(file = './res/step_41_fig_220626/macaque_multiome_ArchR_LSI_clustering_meta_data.rds')
p1 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),],group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = macaque_multiome_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),],group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = meta_data[rownames(macaque_multiome_ArchR@cellColData),],group.by = "ATAC_snn_res.0.7",label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p4 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = meta_data[rownames(macaque_multiome_ArchR@cellColData),],group.by = "ATAC_snn_res.1.7",label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
pdf(file = './res/step_47_fig_220714/macaque_multiome_ArchR_dimplot.pdf',width = 16,height = 12)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()

# macaque integration Seurat dimplot --------------------------------------
macaque_integration_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_integration_Seurat_220713.rds')

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'tech',label = FALSE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
pdf(file = './res/step_47_fig_220714/macaque_integration_Seurat_dimplot.pdf',width = 24,height = 8)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

pdf(file = './res/step_47_fig_220714/macaque_integration_Seurat_cluster_cell_type_confusion_heatmap.pdf',width = 10,height = 7.5)
scibet::Confusion_heatmap(ori = macaque_integration_Seurat$seurat_clusters,prd = macaque_integration_Seurat$sub_cell_type) + theme(aspect.ratio = 19/27) + 
  theme(axis.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + xlab('Clusters') + ylab('cell type')
dev.off()

# mouse multiome RNA sequencing depth -------------------------------------
mouse_multiome_Seurat <- readRDS(file = './res/step_42_fig_220707/raw_mouse_multiome_Seurat_220710.rds')
p1 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('nCount_RNA'),max.cutoff = 'q99') + theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = mouse_multiome_Seurat,features = c('nFeature_RNA'),max.cutoff = 'q99') + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/mouse_multiome_Seurat_sequencing_depth_featureplot.pdf',width = 6,height = 8)
p1+p2+plot_layout(ncol = 1)
dev.off()

pdf(file = './res/step_47_fig_220714/mouse_multiome_Seurat_sequencing_depth_vlnplot.pdf',width = 8,height = 8)
VlnPlot(object = mouse_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA'),pt.size = 0,group.by = 'RNA_snn_res.0.5',ncol = 1)
dev.off()

# macaque multiome data UMAP and tsne -------------------------------------
macaque_multiome_Seurat <- RunTSNE(object = macaque_multiome_Seurat,dims = 1:31)
p1 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_multiome_Seurat,reduction = 'tsne',group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_dimplot_UMAP_vs_tsne.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

# mouse filted data marker featureplot ------------------------------------
meta_data <- readRDS(file = './res/step_42_fig_220707/mouse_multiome_data_rough_annotate_on_RNA_meta_data.rds')
mouse_multiome_Seurat <- mouse_multiome_Seurat[,rownames(meta_data)]
mouse_multiome_Seurat <- my_process_seurat(object = mouse_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
ndims <- 22
mouse_multiome_Seurat <- RunUMAP(object = mouse_multiome_Seurat,dims = 1:ndims,reduction = 'pca')
DimPlot(object = mouse_multiome_Seurat)

pdf(file = './res/step_47_fig_220714/macaque_multiome_Seurat_classic_marker_featureplot.pdf',width = 16,height = 12)
FeaturePlot(object = mouse_multiome_Seurat,features = c('Top2a','Mki67','Hes5','Pax6','Eomes','Ppp1r17','Neurod2','Neurod6','Gad2'))
dev.off()

# macaque Greenleaf alignment ---------------------------------------------
Brain_RNA_Seurat <- readRDS(file = './processed_data/220711_summary/macaque_integration_Seurat_Greenleaf_RNA_Seurat_harmony_integration_Seurat.rds')
pdf(file = './res/step_47_fig_220714/macaque_Greenleaf_data_dimplot.pdf',width = 12,height = 6)
DimPlot(object = Brain_RNA_Seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,split.by = 'dataset') + theme(aspect.ratio = 1)
dev.off()

Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/220711_summary/Greenleaf_RNA_Seurat_220714.rds')

p1 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
p2 <- DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'macaque_cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP') + theme(aspect.ratio = 1) + scale_color_manual(values = as.character(c(ArchRPalettes$stallion,ArchRPalettes$bear)))
pdf(file = './res/step_47_fig_220714/Greenleaf_RNA_Seurat_dimplot_label_transfered.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()