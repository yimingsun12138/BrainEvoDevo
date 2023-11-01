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

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
macaque_integration_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_multiome_Seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']

DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

# data normalize ----------------------------------------------------------
macaque_multiome_cpm <- cpm(counts = macaque_multiome_Seurat@assays$RNA@counts)
macaque_multiome_norm.dat <- log2(macaque_multiome_cpm + 1)
macaque_multiome_norm.dat <- Matrix(macaque_multiome_norm.dat, sparse = TRUE)

# set cluster parameter ---------------------------------------------------
de.param <- de_param(padj.th = 0.05, 
                     lfc.th = 1, 
                     low.th = 1, 
                     q1.th = 0.3,
                     q2.th = NULL,
                     q.diff.th = 0.7, 
                     de.score.th = 150,
                     min.cells = 4)

# Dimension Filtering -----------------------------------------------------
gene.counts <- colSums(macaque_multiome_norm.dat > 0)
rm.eigen <- matrix(log2(gene.counts), ncol = 1)
row.names(rm.eigen) <- names(gene.counts)
colnames(rm.eigen) <- "log2GeneCounts"

# Clustering --------------------------------------------------------------

## Coarse-level clustering -------------------------------------------------
strict.param <- de_param(de.score.th = 500)
onestep.result <- onestep_clust(macaque_multiome_norm.dat, 
                                select.cells = colnames(macaque_multiome_norm.dat), 
                                dim.method = "pca", 
                                de.param = strict.param, 
                                rm.eigen = rm.eigen)

display.result <- display_cl(onestep.result$cl, macaque_multiome_norm.dat, plot = TRUE, de.param = de.param)

#see the one_step cluster results
macaque_multiome_Seurat$temp <- onestep.result$cl[colnames(macaque_multiome_Seurat)]
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'temp',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p1+p2+p3+plot_layout(ncol = 3)

## iterative cluster -------------------------------------------------------
iter.result <- iter_clust(macaque_multiome_norm.dat, 
                          select.cells = colnames(macaque_multiome_Seurat), 
                          dim.method = "pca", 
                          de.param = de.param, 
                          rm.eigen = rm.eigen,
                          result = onestep.result)

macaque_multiome_Seurat$temp <- iter.result$cl[colnames(macaque_multiome_Seurat)]
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'temp',label = FALSE,repel = TRUE,reduction = 'harmonyUMAP') + theme(legend.position = 'none')
p1+p2+p3+plot_layout(ncol = 3)

scibet::Confusion_heatmap(ori = macaque_multiome_Seurat$cell_type,prd = macaque_multiome_Seurat$temp)

# merge cluster -----------------------------------------------------------
macaque_multiome_rd.data <- t(macaque_multiome_norm.dat[iter.result$markers, colnames(macaque_multiome_Seurat)])

merge.param <- de_param(de.score.th = 300) # The original value was 150.

merge.result <- merge_cl(macaque_multiome_norm.dat, 
                         cl = iter.result$cl, 
                         rd.dat = macaque_multiome_rd.data,
                         de.param = merge.param)

macaque_multiome_Seurat$temp <- merge.result$cl[colnames(macaque_multiome_Seurat)]
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'temp',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

pdf(file = './res/step_20_fig_220315/macaque_multiome_Seurat_merged_cluster.pdf',width = 21,height = 7)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

# Bootstrapping and Consensus ---------------------------------------------
set.seed(12345)
result <- run_consensus_clust(macaque_multiome_norm.dat, 
                              select.cells = colnames(macaque_multiome_Seurat),
                              niter = 100, 
                              de.param = de.param, 
                              rm.eigen = rm.eigen, 
                              dim.method = "pca", 
                              output_dir = "./res/step_20_fig_220315/subsample_PCA",
                              mc.cores = 4)
my_send_sms('cluster done!')

macaque_multiome_Seurat$temp <- result$cl.result$cl[colnames(macaque_multiome_Seurat)]
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'temp',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP') + theme(legend.position = 'none')
pdf(file = './res/step_20_fig_220315/macaque_multiome_Seurat_consensus_cluster.pdf',width = 21,height = 7)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()
hist(table(result$cl.result$cl),breaks = 100)

#fucking dead end!