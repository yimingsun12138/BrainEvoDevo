#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: annotate mouse multiome RNA data                                ##
## Data: 2022.09.12                                                                ##
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
library(SingleCellExperiment)
library(sctransform)
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
library(patchwork)
library(networkD3)
library(htmlwidgets)
library(circlize)
library(harmony)
library(ArchR)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# load data ---------------------------------------------------------------

#express matrix
E145_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E145_1/outs/filtered_feature_bc_matrix/')
E145_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E145_2/outs/filtered_feature_bc_matrix/')
E155_1 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_1/outs/filtered_feature_bc_matrix/')
E155_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220706/E155_2/outs/filtered_feature_bc_matrix/')
E155_3 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E155_3/outs/filtered_feature_bc_matrix/')
E165_2 <- Seurat::Read10X(data.dir = './data/mouse_multiome/220907/E165_2/outs/filtered_feature_bc_matrix/')

E145_1 <- E145_1$`Gene Expression`
E145_2 <- E145_2$`Gene Expression`
E155_1 <- E155_1$`Gene Expression`
E155_2 <- E155_2$`Gene Expression`
E155_3 <- E155_3$`Gene Expression`
E165_2 <- E165_2$`Gene Expression`

#rename express matrix
colnames(E145_1) <- paste('E145_1#',colnames(E145_1),sep = '')
colnames(E145_2) <- paste('E145_2#',colnames(E145_2),sep = '')
colnames(E155_1) <- paste('E155_1#',colnames(E155_1),sep = '')
colnames(E155_2) <- paste('E155_2#',colnames(E155_2),sep = '')
colnames(E155_3) <- paste('E155_3#',colnames(E155_3),sep = '')
colnames(E165_2) <- paste('E165_2#',colnames(E165_2),sep = '')

E145_1 <- CreateSeuratObject(counts = E145_1,project = 'E145_1',assay = 'RNA',min.cells = 0,min.features = 0)
E145_2 <- CreateSeuratObject(counts = E145_2,project = 'E145_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_1 <- CreateSeuratObject(counts = E155_1,project = 'E155_1',assay = 'RNA',min.cells = 0,min.features = 0)
E155_2 <- CreateSeuratObject(counts = E155_2,project = 'E155_2',assay = 'RNA',min.cells = 0,min.features = 0)
E155_3 <- CreateSeuratObject(counts = E155_3,project = 'E155_3',assay = 'RNA',min.cells = 0,min.features = 0)
E165_2 <- CreateSeuratObject(counts = E165_2,project = 'E165_2',assay = 'RNA',min.cells = 0,min.features = 0)

#merge
MT_gene_list <- rownames(E145_1@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(E145_1@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat <- cbind(E145_1@assays$RNA@counts,E145_2@assays$RNA@counts,E155_1@assays$RNA@counts,E155_2@assays$RNA@counts,E155_3@assays$RNA@counts,E165_2@assays$RNA@counts)
mouse_multiome_Seurat <- CreateSeuratObject(counts = mouse_multiome_Seurat,project = 'mouse',assay = 'RNA',min.cells = 0,min.features = 0)
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene_list)

#filter cells
meta_data <- readRDS(file = './ArchR/res/step_13_fig_220910/mouse_multiome_Seurat_meta_data.rds')