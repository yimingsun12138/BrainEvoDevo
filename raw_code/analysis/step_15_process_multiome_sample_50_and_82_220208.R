#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: process multiome sample 50 and 82                               ##
## Data: 2022.02.08                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# load data ---------------------------------------------------------------
A50A_multiome_Seurat <- Read10X(data.dir = '/data/User/sunym/project/Brain/data/Multiome/A50A/outs/filtered_feature_bc_matrix/')
A50B_multiome_Seurat <- Read10X(data.dir = '/data/User/sunym/project/Brain/data/Multiome/A50B/outs/filtered_feature_bc_matrix/')
A82A_multiome_Seurat <- Read10X(data.dir = '/data/User/sunym/project/Brain/data/Multiome/A82A/outs/filtered_feature_bc_matrix/')
A82B_multiome_Seurat <- Read10X(data.dir = '/data/User/sunym/project/Brain/data/Multiome/A82B/outs/filtered_feature_bc_matrix/')

#create seurat object
A50A_multiome_Seurat <- CreateSeuratObject(counts = A50A_multiome_Seurat$`Gene Expression`,project = 'A50A',assay = 'RNA',min.cells = 0,min.features = 0)
A50B_multiome_Seurat <- CreateSeuratObject(counts = A50B_multiome_Seurat$`Gene Expression`,project = 'A50B',assay = 'RNA',min.cells = 0,min.features = 0)
A82A_multiome_Seurat <- CreateSeuratObject(counts = A82A_multiome_Seurat$`Gene Expression`,project = 'A82A',assay = 'RNA',min.cells = 0,min.features = 0)
A82B_multiome_Seurat <- CreateSeuratObject(counts = A82B_multiome_Seurat$`Gene Expression`,project = 'A82B',assay = 'RNA',min.cells = 0,min.features = 0)

# QC ----------------------------------------------------------------------
MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(A50A_multiome_Seurat))

A50A_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = A50A_multiome_Seurat,features = MT_gene_list)
A50B_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = A50B_multiome_Seurat,features = MT_gene_list)
A82A_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = A82A_multiome_Seurat,features = MT_gene_list)
A82B_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = A82B_multiome_Seurat,features = MT_gene_list)

#A50A
VlnPlot(A50A_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(A50A_multiome_Seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A50A_multiome_Seurat <- subset(A50A_multiome_Seurat,subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & nFeature_RNA > 200)

#A50B
VlnPlot(A50B_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(A50B_multiome_Seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A50B_multiome_Seurat <- subset(A50B_multiome_Seurat,subset = nFeature_RNA < 6000 & nCount_RNA < 15000 & nFeature_RNA > 200)

#A82A
VlnPlot(A82A_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(A82A_multiome_Seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A82A_multiome_Seurat <- subset(A82A_multiome_Seurat,subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & nFeature_RNA > 200)

#A82B
VlnPlot(A82B_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA','percent.mt'),ncol = 3)
FeatureScatter(A82B_multiome_Seurat,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A82B_multiome_Seurat <- subset(A82B_multiome_Seurat,subset = nFeature_RNA < 6000 & nCount_RNA < 30000 & nFeature_RNA > 200)

# merge -------------------------------------------------------------------
A50A_multiome_Seurat <- RenameCells(object = A50A_multiome_Seurat,new.names = paste('A50A',colnames(A50A_multiome_Seurat),sep = '_'))
A50B_multiome_Seurat <- RenameCells(object = A50B_multiome_Seurat,new.names = paste('A50B',colnames(A50B_multiome_Seurat),sep = '_'))
A82A_multiome_Seurat <- RenameCells(object = A82A_multiome_Seurat,new.names = paste('A82A',colnames(A82A_multiome_Seurat),sep = '_'))
A82B_multiome_Seurat <- RenameCells(object = A82B_multiome_Seurat,new.names = paste('A82B',colnames(A82B_multiome_Seurat),sep = '_'))

table(duplicated(c(colnames(A50A_multiome_Seurat),colnames(A50B_multiome_Seurat),colnames(A82A_multiome_Seurat),colnames(A82B_multiome_Seurat))))
table(rownames(A50A_multiome_Seurat) == rownames(A82A_multiome_Seurat))

macaque_multiome_Seurat <- cbind(A50A_multiome_Seurat@assays$RNA@counts,A50B_multiome_Seurat@assays$RNA@counts,A82A_multiome_Seurat@assays$RNA@counts,A82B_multiome_Seurat@assays$RNA@counts)
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
temp <- base::lapply(X = colnames(macaque_multiome_Seurat),FUN = function(x){
  return(strsplit(x = x,split = '_',fixed = TRUE)[[1]][1])
})

temp <- unlist(temp)
macaque_multiome_Seurat$donor <- temp

# process macaque multiome object -----------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

# label transfer ----------------------------------------------------------
temp <- macaque_multiome_Seurat
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')

feature_list <- SelectIntegrationFeatures(object.list = list(temp,macaque_RNA_seurat),verbose = TRUE)
temp <- RunPCA(object = temp,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)

anchors <- FindTransferAnchors(reference = macaque_RNA_seurat, query = temp, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = macaque_RNA_seurat$cell_type, dims = 1:30)
temp <- AddMetaData(temp,metadata = predictions)
temp@meta.data[1:3,]
macaque_multiome_Seurat@meta.data[colnames(temp),'predicted_label'] <- temp$predicted.id

DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)

# filter doublet ----------------------------------------------------------
#A50A
A50A_doublet <- read.csv(file = './res/step_15_fig_220208/A50A_scrublet_220208.txt')
A50A_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A50A_doublet$barcodes))
A50A_doublet$barcodes <- paste('A50A',as.character(A50A_doublet$barcodes),sep = '_')
table(duplicated(A50A_doublet$barcodes))
rownames(A50A_doublet) <- as.character(A50A_doublet$barcodes)
table(rownames(A50A_doublet) %in% colnames(macaque_multiome_Seurat))

#A50B
A50B_doublet <- read.csv(file = './res/step_15_fig_220208/A50B_scrublet_220208.txt')
A50B_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A50B_doublet$barcodes))
A50B_doublet$barcodes <- paste('A50B',as.character(A50B_doublet$barcodes),sep = '_')
table(duplicated(A50B_doublet$barcodes))
rownames(A50B_doublet) <- as.character(A50B_doublet$barcodes)
table(rownames(A50B_doublet) %in% colnames(macaque_multiome_Seurat))

#A82A
A82A_doublet <- read.csv(file = './res/step_15_fig_220208/A82A_scrublet_220208.txt')
A82A_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A82A_doublet$barcodes))
A82A_doublet$barcodes <- paste('A82A',as.character(A82A_doublet$barcodes),sep = '_')
table(duplicated(A82A_doublet$barcodes))
rownames(A82A_doublet) <- as.character(A82A_doublet$barcodes)
table(rownames(A82A_doublet) %in% colnames(macaque_multiome_Seurat))

#A82B
A82B_doublet <- read.csv(file = './res/step_15_fig_220208/A82B_scrublet_220208.txt')
A82B_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,as.character(A82B_doublet$barcodes))
A82B_doublet$barcodes <- paste('A82B',as.character(A82B_doublet$barcodes),sep = '_')
table(duplicated(A82B_doublet$barcodes))
rownames(A82B_doublet) <- as.character(A82B_doublet$barcodes)
table(rownames(A82B_doublet) %in% colnames(macaque_multiome_Seurat))

#doublet list
double_list <- rbind(A50A_doublet,A50B_doublet,A82A_doublet,A82B_doublet)
double_list <- double_list[,-1]
colnames(double_list) <- c('doublet_score','doublet')
double_list <- double_list[colnames(macaque_multiome_Seurat),]
macaque_multiome_Seurat@meta.data <- cbind(macaque_multiome_Seurat@meta.data,double_list)

p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'doublet',label = FALSE,cols = c('lightgrey','red'))
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)

p1+p2+p3+plot_layout(ncol = 3)

cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'doublet',split.by = 'seurat_clusters')
temp_col <- c('lightgrey','black')
names(temp_col) <- c('False','True')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = temp_col) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6) + 
  geom_hline(yintercept = 0.5,color = 'red')

#cluster 25 3 30 37 39 41 may be doublet
p1+p2+p3+plot_layout(ncol = 3)
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',features = c('nCount_RNA','nFeature_RNA'),pt.size = 0)
#filter cluster 30 37
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$seurat_clusters %in% c('30','37'))]
DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)

# donor specific filter ---------------------------------------------------
cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#cluster 22 24 25 26 29 32 40 5 9 may be donor sepcific
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'donor',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

#filter cluster 29
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$seurat_clusters %in% c('29'))]
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)


# InPSB filter ------------------------------------------------------------
temp <- macaque_multiome_Seurat@meta.data
temp$InPSB <- 'NO'
temp[temp$predicted_label == 'InPSB',"InPSB"] <- 'YES'

cell_proportion_table <- My_Cell_Proportion(meta_data = temp,group.by = 'InPSB',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=InPSB)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#cluster 12 17 28 33 41
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#filter 12 17 28 33
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$seurat_clusters %in% c('12','17','28','33'))]
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

# low quality filter ------------------------------------------------------
p1 <- VlnPlot(macaque_multiome_Seurat,features = c('nCount_RNA'),pt.size = 0)
p2 <- VlnPlot(macaque_multiome_Seurat,features = c('nFeature_RNA'),pt.size = 0)
p1+p2+plot_layout(ncol = 1)

#cluster 2 8 11 24 all In, not filter
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#cluster 4 5 16 all meaningful, not filter
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

FeaturePlot(macaque_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA'),label = TRUE,repel = TRUE)

p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#redo cluster to filter the doublet
macaque_multiome_Seurat <- FindNeighbors(macaque_multiome_Seurat,dims = 1:35,reduction = 'pca')
macaque_multiome_Seurat <- FindClusters(object = macaque_multiome_Seurat,resolution = 1.5)

p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

p1 <- VlnPlot(macaque_multiome_Seurat,features = c('nCount_RNA'),pt.size = 0)
p2 <- VlnPlot(macaque_multiome_Seurat,features = c('nFeature_RNA'),pt.size = 0)
p1+p2+plot_layout(ncol = 1)

#filter cluster 33
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$seurat_clusters %in% c('33'))]
p1 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

# redo process ------------------------------------------------------------
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("donor","doublet_score","doublet")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 200)
dim(macaque_multiome_Seurat)
gc()

## preprocess --------------------------------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 35,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

## label transfer ----------------------------------------------------------
temp <- macaque_multiome_Seurat
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')

feature_list <- SelectIntegrationFeatures(object.list = list(temp,macaque_RNA_seurat),verbose = TRUE)
temp <- RunPCA(object = temp,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)

anchors <- FindTransferAnchors(reference = macaque_RNA_seurat, query = temp, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = macaque_RNA_seurat$cell_type, dims = 1:30)
temp <- AddMetaData(temp,metadata = predictions)
macaque_multiome_Seurat@meta.data[colnames(temp),'predicted_label'] <- temp$predicted.id

my_send_sms('label transfer done!')
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)
#cluster 11 seems doublet

## filter donor specific ---------------------------------------------------
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

#cluster 10 16 17 20 27 34 36 may be donor specific
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

## filter doublet ----------------------------------------------------------
cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'doublet',split.by = 'seurat_clusters')
temp_col <- c('lightgrey','black')
names(temp_col) <- c('False','True')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = temp_col) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6) + 
  geom_hline(yintercept = 0.5,color = 'red')

#cluster 11 31 34 37 4 may be doublet
DimPlot(macaque_multiome_Seurat,cells.highlight = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$seurat_clusters %in% c('11')])
DimPlot(macaque_multiome_Seurat,cells.highlight = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$seurat_clusters %in% c('31')])
DimPlot(macaque_multiome_Seurat,cells.highlight = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$seurat_clusters %in% c('34')])
DimPlot(macaque_multiome_Seurat,cells.highlight = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$seurat_clusters %in% c('37')])
DimPlot(macaque_multiome_Seurat,cells.highlight = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$seurat_clusters %in% c('4')])

#filter cluster 11 31 34
cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$seurat_clusters %in% c('11','31','34'))]

p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

## low quality filter ------------------------------------------------------
FeaturePlot(macaque_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA'),label = TRUE)
p1 <- VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = 'nCount_RNA')
p2 <- VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = 'nFeature_RNA')
p1+p2+plot_layout(ncol = 2)

# redo filter -------------------------------------------------------------
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("donor","doublet_score","doublet")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 200)
dim(macaque_multiome_Seurat)
gc()

## preprocess --------------------------------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 24,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

## label transfer ----------------------------------------------------------
temp <- macaque_multiome_Seurat
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')

feature_list <- SelectIntegrationFeatures(object.list = list(temp,macaque_RNA_seurat),verbose = TRUE)
temp <- RunPCA(object = temp,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)

anchors <- FindTransferAnchors(reference = macaque_RNA_seurat, query = temp, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = macaque_RNA_seurat$cell_type, dims = 1:30)
temp <- AddMetaData(temp,metadata = predictions)
macaque_multiome_Seurat@meta.data[colnames(temp),'predicted_label'] <- temp$predicted.id

my_send_sms('label transfer done!')
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

## donor specific filter ---------------------------------------------------
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

## doublet filter ----------------------------------------------------------
FeaturePlot(macaque_multiome_Seurat,features = c('doublet_score'),label = TRUE)
VlnPlot(object = macaque_multiome_Seurat,features = 'doublet_score',group.by = 'seurat_clusters',pt.size = 0)

#cluster 1 25 33 may be doublet
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

DimPlot(macaque_multiome_Seurat,cells.highlight = colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$seurat_clusters == '1'])

#filter cluster 1
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$seurat_clusters %in% c('1'))]

p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)


## low quality filter ------------------------------------------------------
FeaturePlot(macaque_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA'),label = TRUE)

# final filter ------------------------------------------------------------
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("donor","doublet_score","doublet")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 200)
dim(macaque_multiome_Seurat)
gc()

## preprocess --------------------------------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 25,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)
## label transfer ----------------------------------------------------------
temp <- macaque_multiome_Seurat
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')

feature_list <- SelectIntegrationFeatures(object.list = list(temp,macaque_RNA_seurat),verbose = TRUE)
temp <- RunPCA(object = temp,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)

anchors <- FindTransferAnchors(reference = macaque_RNA_seurat, query = temp, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = macaque_RNA_seurat$cell_type, dims = 1:30)
temp <- AddMetaData(temp,metadata = predictions)
macaque_multiome_Seurat@meta.data[colnames(temp),'predicted_label'] <- temp$predicted.id

my_send_sms('label transfer done!')

p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 24 and cluster 25 so weird
#cluster 34 26 is weird
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = c('nCount_RNA'))
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = c('nFeature_RNA'))
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = c('doublet_score'))

cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

table(macaque_multiome_Seurat[,macaque_multiome_Seurat$seurat_clusters == '34']$predicted_label)
table(macaque_multiome_Seurat[,macaque_multiome_Seurat$seurat_clusters == '26']$predicted_label)

#filter cluster 34 (doublet)
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!macaque_multiome_Seurat$seurat_clusters %in% c('34')]

p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 24 25 marker
marker_list <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('24'),group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(marker_list)[1:12],pt.size = 0,group.by = 'seurat_clusters')
table(marker_list$pct.1 > 0.5 & marker_list$pct.2 < 0.3)

marker_list <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = c('25'),group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE,verbose = TRUE)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(marker_list)[1:12],pt.size = 0,group.by = 'seurat_clusters')
table(marker_list$pct.1 > 0.5 & marker_list$pct.2 < 0.3)

## marker gene express -----------------------------------------------------
dotplot_matrix <- my_dotplot(macaque_multiome_Seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Ast=c('AQP4','APOE'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2'),
                                             PSB=c('MEIS2','ETV1')),
                             group.by = 'seurat_clusters', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

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

#filter cluster 24 25
macaque_multiome_Seurat <- macaque_multiome_Seurat[,!(macaque_multiome_Seurat$seurat_clusters %in% c('24','25'))]
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+p3+plot_layout(ncol = 3)

# final final filter ------------------------------------------------------
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("donor","doublet_score","doublet")]
macaque_multiome_Seurat <- macaque_multiome_Seurat@assays$RNA@counts
macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 200)
dim(macaque_multiome_Seurat)
gc()

## preprocess --------------------------------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 34,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

DimPlot(object = macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
DimPlot(object = macaque_multiome_Seurat,cells.highlight = colnames(macaque_multiome_Seurat[,macaque_multiome_Seurat$seurat_clusters == '28']))

#cluster 28 29 32 seems weird

## donor specific filter ---------------------------------------------------
cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

## low quality filter ------------------------------------------------------
FeaturePlot(macaque_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA'),label = TRUE)
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = 'nCount_RNA') + 
  VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = 'nFeature_RNA') + 
  plot_layout(ncol = 1)

## doublet filter ----------------------------------------------------------
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',features = 'doublet_score',pt.size = 0)

## classic marker filter ---------------------------------------------------
dotplot_matrix <- my_dotplot(macaque_multiome_Seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Ast=c('AQP4','APOE'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2'),
                                             PSB=c('MEIS2','ETV1')),
                             group.by = 'seurat_clusters', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

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

#cluster 28 may by InPSB
marker_list <- FindMarkers(object = macaque_multiome_Seurat,ident.1 = '28',group.by = 'seurat_clusters',assay = 'RNA',slot = 'data',test.use = 'bimod',verbose = TRUE,only.pos = TRUE)
VlnPlot(object = macaque_multiome_Seurat,features = rownames(marker_list)[1:12],pt.size = 0)

## label transfer ----------------------------------------------------------
temp <- macaque_multiome_Seurat
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')

feature_list <- SelectIntegrationFeatures(object.list = list(temp,macaque_RNA_seurat),verbose = TRUE)
temp <- RunPCA(object = temp,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'RNA',features = feature_list,npcs = 50,verbose = TRUE)

anchors <- FindTransferAnchors(reference = macaque_RNA_seurat, query = temp, dims = 1:30, reference.reduction = "pca")
predictions <- TransferData(anchorset = anchors, refdata = macaque_RNA_seurat$cell_type, dims = 1:30)
temp <- AddMetaData(temp,metadata = predictions)
macaque_multiome_Seurat@meta.data[colnames(temp),'predicted_label'] <- temp$predicted.id

my_send_sms('label transfer done!')
p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'predicted_label',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

# save data ---------------------------------------------------------------
saveRDS(macaque_multiome_Seurat@assays$RNA@counts,file = './res/step_15_fig_220208/macaque_multiome_counts.rds')
saveRDS(macaque_multiome_Seurat@meta.data,file = './res/step_15_fig_220208/macaque_multiome_meta_data.rds')

# multiome ATAC_seq filter ------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './res/step_15_fig_220208/macaque_multiome_counts.rds')
meta_data <- readRDS(file = './res/step_15_fig_220208/macaque_multiome_meta_data.rds')
meta_data <- meta_data[,c("donor","doublet_score","doublet")]

macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
macaque_multiome_Seurat <- RenameCells(object = macaque_multiome_Seurat,new.names = sub(pattern = '_',replacement = '#',fixed = TRUE,x = colnames(macaque_multiome_Seurat)))

#ATAC_seq filter
meta_data <- readRDS(file = './ArchR/res/step_5_fig_220215/macaque_multiome_ArchR_meta_data.rds')
table(rownames(meta_data) %in% colnames(macaque_multiome_Seurat))
macaque_multiome_Seurat <- macaque_multiome_Seurat[,rownames(meta_data)]
dim(macaque_multiome_Seurat)
gc()

## preprocess --------------------------------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 32,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

## donor specific filter ---------------------------------------------------
cell_proportion_table <- My_Cell_Proportion(meta_data = macaque_multiome_Seurat@meta.data,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(cell_proportion_table,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

p1 <- DimPlot(macaque_multiome_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

## low quality filter ------------------------------------------------------
FeaturePlot(macaque_multiome_Seurat,features = c('nCount_RNA','nFeature_RNA'),label = TRUE)
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = 'nCount_RNA') + 
  VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',pt.size = 0,features = 'nFeature_RNA') + 
  plot_layout(ncol = 1)

## doublet filter ----------------------------------------------------------
VlnPlot(macaque_multiome_Seurat,group.by = 'seurat_clusters',features = 'doublet_score',pt.size = 0)


## combined with scRNA_seq -------------------------------------------------
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('InPSB'))]
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Astrocyte',"cell_type"] <- 'Ependymal'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#merge
gene_list <- dplyr::intersect(rownames(macaque_RNA_seurat@assays$RNA@counts),rownames(macaque_multiome_Seurat@assays$RNA@counts))
macaque_integrated_Seurat <- cbind(macaque_RNA_seurat@assays$RNA@counts[gene_list,],macaque_multiome_Seurat@assays$RNA@counts[gene_list,])
macaque_integrated_Seurat <- CreateSeuratObject(counts = macaque_integrated_Seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)

#add meta data
macaque_integrated_Seurat$donor <- NA
macaque_integrated_Seurat$batch <- NA
macaque_integrated_Seurat$tech <- NA
macaque_integrated_Seurat$cell_type <- NA

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"donor"] <- as.character(macaque_RNA_seurat$donor)
macaque_integrated_Seurat@meta.data[colnames(macaque_multiome_Seurat),"donor"] <- as.character(macaque_multiome_Seurat$donor)

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"batch"] <- as.character(macaque_RNA_seurat$batch)
macaque_integrated_Seurat@meta.data[colnames(macaque_multiome_Seurat),"batch"] <- '220225'

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"tech"] <- 'scRNA'
macaque_integrated_Seurat@meta.data[colnames(macaque_multiome_Seurat),"tech"] <- 'multiome'

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)

#process
macaque_integrated_Seurat <- my_process_seurat(object = macaque_integrated_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor','batch','tech'),npcs = 50,preprocess = TRUE)
macaque_integrated_Seurat <- my_process_seurat(object = macaque_integrated_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

p1 <- DimPlot(macaque_integrated_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_integrated_Seurat,group.by = 'batch',label = FALSE)
p3 <- DimPlot(macaque_integrated_Seurat,group.by = 'tech',label = FALSE)
p4 <- DimPlot(macaque_integrated_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+p3+p4+plot_layout(ncol = 2)

#label transfer
ref <- macaque_integrated_Seurat[,macaque_integrated_Seurat$tech == 'scRNA']
query <- macaque_integrated_Seurat[,macaque_integrated_Seurat$tech == 'multiome']
anchors <- my_FindTransferAnchors(reference = ref,query = query,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'RNA',query_assay = 'RNA',l2.norm = FALSE,dims = 1:30,features = NULL,verbose = TRUE)

predictions <- TransferData(anchorset = anchors,refdata = ref$cell_type,weight.reduction = 'pcaproject',l2.norm = FALSE,verbose = TRUE,dims = 1:30)
query <- AddMetaData(object = query,metadata = predictions)
DimPlot(query,group.by = 'predicted.id',reduction = 'umap',label = TRUE,repel = TRUE)
hist(query$prediction.score.max)

macaque_integrated_Seurat@meta.data[colnames(query),"cell_type"] <- query$predicted.id

pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_merge_dim_30_label_transfer_dimplot.pdf',width = 16,height = 8)
DimPlot(macaque_integrated_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'tech')
dev.off()

p1 <- DimPlot(macaque_integrated_Seurat,group.by = 'batch',label = FALSE)
p2 <- DimPlot(macaque_integrated_Seurat,group.by = 'tech',label = FALSE)
p3 <- DimPlot(macaque_integrated_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p4 <- DimPlot(macaque_integrated_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_merge_dim_30_dimplot.pdf',width = 14,height = 10)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()

#marker plot
dotplot_matrix <- my_dotplot(macaque_integrated_Seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('AQP4','FOXJ1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))

pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_merge_dim_30_marker_dotplot.pdf',width = 16,height = 8)
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

dotplot_matrix <- my_dotplot(macaque_integrated_Seurat[,macaque_integrated_Seurat$tech == 'multiome'],assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('AQP4','FOXJ1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))

pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_merge_dim_30_only_multiome_marker_dotplot.pdf',width = 16,height = 8)
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

#save data
meta_data <- macaque_integrated_Seurat@meta.data
saveRDS(meta_data,file = './res/step_15_fig_220208/macaque_integrated_seurat_merge_dim_30_label_transfer_metadata.rds')

# harmony integration with scRNA_seq --------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './res/step_15_fig_220208/macaque_multiome_counts.rds')
meta_data <- readRDS(file = './res/step_15_fig_220208/macaque_multiome_meta_data.rds')
meta_data <- meta_data[,c("donor","doublet_score","doublet")]

macaque_multiome_Seurat <- CreateSeuratObject(counts = macaque_multiome_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
macaque_multiome_Seurat <- RenameCells(object = macaque_multiome_Seurat,new.names = sub(pattern = '_',replacement = '#',fixed = TRUE,x = colnames(macaque_multiome_Seurat)))

#ATAC_seq filter
meta_data <- readRDS(file = './ArchR/res/step_5_fig_220215/macaque_multiome_ArchR_meta_data.rds')
table(rownames(meta_data) %in% colnames(macaque_multiome_Seurat))
macaque_multiome_Seurat <- macaque_multiome_Seurat[,rownames(meta_data)]
dim(macaque_multiome_Seurat)
gc()

## preprocess --------------------------------------------------------------
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_Seurat <- my_process_seurat(object = macaque_multiome_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 32,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

## combined with scRNA_seq -------------------------------------------------
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('InPSB'))]
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Astrocyte',"cell_type"] <- 'Ependymal'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#check clusters 28
VlnPlot(object = macaque_RNA_seurat,features = c('nCount_RNA'),group.by = 'seurat_clusters',pt.size = 0)
VlnPlot(object = macaque_RNA_seurat,features = c('nFeature_RNA'),group.by = 'seurat_clusters',pt.size = 0)

dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('AQP4','FOXJ1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'seurat_clusters', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

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

#filter cluster 28
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('28'))]

#harmony
macaque_multiome_Seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

table(VariableFeatures(macaque_RNA_seurat) %in% rownames(macaque_multiome_Seurat))
gene_list <- VariableFeatures(macaque_RNA_seurat)
dim_used <- 33
#try 32 33

macaque_integrated_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_Seurat),
                                                    assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                          multiome=c('nCount_RNA','donor')),
                                                    npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                    UMAP_dim = dim_used,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                                    yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integrated_Seurat$cell_type <- NA
macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$cell_type

p1 <- DimPlot(object = macaque_integrated_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integrated_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integrated_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#label transfer
ref <- macaque_integrated_Seurat[,macaque_integrated_Seurat$dataset == 'scRNA']
query <- macaque_integrated_Seurat[,macaque_integrated_Seurat$dataset == 'multiome']

anchors <- my_FindTransferAnchors(reference = ref,query = query,ref_reduction = 'pca',query_reduction = 'pca',
                                  ref_assay = 'integration',query_assay = 'integration',l2.norm = FALSE,dims = 1:dim_used,
                                  features = NULL,verbose = TRUE)

predictions <- TransferData(anchorset = anchors,refdata = ref$cell_type,weight.reduction = 'pcaproject',
                            l2.norm = FALSE,dims = 1:dim_used,verbose = TRUE)

query <- AddMetaData(object = query,metadata = predictions)
DimPlot(query,group.by = 'predicted.id',label = TRUE,repel = TRUE)
hist(query$prediction.score.max)

macaque_integrated_Seurat@meta.data[colnames(query),"cell_type"] <- query$predicted.id

#recreate macaque_integration_Seurat
temp <- macaque_integrated_Seurat
macaque_integrated_Seurat <- CreateSeuratObject(counts = temp@assays$integration@counts,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)

#add metadata
macaque_integrated_Seurat$donor <- NA
macaque_integrated_Seurat$batch <- NA
macaque_integrated_Seurat$tech <- NA
macaque_integrated_Seurat$cell_type <- NA

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"donor"] <- as.character(macaque_RNA_seurat$donor)
macaque_integrated_Seurat@meta.data[colnames(macaque_multiome_Seurat),"donor"] <- as.character(macaque_multiome_Seurat$donor)

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"batch"] <- as.character(macaque_RNA_seurat$batch)
macaque_integrated_Seurat@meta.data[colnames(macaque_multiome_Seurat),"batch"] <- '220225'

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"tech"] <- 'scRNA'
macaque_integrated_Seurat@meta.data[colnames(macaque_multiome_Seurat),"tech"] <- 'multiome'

macaque_integrated_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)
macaque_integrated_Seurat@meta.data[colnames(query),"cell_type"] <- as.character(query$predicted.id)

#process
macaque_integrated_Seurat <- my_process_seurat(object = macaque_integrated_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor','batch','tech'),npcs = 50,preprocess = TRUE)
macaque_integrated_Seurat <- my_process_seurat(object = macaque_integrated_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 30,resolution = c(0.5,1,1.5),group.by = 'seurat_clusters',label = TRUE)

p1 <- DimPlot(macaque_integrated_Seurat,group.by = 'donor',label = FALSE)
p2 <- DimPlot(macaque_integrated_Seurat,group.by = 'batch',label = FALSE)
p3 <- DimPlot(macaque_integrated_Seurat,group.by = 'tech',label = FALSE)
p4 <- DimPlot(macaque_integrated_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+p3+p4+plot_layout(ncol = 2)

macaque_integrated_Seurat[['harmonyPCA']] <- temp@reductions$pca
macaque_integrated_Seurat@reductions$harmonyUMAP <- temp@reductions$umap

#plot
pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_harmony_dim_33_label_transfer_dimplot.pdf',width = 16,height = 8)
DimPlot(macaque_integrated_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,split.by = 'tech',reduction = 'harmonyUMAP')
dev.off()

p1 <- DimPlot(macaque_integrated_Seurat,group.by = 'batch',label = FALSE,reduction = 'harmonyUMAP')
p2 <- DimPlot(macaque_integrated_Seurat,group.by = 'tech',label = FALSE,reduction = 'harmonyUMAP')
p3 <- DimPlot(macaque_integrated_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p4 <- DimPlot(macaque_integrated_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_harmony_dim_33_dimplot.pdf',width = 14,height = 10)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()

#marker plot
dotplot_matrix <- my_dotplot(macaque_integrated_Seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('AQP4','FOXJ1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))

pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_harmony_dim_33_marker_dotplot.pdf',width = 16,height = 8)
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

dotplot_matrix <- my_dotplot(macaque_integrated_Seurat[,macaque_integrated_Seurat$tech == 'multiome'],assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('AQP4','FOXJ1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))

pdf(file = './res/step_15_fig_220208/macaque_integrated_seurat_harmony_dim_33_only_multiome_marker_dotplot.pdf',width = 16,height = 8)
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

#save data
saveRDS(macaque_integrated_Seurat,file = './processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_220226.rds')
