#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: merge sample from 68 and 84                                     ##
## Data: 2021.06.10                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

#load data
A68A_RNA <- Read10X(data.dir = './data/scRNA-seq/200919/A-68A/outs/filtered_feature_bc_matrix/')
A68B_RNA <- Read10X(data.dir = './data/scRNA-seq/200919/A-68B/outs/filtered_feature_bc_matrix/')
A84B_RNA <- Read10X(data.dir = './data/scRNA-seq/200919/A-84B/outs/filtered_feature_bc_matrix/')
A84C_RNA <- Read10X(data.dir = './data/scRNA-seq/200919/A-84C/outs/filtered_feature_bc_matrix/')

#create seurat object
A68A_RNA <- CreateSeuratObject(counts = A68A_RNA,project = 'A68A',min.cells = 0,min.features = 200)
A68B_RNA <- CreateSeuratObject(counts = A68B_RNA,project = 'A68B',min.cells = 0,min.features = 200)
A84B_RNA <- CreateSeuratObject(counts = A84B_RNA,project = 'A84B',min.cells = 0,min.features = 200)
A84C_RNA <- CreateSeuratObject(counts = A84C_RNA,project = 'A84C',min.cells = 0,min.features = 200)
table(rownames(A68A_RNA) == rownames(A68B_RNA))

#QC
MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(A68A_RNA))

A68A_RNA[['percent.mt']] <- PercentageFeatureSet(object = A68A_RNA,features = MT_gene_list)
A68B_RNA[['percent.mt']] <- PercentageFeatureSet(object = A68B_RNA,features = MT_gene_list)
A84B_RNA[['percent.mt']] <- PercentageFeatureSet(object = A84B_RNA,features = MT_gene_list)
A84C_RNA[['percent.mt']] <- PercentageFeatureSet(object = A84C_RNA,features = MT_gene_list)

VlnPlot(A68A_RNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(A68A_RNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A68A_RNA <- subset(A68A_RNA,subset = nFeature_RNA < 5000 & nCount_RNA < 15000)

VlnPlot(A68B_RNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(A68B_RNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A68B_RNA <- subset(A68B_RNA,subset = nFeature_RNA < 5000 & nCount_RNA < 15000)

VlnPlot(A84B_RNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(A84B_RNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A84B_RNA <- subset(A84B_RNA,subset = nFeature_RNA < 5000 & nCount_RNA < 15000)

VlnPlot(A84C_RNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(A84C_RNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A84C_RNA <- subset(A84C_RNA,subset = nFeature_RNA < 5000 & nCount_RNA < 15000)

#merge
A68A_RNA <- RenameCells(A68A_RNA,new.names = paste('A68A',colnames(A68A_RNA),sep = '_'))
A68B_RNA <- RenameCells(A68B_RNA,new.names = paste('A68B',colnames(A68B_RNA),sep = '_'))
A84B_RNA <- RenameCells(A84B_RNA,new.names = paste('A84B',colnames(A84B_RNA),sep = '_'))
A84C_RNA <- RenameCells(A84C_RNA,new.names = paste('A84C',colnames(A84C_RNA),sep = '_'))

table(rownames(A68A_RNA) == rownames(A68B_RNA))
table(rownames(A68A_RNA) == rownames(A84B_RNA))
table(rownames(A68A_RNA) == rownames(A84C_RNA))

macaque_RNA <- cbind(A68A_RNA@assays$RNA@counts,A68B_RNA@assays$RNA@counts,A84B_RNA@assays$RNA@counts,A84C_RNA@assays$RNA@counts)

macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA,project = 'macaque',min.cells = 0,min.features = 0)
macaque_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_seurat,features = MT_gene_list)
remove(list = c('A68A_RNA','A68B_RNA','A84B_RNA','A84C_RNA','macaque_RNA'))

#normalize data
macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)

#sctransform for dim reduction
macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

#PCA
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

#find clusters
macaque_RNA_seurat <- FindNeighbors(macaque_RNA_seurat,dims = 1:32)
macaque_RNA_seurat <- FindClusters(macaque_RNA_seurat,resolution = 2)

#UMAP
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:32,reduction = 'pca')
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'seurat_clusters',label = TRUE,repel = TRUE)

#add annotation
macaque_RNA_seurat$species <- 'macaque'
temp <- rownames(macaque_RNA_seurat@meta.data)
temp <- base::lapply(strsplit(x = temp,split = '_'),function(x){
  return(x[1])
})
temp <- unlist(temp)
temp <- as.character(temp)

macaque_RNA_seurat$donor <- temp

DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',split.by = 'donor',label = TRUE,repel = TRUE)

#using annotation from liuyt to see the cluster representation
liuyt_meta <- readRDS(file = './processed_data/macaque_scRNA_seq_annotation_from_liuyt_210612.rds')
rownames(liuyt_meta) <- paste('A',rownames(liuyt_meta),'-1',sep = '')
table(colnames(macaque_RNA_seurat) %in% rownames(liuyt_meta))

cell_list <- dplyr::intersect(rownames(liuyt_meta),colnames(macaque_RNA_seurat))
cluster_representation <- data.frame(cell_id=cell_list,
                                     sunym_cluster=as.character(macaque_RNA_seurat@meta.data[cell_list,"seurat_clusters"]),
                                     liuyt_cluster=as.character(liuyt_meta[cell_list,"SCT_snn_res.1.8"]))

temp <- My_confusion_heatmap(ori = cluster_representation$liuyt_cluster,prd = cluster_representation$sunym_cluster,return_data = TRUE)
temp_matrix <- matrix(nrow = length(unique(temp$prd)),ncol = length(unique(temp$ori)))
colnames(temp_matrix) <- unique(temp$ori)
rownames(temp_matrix) <- unique(temp$prd)

for (i in 1:dim(temp_matrix)[1]) {
  for (j in 1:dim(temp_matrix)[2]) {
    temp_matrix[i,j] <- as.numeric(temp[temp$ori == colnames(temp_matrix)[j] & temp$prd == rownames(temp_matrix)[i],"Prob"])
  }
}

temp_matrix <- temp_matrix[names(sort(base::apply(temp_matrix,1,max),decreasing = TRUE)),
                           names(sort(base::apply(temp_matrix,2,max),decreasing = TRUE))]

pdf(file = './res/step_1_fig_210612/cluster_representation.pdf',width = 10,height = 10)
Heatmap(temp_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_row_names = TRUE,show_column_names = TRUE,
        column_title = 'liuyt cluster',row_title = 'sunym cluster',
        height = 4.9,width = 5.2,
        heatmap_width = unit(8,'inches'),heatmap_height = unit(8,'inches'),
        name = 'Prop')
dev.off()

# find doublet
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
#A68A
A68A_doublet <- read.csv(file = './res/step_1_fig_210612/A68A_scrublet_210613.txt',as.is = TRUE)
A68A_doublet$barcodes <- sub(pattern = '.',fixed = TRUE,replacement = '-',A68A_doublet$barcodes)
A68A_doublet$barcodes <- paste('A68A',A68A_doublet$barcodes,sep = '_')
rownames(A68A_doublet) <- A68A_doublet$barcodes
A68A_doublet <- A68A_doublet[rownames(A68A_doublet) %in% colnames(macaque_RNA_seurat),]
#A68B
A68B_doublet <- read.csv(file = './res/step_1_fig_210612/A68B_scrublet_210613.txt',as.is = TRUE)
A68B_doublet$barcodes <- sub(pattern = '.',fixed = TRUE,replacement = '-',A68B_doublet$barcodes)
A68B_doublet$barcodes <- paste('A68B',A68B_doublet$barcodes,sep = '_')
rownames(A68B_doublet) <- A68B_doublet$barcodes
A68B_doublet <- A68B_doublet[rownames(A68B_doublet) %in% colnames(macaque_RNA_seurat),]
#A84B
A84B_doublet <- read.csv(file = './res/step_1_fig_210612/A84B_scrublet_210613.txt',as.is = TRUE)
A84B_doublet$barcodes <- sub(pattern = '.',fixed = TRUE,replacement = '-',A84B_doublet$barcodes)
A84B_doublet$barcodes <- paste('A84B',A84B_doublet$barcodes,sep = '_')
rownames(A84B_doublet) <- A84B_doublet$barcodes
A84B_doublet <- A84B_doublet[rownames(A84B_doublet) %in% colnames(macaque_RNA_seurat),]
#A84C
A84C_doublet <- read.csv(file = './res/step_1_fig_210612/A84C_scrublet_210613.txt',as.is = TRUE)
A84C_doublet$barcodes <- sub(pattern = '.',fixed = TRUE,replacement = '-',A84C_doublet$barcodes)
A84C_doublet$barcodes <- paste('A84C',A84C_doublet$barcodes,sep = '_')
rownames(A84C_doublet) <- A84C_doublet$barcodes
A84C_doublet <- A84C_doublet[rownames(A84C_doublet) %in% colnames(macaque_RNA_seurat),]

doublet_score <- rbind(A68A_doublet,A68B_doublet,A84B_doublet,A84C_doublet)
doublet_score <- doublet_score[rownames(macaque_RNA_seurat@meta.data),]
macaque_RNA_seurat$doublet_score <- doublet_score$score
macaque_RNA_seurat$doublet <- doublet_score$prediction

remove(list = c('A68A_doublet','A68B_doublet','A84B_doublet','A84C_doublet','doublet_score'))

#remove the doublet 
p1 <- DimPlot(macaque_RNA_seurat,group.by = 'doublet',cols = c('grey','black'))
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p1+p2
sort(table(macaque_RNA_seurat[,macaque_RNA_seurat$doublet == 'True']$seurat_clusters)/table(macaque_RNA_seurat$seurat_clusters),decreasing = TRUE)
#filter cluster 27 and 31
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('27','31'))]
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
#validate by liuyt meta
liuyt_meta <- readRDS(file = './processed_data/macaque_scRNA_seq_annotation_from_liuyt_210612.rds')
rownames(liuyt_meta) <- paste('A',rownames(liuyt_meta),'-1',sep = '')
temp <- macaque_RNA_seurat
cell_list <- dplyr::intersect(colnames(temp),rownames(liuyt_meta))
temp <- temp[,cell_list]
liuyt_meta <- liuyt_meta[cell_list,]
temp$doublet_liuyt <- liuyt_meta$DoubletPrediction
DimPlot(temp,label = TRUE,repel = TRUE,group.by = 'doublet_liuyt',cols = c('black','grey'))
sort(table(temp[,temp$doublet_liuyt == 'Doublet']$seurat_clusters)/table(temp$seurat_clusters),decreasing = TRUE)

#remove donor specific
#note in liuyt results, the donor specific cluster is:21 25 31 35 38 39
Heatmap(temp_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_row_names = TRUE,show_column_names = TRUE,
        column_title = 'liuyt cluster',row_title = 'sunym cluster',
        height = 4.9,width = 5.2,
        heatmap_width = unit(8,'inches'),heatmap_height = unit(8,'inches'),
        name = 'Prop')
#liuyt donor specific cluster refers to cluster in my data:
#21,25,37,23,41,42
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
#de novo
donor_cluster_proportion <- My_Cell_Proportion(macaque_RNA_seurat,split.by = 'seurat_clusters',group.by = 'donor')
ggplot(data = donor_cluster_proportion, mapping = aes(x = seurat_clusters,y = Proportion,fill = donor))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of donor contribution to each cluster',fill = 'donor')+
  scale_fill_manual(values = c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF'))+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()
#may have problem:
#21,23,25,36,37,41,42,43,40
#to be disscussed
#23 only rare in A68A
#40 only rare in A68A
#36 most in A68A
#43 small group, no problem
cluster_36_marker <- FindMarkers(macaque_RNA_seurat,group.by = 'seurat_clusters',
                                 ident.1 = '36',assay = 'RNA',slot = 'data',only.pos = TRUE,
                                 logfc.threshold = 0,min.pct = 0.2)
#keep cluster 23 36 43
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('21','25','37','41','42'))]
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)


###################################################################
##re do the process to finally filter all donor specific clusters##
###################################################################

#reload data
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts

#preprocess
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',min.cells = 0,min.features = 0)

MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(macaque_RNA_seurat))
macaque_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_seurat,features = MT_gene_list)
VlnPlot(macaque_RNA_seurat,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)

macaque_RNA_seurat$species <- 'macaque'
donor_list <- rownames(macaque_RNA_seurat@meta.data)
donor_list <- base::lapply(donor_list,function(x){
  temp <- strsplit(x,split = '_',fixed = TRUE)
  return(temp[[1]][1])
})
donor_list <- unlist(donor_list)
macaque_RNA_seurat$donor <- donor_list
table(macaque_RNA_seurat$donor)

macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

#reduce dimension
use_dims <- 41
macaque_RNA_seurat <- FindNeighbors(macaque_RNA_seurat,dims = 1:use_dims)
macaque_RNA_seurat <- FindClusters(macaque_RNA_seurat,resolution = 2)
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:use_dims,reduction = 'pca')
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'seurat_clusters',split.by = 'donor',label = TRUE,repel = TRUE)

#show donor specific clusters
donor_cluster_proportion <- My_Cell_Proportion(macaque_RNA_seurat,split.by = 'seurat_clusters',group.by = 'donor')
ggplot(data = donor_cluster_proportion, mapping = aes(x = seurat_clusters,y = Proportion,fill = donor))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of donor contribution to each cluster',fill = 'donor')+
  scale_fill_manual(values = c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF'))+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()
# cluster 36
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('36'))]

###################################################################
##re do the process to finally filter all donor specific clusters##
###################################################################
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts

#preprocess
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',min.cells = 0,min.features = 0)

MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(macaque_RNA_seurat))
macaque_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_seurat,features = MT_gene_list)
VlnPlot(macaque_RNA_seurat,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)

macaque_RNA_seurat$species <- 'macaque'
donor_list <- rownames(macaque_RNA_seurat@meta.data)
donor_list <- base::lapply(donor_list,function(x){
  temp <- strsplit(x,split = '_',fixed = TRUE)
  return(temp[[1]][1])
})
donor_list <- unlist(donor_list)
macaque_RNA_seurat$donor <- donor_list
table(macaque_RNA_seurat$donor)

macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

#reduce dimension
use_dims <- 37
macaque_RNA_seurat <- FindNeighbors(macaque_RNA_seurat,dims = 1:use_dims)
macaque_RNA_seurat <- FindClusters(macaque_RNA_seurat,resolution = 2)
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:use_dims,reduction = 'pca')
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'seurat_clusters',split.by = 'donor',label = TRUE,repel = TRUE)

# filter donor specific clusters
donor_cluster_proportion <- My_Cell_Proportion(macaque_RNA_seurat,split.by = 'seurat_clusters',group.by = 'donor')
ggplot(data = donor_cluster_proportion, mapping = aes(x = seurat_clusters,y = Proportion,fill = donor))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6)+
  labs(title = 'Proportion of donor contribution to each cluster',fill = 'donor')+
  scale_fill_manual(values = c('#FED439FF','#709AE1FF','#8A9197FF','#D2AF81FF'))+
  theme_cowplot()+
  theme(axis.text.x = element_text(vjust = 0.5, hjust = 0.5, angle = 45),
        plot.title = element_text(),
        legend.title = element_text(),
        legend.text = element_text(),
        panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
        axis.line = element_blank())+
  xlab('')+
  CenterTitle()

#save express matrix
saveRDS(macaque_RNA_seurat@assays$RNA@counts,file = './processed_data/macaque_RNA_filted_express_matrix_210618.rds')
