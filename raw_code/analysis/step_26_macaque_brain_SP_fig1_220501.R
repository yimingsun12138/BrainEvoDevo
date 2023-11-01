#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque brain SP fig1                                           ##
## Data: 2022.05.01                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# create color param ------------------------------------------------------
# color_param <- data.frame(term = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','Cyc-S','Cyc-G2M','RG-1','RG-2','RG-3','OPC','Ependymal','Mic','End','Per'),
#                           color = c('#916848','#f5b390','#37526d','#0d74b6','#6ebcbc','#fbdf72','#41a30d','#faa818','#367d7d','#725ca5','#d33502','#a6d9ee','#60824f','#bed678','#e0598b','#342739'))
# write.csv(color_param,file = './data/parameter/scRNAseq_color_param.csv')

# load data ---------------------------------------------------------------

#integrated RNA data
macaque_integrated_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
#scRNA data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220504_summary/macaque_RNA_Seurat_220507.rds')
#scATAC data
macaque_ATAC_ArchR <- readRDS(file = './ArchR/processed_data/macaque_ATAC_ArchR_220506/Save-ArchR-Project.rds')

# fig ---------------------------------------------------------------------
## scRNA-seq UMAP ----------------------------------------------------------
# macaque_RNA_Seurat <- macaque_integrated_Seurat[,macaque_integrated_Seurat$tech == 'scRNA']
# meta_data <- macaque_RNA_Seurat@meta.data
# meta_data <- meta_data[,c('batch','donor','tech','cell_type')]
# macaque_RNA_Seurat <- macaque_RNA_Seurat@assays$RNA@counts
# macaque_RNA_Seurat <- CreateSeuratObject(counts = macaque_RNA_Seurat,project = 'macaque',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
# 
# macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 3000,vars.to.regress = c('nCount_RNA','batch','donor'),npcs = 50,preprocess = TRUE)
# macaque_RNA_Seurat <- my_process_seurat(object = macaque_RNA_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 27,resolution = c(0.5,1,1.5),group.by = 'cell_type',label = TRUE)
# macaque_RNA_Seurat$old_cell_type <- macaque_RNA_Seurat$cell_type
# macaque_RNA_Seurat$cell_type <- NA
# 
# #refine annotate
# DimPlot(macaque_RNA_Seurat,group.by = 'RNA_snn_res.1.5',label = TRUE,repel = TRUE) + 
#   DimPlot(macaque_RNA_Seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
# 
# #cluster 19 34 OPC
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('19','34'),"cell_type"] <- 'OPC'
# 
# #cluster 27 RG-3
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('27'),"cell_type"] <- 'RG-3'
# 
# #cluster 17 RG-2
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('17'),"cell_type"] <- 'RG-2'
# 
# #cluster 10 RG-1
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('10'),"cell_type"] <- 'RG-1'
# 
# #cluster 32 Ependymal
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('32'),"cell_type"] <- 'Ependymal'
# 
# #cluster 9 Cycling
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('9'),"cell_type"] <- 'Cycling'
# 
# #cluster 11 IP
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('11'),"cell_type"] <- 'IP'
# 
# #cluster 4 6 8 13 14 18 Ex-1
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('4','6','8','13','14','18'),"cell_type"] <- 'Ex-1'
# 
# #cluster 2 3 12 26
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('2','3','12','26'),"cell_type"] <- 'Ex-2'
# 
# #cluster 15 20 23 33 Ex-3
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('15','20','23','33'),"cell_type"] <- 'Ex-3'
# 
# #cluster 0 5 21 InCGE
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('0','5','21'),"cell_type"] <- 'InCGE'
# 
# #cluster 1 7 16 22 25 29
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('1','7','16','22','25','29'),"cell_type"] <- 'InMGE'
# 
# #cluster 30 31 Per
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('30','31'),"cell_type"] <- 'Per'
# 
# #cluster 24 End
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('24'),"cell_type"] <- 'End'
# 
# #cluster 28 Mic
# macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$RNA_snn_res.1.5 %in% c('28'),"cell_type"] <- 'Mic'
# 
# DimPlot(macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# 
# #sub cluster Cycling
# temp <- macaque_RNA_Seurat
# Idents(temp) <- 'cell_type'
# temp <- FindSubCluster(object = temp,cluster = 'Cycling',graph.name = 'RNA_snn',resolution = 0.1)
# DimPlot(temp,group.by = 'sub.cluster',label = TRUE,repel = TRUE)
# 
# macaque_RNA_Seurat@meta.data[temp$sub.cluster %in% c('Cycling_1'),"cell_type"] <- 'Cyc-G2M'
# macaque_RNA_Seurat@meta.data[temp$sub.cluster %in% c('Cycling_0'),"cell_type"] <- 'Cyc-S'
# 
# #the consistence with old cell type
# scibet::Confusion_heatmap(ori = macaque_RNA_Seurat$old_cell_type,prd = macaque_RNA_Seurat$cell_type)
# macaque_RNA_Seurat@meta.data <- macaque_RNA_Seurat@meta.data[,!(colnames(macaque_RNA_Seurat@meta.data) %in% c('old_cell_type'))]
# 
# #save data
# saveRDS(macaque_RNA_Seurat,file = './processed_data/220504_summary/macaque_RNA_Seurat_220504.rds')

# #re-do umap
# macaque_RNA_Seurat <- RunUMAP(object = macaque_RNA_Seurat,dims = 1:27,reduction = 'pca',n.neighbors = 50,metric = 'cosine',min.dist = 0.6)
# #save data
# saveRDS(macaque_RNA_Seurat,file = './processed_data/220504_summary/macaque_RNA_Seurat_220507.rds')

#color param
color_param <- read.csv(file = './data/parameter/scRNAseq_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

#plot
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_UMAP.pdf',width = 7,height = 7)
DimPlot(macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'umap') +
  theme_cowplot() +
  scale_color_manual(values = col_value) +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'snRNA-seq')
dev.off()

## scATAC-seq UMAP ---------------------------------------------------------

# #color param
# color_param <- data.frame(term = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG','OPC','Mic','End/Per'),
#                           color = c('#916848','#f5b390','#37526d','#0d74b6','#6ebcbc','#fbdf72','#f47d2b','#a6d9ee','#bed678','#7e1416'))
# write.csv(color_param,file = './data/parameter/scATACseq_color_param.csv')

#plot
color_param <- read.csv(file = './data/parameter/scATACseq_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_UMAP.pdf',width = 7,height = 7)
my_dimplot(embedding = macaque_ATAC_ArchR@embeddings$harmonyUMAP$df,
           meta_data = as.data.frame(macaque_ATAC_ArchR@cellColData),
           group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_cowplot() + 
  scale_color_manual(values = col_value) + 
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'scATAC-seq')
dev.off()

#plot label transfer based annotation UMAP
color_param <- read.csv(file = './data/parameter/scATACseq_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

p1 <- my_dimplot(embedding = macaque_ATAC_ArchR@embeddings$harmonyUMAP$df,
                 meta_data = as.data.frame(macaque_ATAC_ArchR@cellColData),
                 group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_cowplot() + 
  scale_color_manual(values = col_value) + 
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'scATAC-seq annotation')

color_param <- read.csv(file = './data/parameter/scRNAseq_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

p2 <- my_dimplot(embedding = macaque_ATAC_ArchR@embeddings$harmonyUMAP$df,
                 meta_data = as.data.frame(macaque_ATAC_ArchR@cellColData),
                 group.by = 'predictedGroup',label = TRUE,repel = TRUE) + 
  theme_cowplot() + 
  scale_color_manual(values = col_value) + 
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'scATAC-seq label transfer based annotation')

pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_label_transfer_based_annotation.pdf',width = 14,height = 7)
p2+p1+plot_layout(ncol = 2)
dev.off()

#cluster based anotation compare to label transfer based annotation
temp <- function(ori,prd){
  cross.validation.filt <- tibble(ori = ori, prd = prd) %>% 
    dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
  cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
  return(cross.validation.filt)
}

temp <- temp(ori = macaque_ATAC_ArchR$cell_type,prd = macaque_ATAC_ArchR$predictedGroup)
temp$ori <- factor(temp$ori,levels = c('End/Per','Mic','OPC','RG','IP','Ex-1','Ex-2','Ex-3','InCGE','InMGE'))
temp$prd <- factor(temp$prd,levels = c('End','Per','Mic','OPC','RG-3','RG-2','RG-1','IP','Ex-1','Ex-2','Ex-3','InCGE','InMGE'))

pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_label_transfer_clustering_annotation_confusion_heatmap.pdf',width = 6,height = 6)
ggplot(data = temp,aes(ori, prd, fill = Prob)) + 
  geom_tile() + 
  theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), 
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_viridis() + 
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 14,face = 'bold',colour = 'black')) + 
  ylab('label transfer based annotation') + xlab('clustering based annotation')
dev.off()

## vlnplot for snRNA-seq quality -------------------------------------------
#create color param
# color_param <- data.frame(term = c('A68A','A68B','A84B','A84C','A50A','A50B','A82A','A82B'),
#                           color = c('#FFB300','#803E75','#FF6800','#A6BDD7','#00538A','#CEA262','#FF8E00','#B32851'))
# write.csv(color_param,file = './data/parameter/donor_color_param.csv')

#MT ratio
MT_gene <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
macaque_RNA_Seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_Seurat,features = MT_gene,assay = 'RNA')

#plot
color_param<- read.csv(file = './data/parameter/donor_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

meta_data <- macaque_RNA_Seurat@meta.data
meta_data$donor <- factor(meta_data$donor,levels = c('A68A','A68B','A84B','A84C','A50A','A82B'))

p1 <- ggplot(meta_data,aes(x = donor,y = nCount_RNA,fill = donor,color = donor)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + ylim(c(0,15000)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1)) + 
  ylab('Reads')

p2 <- ggplot(meta_data,aes(x = donor,y = nFeature_RNA,fill = donor,color = donor)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1)) + 
  ylab('Genes')

p3 <- ggplot(meta_data,aes(x = donor,y = percent.mt,fill = donor,color = donor)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1),
        axis.ticks = element_blank()) + 
  ylab('MT reads %')

pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_QC_vlnplot.pdf',width = 4,height = 6)
p1+p2+p3+plot_layout(ncol = 1)
dev.off()

## vlnplot for scATAC-seq quality ------------------------------------------
meta_data <- as.data.frame(macaque_ATAC_ArchR@cellColData)
meta_data$Sample <- factor(meta_data$Sample,levels = c('A68A','A68B','A84B','A84C','A50A','A82A','A82B'))
meta_data$lognFrag <- log(meta_data$nFrags)

#plot
color_param<- read.csv(file = './data/parameter/donor_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

p1 <- ggplot(meta_data,aes(x = Sample,y = lognFrag,fill = Sample,color = Sample)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1)) + 
  ylab('log(nFrag)')

p2 <- ggplot(meta_data,aes(x = Sample,y = TSSEnrichment,fill = Sample,color = Sample)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + ylim(c(0,35)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1)) + 
  ylab('TSS enrichment')

p3 <- ggplot(meta_data,aes(x = Sample,y = FRIP,fill = Sample,color = Sample)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value) + 
  scale_color_manual(values = col_value) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1),
        axis.ticks = element_blank()) + 
  ylab('FRIP')

pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_QC_vlnplot.pdf',width = 4,height = 6)
p1+p2+p3+plot_layout(ncol = 1)
dev.off()

## correlation of all snRNA-seq sample -------------------------------------
#create pseudobulk matrix
macaque_RNA_Seurat@meta.data[1:3,]
temp <- Seurat::AggregateExpression(object = macaque_RNA_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'donor',slot = 'counts',verbose = TRUE)
temp <- temp$RNA
temp <- temp[rowSums(temp) > 1,]

#create correlation matrix
cor_matrix <- matrix(data = NA,nrow = dim(temp)[2],ncol = dim(temp)[2])
cor_matrix <- as.data.frame(cor_matrix)
rownames(cor_matrix) <- colnames(temp)
colnames(cor_matrix) <- colnames(temp)
for (i in rownames(cor_matrix)) {
  for (j in colnames(cor_matrix)) {
    cor_matrix[i,j] <- cor(temp[,i],temp[,j])
  }
}

#plot
group_list <- factor(c('batch_2','batch_1','batch_1','batch_2','batch_1','batch_1'),levels = c('batch_1','batch_2'))
annotation_top <- HeatmapAnnotation(batch = c('batch_2','batch_1','batch_1','batch_2','batch_1','batch_1'),
                                    col = list(batch = c('batch_1' = '#34a047','batch_2' = '#f47e1f')),
                                    which = 'column')
annotation_row <- HeatmapAnnotation(batch = c('batch_2','batch_1','batch_1','batch_2','batch_1','batch_1'),
                                    col = list(batch = c('batch_1' = '#34a047','batch_2' = '#f47e1f')),
                                    which = 'row')
col_fun <- colorRamp2(breaks = c(0.8,0.9,1),colors = c('#3361A5','#C1D5DC','#A31D1D'))

pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_donor_correlation.pdf',width = 6,height = 4.5)
Heatmap(matrix = cor_matrix,show_row_names = TRUE,show_column_names = TRUE,
        cluster_columns = cluster_within_group(mat = cor_matrix,factor = group_list),
        cluster_rows = cluster_within_group(mat = t(cor_matrix),factor = group_list),
        top_annotation = annotation_top,left_annotation = annotation_row,col = col_fun,
        height = unit(2.5,'inches'),width = unit(2.5,'inches'),name = 'correlation')
dev.off()

## correlation of all scATAC-seq sample ------------------------------------
#get tile matrix
getAvailableMatrices(ArchRProj = macaque_ATAC_ArchR)
temp <- getMatrixFromProject(ArchRProj = macaque_ATAC_ArchR,useMatrix = 'TileMatrix',verbose = TRUE,binarize = TRUE)
gene_list <- paste(temp@elementMetadata$seqnames,temp@elementMetadata$start,sep = '_')
temp <- temp@assays@data$TileMatrix
rownames(temp) <- gene_list

#create pseudo bulk
temp <- CreateSeuratObject(counts = temp,project = 'temp',assay = 'RNA',
                           meta.data = as.data.frame(macaque_ATAC_ArchR@cellColData),
                           min.cells = 0,min.features = 0)
temp <- Seurat::AggregateExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'donor',slot = 'counts',verbose = TRUE)
temp <- temp$RNA
temp <- temp[rowSums(temp) > 1,]

#create correlation matrix
cor_matrix <- matrix(data = NA,nrow = dim(temp)[2],ncol = dim(temp)[2])
cor_matrix <- as.data.frame(cor_matrix)
rownames(cor_matrix) <- colnames(temp)
colnames(cor_matrix) <- colnames(temp)
for (i in rownames(cor_matrix)) {
  for (j in colnames(cor_matrix)) {
    cor_matrix[i,j] <- cor(temp[,i],temp[,j])
  }
}

#plot
group_list <- factor(c('batch_2','batch_1','batch_1','batch_2','batch_2','batch_1','batch_1'),levels = c('batch_1','batch_2'))
annotation_top <- HeatmapAnnotation(batch = c('batch_2','batch_1','batch_1','batch_2','batch_2','batch_1','batch_1'),
                                    col = list(batch = c('batch_1' = '#34a047','batch_2' = '#f47e1f')),
                                    which = 'column')
annotation_row <- HeatmapAnnotation(batch = c('batch_2','batch_1','batch_1','batch_2','batch_2','batch_1','batch_1'),
                                    col = list(batch = c('batch_1' = '#34a047','batch_2' = '#f47e1f')),
                                    which = 'row')
col_fun <- colorRamp2(breaks = c(0.8,0.9,1),colors = c('#3361A5','#C1D5DC','#A31D1D'))

pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_donor_correlation.pdf',width = 6,height = 5)
Heatmap(matrix = cor_matrix,show_row_names = TRUE,show_column_names = TRUE,
        cluster_columns = cluster_within_group(mat = cor_matrix,factor = group_list),
        cluster_rows = cluster_within_group(mat = t(cor_matrix),factor = group_list),
        top_annotation = annotation_top,left_annotation = annotation_row,col = col_fun,
        height = unit(2.92,'inches'),width = unit(2.92,'inches'),name = 'correlation')
dev.off()

## snRNA-seq data UMAP by donor --------------------------------------------
#color param
color_param <- read.csv(file = './data/parameter/donor_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

#plot
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_UMAP_by_donor.pdf',width = 7,height = 7)
DimPlot(macaque_RNA_Seurat,group.by = 'donor',label = FALSE,reduction = 'umap') +
  theme_cowplot() +
  scale_color_manual(values = col_value) +
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'snRNA-seq UMAP by donor')
dev.off()

## scATAC-seq data UMAP by donor -------------------------------------------

#color parameter
color_param <- read.csv(file = './data/parameter/donor_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

#plot
pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_UMAP_by_donor.pdf',width = 7,height = 7)
my_dimplot(embedding = macaque_ATAC_ArchR@embeddings$harmonyUMAP$df,
           meta_data = as.data.frame(macaque_ATAC_ArchR@cellColData),
           group.by = 'donor',label = FALSE) + 
  theme_cowplot() + 
  scale_color_manual(values = col_value) + 
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'scATAC-seq UMAP by donor')
dev.off()

## scATAC-seq fragment length distribution ---------------------------------
p <- plotFragmentSizes(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'donor', 
                       returnDF = FALSE)

#color parameter
color_param <- read.csv(file = './data/parameter/donor_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

#plot
pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_fragment_size_distribution_by_donor.pdf',width = 4,height = 4)
p + theme_cowplot() + scale_color_manual(values = col_value[unique(macaque_ATAC_ArchR$donor)]) + 
  theme(aspect.ratio = 1) + 
  xlab('Fragment size') + ylab('Density')
dev.off()

## scATAC-seq TSS enrichment plot ------------------------------------------
p <- plotTSSEnrichment(ArchRProj = macaque_ATAC_ArchR,
                       groupBy = 'donor',
                       returnDF = FALSE)

#color parameter
color_param <- read.csv(file = './data/parameter/donor_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

#plot
pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_TSS_enrichment_by_donor.pdf',width = 4,height = 4)
p + theme_cowplot() + scale_color_manual(values = col_value[unique(macaque_ATAC_ArchR$donor)]) + 
  theme(aspect.ratio = 1) + 
  xlab('Position from TSS (bp)') + ylab('Normalized insertion profile')
dev.off()

## scATAC-seq nFrag TSS dotplot --------------------------------------------
df <- macaque_ATAC_ArchR@cellColData[,c('nFrags','TSSEnrichment','Sample')]
df$nFrags <- log10(df$nFrags)
colnames(df) <- replace(x = colnames(df),list = c(colnames(df) == 'nFrags'),values = 'log10(nFrags)')
df

for (i in unique(df$Sample)) {
  p <- df[df$Sample == i,]
  assign(x = paste('QC',as.character(i),sep = '_'),
         value = ggPoint(
           x = p[,1], 
           y = p[,2], 
           colorDensity = TRUE,
           continuousSet = "blueYellow",
           pal = viridisLite::viridis(begin = 0,end = 1,option = 'D',n = 9),
           xlabel = "Log10 Unique Fragments",
           ylabel = "TSS Enrichment",
           title = as.character(i),
           xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
           ylim = c(0, quantile(df[,2], probs = 0.99))) + 
           theme(plot.title = element_text(hjust = 0.5)) + 
           geom_hline(yintercept = 4, lty = "dashed") + 
           geom_vline(xintercept = 3, lty = "dashed") + 
           theme(legend.position = 'none'))
}

pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_nFrag_vs_TSS_dotplot_no_legend.pdf',width = 12,height = 6)
QC_A68A+QC_A68B+QC_A84B+QC_A84C+QC_A50A+QC_A82A+QC_A82B+plot_layout(ncol = 4)
dev.off()

for (i in unique(df$Sample)) {
  p <- df[df$Sample == i,]
  assign(x = paste('QC',as.character(i),sep = '_'),
         value = ggPoint(
           x = p[,1], 
           y = p[,2], 
           colorDensity = TRUE,
           continuousSet = "blueYellow",
           pal = viridisLite::viridis(begin = 0,end = 1,option = 'D',n = 9),
           xlabel = "Log10 Unique Fragments",
           ylabel = "TSS Enrichment",
           title = as.character(i),
           xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
           ylim = c(0, quantile(df[,2], probs = 0.99))) + 
           theme(plot.title = element_text(hjust = 0.5)) + 
           geom_hline(yintercept = 4, lty = "dashed") + 
           geom_vline(xintercept = 3, lty = "dashed"))
}

pdf(file = './res/step_26_fig_220501/macaque_ATAC_ArchR_nFrag_vs_TSS_dotplot_with_legend.pdf',width = 12,height = 8)
QC_A68A+QC_A68B+QC_A84B+QC_A84C+QC_A50A+QC_A82A+QC_A82B+plot_layout(ncol = 4)
dev.off()

## classic cell type marker ---------------------------------------------------

#RG
marker_list <- c('ATP1A2','PAX6','SLC1A3')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_RG_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q95',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#vRG
marker_list <- c('PALLD','CRYAB')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_vRG_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q95',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#oRG
marker_list <- c('MOXD1','HOPX','TNC')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_oRG_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q95',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#Cycling
marker_list <- c('TOP2A','MKI67','CLSPN','AURKA')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_Cycling_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q95',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#RG-1
marker_list <- c('NOTCH2','TMEM132D','IQGAP2')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_RG_1_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q99',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#RG-2
marker_list <- c('ATP13A4','SLC4A4','AQP4','MMD2')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_RG_2_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q99',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#RG-3
marker_list <- c('EGFR','SPATA6','CHST9')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_RG_3_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q99',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#OPC
marker_list <- c('SOX10','OLIG2')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_OPC_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q99',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()

#Mic
marker_list <- c('CD84','CX3CR1')
pdf(file = './res/step_26_fig_220501/macaque_RNA_Seurat_Mic_marker_featureplot.pdf',width = 4*length(marker_list),height = 4)
FeaturePlot(object = macaque_RNA_Seurat,features = marker_list,slot = 'data',
            cols = ArchRPalettes$greyMagma,max.cutoff = 'q99',ncol = length(marker_list)) + 
  theme(aspect.ratio = 1)
dev.off()