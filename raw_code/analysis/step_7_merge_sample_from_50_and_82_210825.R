#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: merge sample from 50 and 82                                     ##
## Data: 2021.08.26                                                                ##
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
library(patchwork)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

# load data ---------------------------------------------------------------
A50A_RNA <- Read10X(data.dir = './data/scRNA-seq/210822/A-50A/outs/filtered_feature_bc_matrix/')
A82A_RNA <- Read10X(data.dir = './data/scRNA-seq/210822/A-82A/outs/filtered_feature_bc_matrix/')
A82B_RNA <- Read10X(data.dir = './data/scRNA-seq/210822/A-82B/outs/filtered_feature_bc_matrix/')

# create seurat object
A50A_RNA <- CreateSeuratObject(counts = A50A_RNA,project = 'A50A',min.cells = 0,min.features = 200)
A82A_RNA <- CreateSeuratObject(counts = A82A_RNA,project = 'A82A',min.cells = 0,min.features = 200)
A82B_RNA <- CreateSeuratObject(counts = A82B_RNA,project = 'A82B',min.cells = 0,min.features = 200)

# QC ----------------------------------------------------------------------
MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(A50A_RNA))

A50A_RNA[['percent.mt']] <- PercentageFeatureSet(object = A50A_RNA,features = MT_gene_list)
A82A_RNA[['percent.mt']] <- PercentageFeatureSet(object = A82A_RNA,features = MT_gene_list)
A82B_RNA[['percent.mt']] <- PercentageFeatureSet(object = A82B_RNA,features = MT_gene_list)

VlnPlot(A50A_RNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(A50A_RNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A50A_RNA <- subset(A50A_RNA,subset = nFeature_RNA < 4000 & nCount_RNA < 10000)

VlnPlot(A82A_RNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(A82A_RNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A82A_RNA <- subset(A82A_RNA,subset = nFeature_RNA < 4000 & nCount_RNA < 10000)

VlnPlot(A82B_RNA,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),ncol = 3)
FeatureScatter(A82B_RNA,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
A82B_RNA <- subset(A82B_RNA,subset = nFeature_RNA < 5000 & nCount_RNA < 15000)

# merge -------------------------------------------------------------------
A50A_RNA <- RenameCells(object = A50A_RNA,new.names = paste('A50A',colnames(A50A_RNA),sep = '_'))
A82A_RNA <- RenameCells(object = A82A_RNA,new.names = paste('A82A',colnames(A82A_RNA),sep = '_'))
A82B_RNA <- RenameCells(object = A82B_RNA,new.names = paste('A82B',colnames(A82B_RNA),sep = '_'))

table(rownames(A50A_RNA) == rownames(A82A_RNA))
table(rownames(A50A_RNA) == rownames(A82B_RNA))

macaque_RNA <- cbind(A50A_RNA@assays$RNA@counts,A82A_RNA@assays$RNA@counts,A82B_RNA@assays$RNA@counts)

macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA,project = 'macaque',min.cells = 0,min.features = 0)
macaque_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_seurat,features = MT_gene_list)
remove(list = c('A50A_RNA','A82A_RNA','A82B_RNA','macaque_RNA'))

# process seurat object ---------------------------------------------------
macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)

macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

# filter clusters ---------------------------------------------------------
macaque_RNA_seurat <- FindNeighbors(macaque_RNA_seurat,dims = 1:38)
macaque_RNA_seurat <- FindClusters(macaque_RNA_seurat,resolution = 1)
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:38,reduction = 'pca')
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

#filter  doublet
A50A_doublet <- read.csv(file = './res/step_7_fig_210825/A50A_scrublet_210825.csv')
A50A_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,A50A_doublet$barcodes)
A50A_doublet$barcodes <- paste('A50A',A50A_doublet$barcodes,sep = '_')
table(A50A_doublet$barcodes %in% colnames(macaque_RNA_seurat))
table(duplicated(A50A_doublet$barcodes))
rownames(A50A_doublet) <- A50A_doublet$barcodes
A50A_doublet <- A50A_doublet[A50A_doublet$barcodes %in% colnames(macaque_RNA_seurat),]

A82A_doublet <- read.csv(file = './res/step_7_fig_210825/A82A_scrublet_210825.csv')
A82A_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,A82A_doublet$barcodes)
A82A_doublet$barcodes <- paste('A82A',A82A_doublet$barcodes,sep = '_')
table(A82A_doublet$barcodes %in% colnames(macaque_RNA_seurat))
table(duplicated(A82A_doublet$barcodes))
rownames(A82A_doublet) <- A82A_doublet$barcodes
A82A_doublet <- A82A_doublet[A82A_doublet$barcodes %in% colnames(macaque_RNA_seurat),]

A82B_doublet <- read.csv(file = './res/step_7_fig_210825/A82B_scrublet_210825.csv')
A82B_doublet$barcodes <- sub(pattern = '.',replacement = '-',fixed = TRUE,A82B_doublet$barcodes)
A82B_doublet$barcodes <- paste('A82B',A82B_doublet$barcodes,sep = '_')
table(A82B_doublet$barcodes %in% colnames(macaque_RNA_seurat))
table(duplicated(A82B_doublet$barcodes))
rownames(A82B_doublet) <- A82B_doublet$barcodes
A82B_doublet <- A82B_doublet[A82B_doublet$barcodes %in% colnames(macaque_RNA_seurat),]

doublet_list <- rbind(A50A_doublet,A82A_doublet,A82B_doublet)
doublet_list <- doublet_list[rownames(macaque_RNA_seurat@meta.data),]
macaque_RNA_seurat$doublet_score <- doublet_list$score
macaque_RNA_seurat$doublet <- doublet_list$prediction

p1 <- DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'doublet',label = FALSE,cols = c('lightgrey','red'))

p1+p2


temp <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'doublet',split.by = 'seurat_clusters')
temp_col <- c('lightgrey','black')
names(temp_col) <- c('False','True')

ggplot(temp,aes(x=seurat_clusters,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = temp_col) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

# filter donor specific
temp <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(temp,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)
# filter cluster 0,22,27,28,29
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('0','22','27','28','29'))]


# filter low quality cluster
FeaturePlot(macaque_RNA_seurat,features = c('nFeature_RNA','nCount_RNA','percent.mt'))
VlnPlot(macaque_RNA_seurat,features = c('nFeature_RNA','nCount_RNA','percent.mt'))

# redo filter -------------------------------------------------------------
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',min.cells = 0,min.features = 200)

macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)

macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

macaque_RNA_seurat <- FindNeighbors(macaque_RNA_seurat,dims = 1:44)
macaque_RNA_seurat <- FindClusters(macaque_RNA_seurat,resolution = 1.5)
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:44,reduction = 'pca')
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

doublet_list <- doublet_list[colnames(macaque_RNA_seurat),]
macaque_RNA_seurat$doublet <- doublet_list$prediction

temp <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'doublet',split.by = 'seurat_clusters')
ggplot(temp,aes(x=seurat_clusters,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = temp_col) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)


temp <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'donor',split.by = 'seurat_clusters')
ggplot(temp,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)

# filter cluster 21
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$seurat_clusters %in% c('21'))]

# final process -----------------------------------------------------------
## pre_process -------------------------------------------------------------
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',min.cells = 0,min.features = 200)
macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)


## cluster and umap --------------------------------------------------------
macaque_RNA_seurat <- FindNeighbors(macaque_RNA_seurat,dims = 1:45)
macaque_RNA_seurat <- FindClusters(macaque_RNA_seurat,resolution = 1.5)
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:45,reduction = 'pca')
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'seurat_clusters',label = TRUE,repel = TRUE)

macaque_RNA_seurat$clusters <- macaque_RNA_seurat$seurat_clusters

## QC --------------------------------------------------------------
macaque_RNA_seurat$species <- 'macaque'
temp <- rownames(macaque_RNA_seurat@meta.data)
temp <- base::lapply(strsplit(x = temp,split = '_'),function(x){
  return(x[1])
})
temp <- unlist(temp)
temp <- as.character(temp)

macaque_RNA_seurat$donor <- temp

DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',split.by = 'donor',label = TRUE,repel = TRUE)

doublet_list <- doublet_list[rownames(macaque_RNA_seurat@meta.data),]
macaque_RNA_seurat$doublet_score <- doublet_list$score
macaque_RNA_seurat$doublet <- doublet_list$prediction

pdf(file = './res/step_7_fig_210825/macaque_RNA_seurat_dimplot_split_by_donor.pdf',width = 21,height = 7)
DimPlot(macaque_RNA_seurat,group.by = 'clusters',label = TRUE,repel = TRUE,split.by = 'donor') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = '')
dev.off()

p1 <- DimPlot(macaque_RNA_seurat,group.by = 'clusters',label = TRUE,repel = TRUE) + theme(legend.position = 'none')
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'doublet',label = FALSE,cols = c('lightgrey','red'))
pdf(file = './res/step_7_fig_210825/macaque_RNA_seurat_doublet_dimplot.pdf',width = 14,height = 7)
p1+p2
dev.off()

temp <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'donor',split.by = 'seurat_clusters')
pdf(file = './res/step_7_fig_210825/macaque_RNA_seurat_donor_proportion_by_cluster.pdf',width = 12,height = 6)
ggplot(temp,aes(x=seurat_clusters,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6)
dev.off()

DefaultAssay(macaque_RNA_seurat) <- 'RNA'
MT_gene_list <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene_list <- as.character(MT_gene_list)
table(MT_gene_list %in% rownames(macaque_RNA_seurat))
macaque_RNA_seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_seurat,features = MT_gene_list)

p1 <- VlnPlot(macaque_RNA_seurat,features = 'nCount_RNA',group.by = 'donor',assay = 'RNA',pt.size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6,
        plot.title = element_blank(),
        legend.position = 'none') + 
  ylab('nCount_RNA')

p2 <- VlnPlot(macaque_RNA_seurat,features = 'nFeature_RNA',group.by = 'donor',assay = 'RNA',pt.size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6,
        plot.title = element_blank(),
        legend.position = 'none') + 
  ylab('nFeature_RNA')

p3 <- VlnPlot(macaque_RNA_seurat,features = 'percent.mt',group.by = 'donor',assay = 'RNA',pt.size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6,
        plot.title = element_blank(),
        legend.position = 'none') + 
  ylab('percent.mt')

pdf(file = './res/step_7_fig_210825/QC_plot.pdf',width = 4,height = 6)
p1+p2+p3+plot_layout(nrow = 3,heights = 0.5)
dev.off()

p1 <- VlnPlot(macaque_RNA_seurat,features = 'nCount_RNA',group.by = 'clusters',assay = 'RNA',pt.size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 0.3,
        plot.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank()) + 
  ylab('nCount_RNA')

p2 <- VlnPlot(macaque_RNA_seurat,features = 'nFeature_RNA',group.by = 'clusters',assay = 'RNA',pt.size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 0.3,
        plot.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank()) + 
  ylab('nFeature_RNA')

p3 <- VlnPlot(macaque_RNA_seurat,features = 'percent.mt',group.by = 'clusters',assay = 'RNA',pt.size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 0.3,
        plot.title = element_blank(),
        legend.position = 'none',
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.x = element_blank()) + 
  ylab('percent.mt')

pdf(file = './res/step_7_fig_210825/QC_plot_by_clusters.pdf',width = 12,height = 10)
p1+p2+p3+plot_layout(nrow = 3,heights = 0.5)
dev.off()

p1 <- FeaturePlot(object = macaque_RNA_seurat,features = c('nCount_RNA')) + 
  theme_classic() + 
  theme(aspect.ratio = 1)

p2 <- FeaturePlot(object = macaque_RNA_seurat,features = c('nFeature_RNA')) + 
  theme_classic() + 
  theme(aspect.ratio = 1)

p3 <- FeaturePlot(object = macaque_RNA_seurat,features = c('percent.mt')) + 
  theme_classic() + 
  theme(aspect.ratio = 1)

p4 <- DimPlot(object = macaque_RNA_seurat,group.by = 'clusters',label = TRUE,repel = TRUE) + 
  theme_classic() + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_7_fig_210825/QC_plot_featureplot_form.pdf',width = 12,height = 10)
p4+p1+p2+p3+plot_layout(ncol = 2)
dev.off()


## marker and annotation ---------------------------------------------------
#round 1 annotation
dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Ast=c('AQP4'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1'),
                                             OPC=c('SOX10','OLIG2'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7'),
                                             In=c('DLX5','GAD2','SP8','SST','MEIS2')),
                             group.by = 'clusters', cols = c('#2CA02CFF','white','#D62728FF'),
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

#cluster 21 End
#cluster 27 Per
#cluster 25 Mic
#cluster 14 Cycling
#cluster 32 Ip
#cluster 17 OPC
macaque_RNA_seurat$cell_type <- macaque_RNA_seurat$clusters
macaque_RNA_seurat@meta.data[,"cell_type"] <- as.character(macaque_RNA_seurat@meta.data[,"cell_type"])

macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters == '21',"cell_type"] <- 'End'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters == '27',"cell_type"] <- 'Per'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters == '25',"cell_type"] <- 'Mic'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters == '14',"cell_type"] <- 'Cycling'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters == '32',"cell_type"] <- 'IP'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters == '17',"cell_type"] <- 'OPC'

DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)


#round 2 annotation
dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Ast=c('AQP4'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7'),
                                             In=c('DLX5','GAD2','SP8','SST','MEIS2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
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

#cluster 24 oIP
#cluster 13,8 RG
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters == '24',"cell_type"] <- 'oIP'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('13','8'),"cell_type"] <- 'RG'

DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)


#round 3 annotation
dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
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
                                             In=c('DLX5','GAD2','SP8','SST','MEIS2','GAD1','DLX2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
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

#cluster 6,4,3,22,2,18,16,1,9,30,19,15 Ex
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('6','4','3','22','2','18','16','1','9','30','19','15'),"cell_type"] <- 'Ex'

DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#round 4 annotation
dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
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
                                             In=c('DLX5','GAD2','SP8','SST','MEIS2','GAD1','DLX2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
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

#cluster 5,26,20,11,10,0,12,7,23 In
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('5','26','20','11','10','0','12','7','23'),"cell_type"] <- 'In'

DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#the remain unknown
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('28'),"cell_type"] <- 'Unknown-1'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('29'),"cell_type"] <- 'Unknown-2'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('31'),"cell_type"] <- 'Unknown-3'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('34'),"cell_type"] <- 'Unknown-4'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$clusters %in% c('33'),"cell_type"] <- 'Unknown-5'
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

## label transfer ----------------------------------------------------------
#load data
macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')

macaque_old_seurat@meta.data[macaque_old_seurat$cell_type %in% c('oRG','vRG'),"cell_type"] <- 'RG'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('30','9','1','4','22','2','12'),'cell_type'] <- 'Ex-1'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('14','8'),'cell_type'] <- 'Ex-2'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('16','15','24','6','18'),'cell_type'] <- 'Ex-3'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('19','20'),'cell_type'] <- 'Ex-4'

macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'oRG',"sub_cell_type"] <- 'RG-1'
macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'vRG-1',"sub_cell_type"] <- 'RG-2'
macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'vRG-2',"sub_cell_type"] <- 'RG-3'

macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-1',"sub_cell_type"] <- 'Ex-1'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-2',"sub_cell_type"] <- 'Ex-2'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-3',"sub_cell_type"] <- 'Ex-3'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-4',"sub_cell_type"] <- 'Ex-4'

DimPlot(macaque_old_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
DimPlot(macaque_old_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#label transfer
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',nfeatures = 3000,npcs = 50,preprocess = TRUE)
macaque_old_seurat <- my_process_seurat(object = macaque_old_seurat,assay = 'RNA',nfeatures = 3000,npcs = 50,preprocess = TRUE)

anchor <- FindTransferAnchors(reference = macaque_old_seurat,query = macaque_RNA_seurat,dims = 1:30)
predictions <- TransferData(anchorset = anchor, refdata = macaque_old_seurat$sub_cell_type, dims = 1:30)
macaque_RNA_seurat <- AddMetaData(object = macaque_RNA_seurat,metadata = predictions)

scibet::Confusion_heatmap(ori = macaque_RNA_seurat$cell_type,prd = macaque_RNA_seurat$predicted.id)

DimPlot(macaque_RNA_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'umap')

## re-annotate -------------------------------------------------------------
p1 <- DimPlot(macaque_RNA_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'umap')
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'umap')
p1+p2

#unknown 1 Astrocyte
#round 4 annotation
dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
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
                                             In=c('DLX5','GAD2','SP8','SST','MEIS2','GAD1','DLX2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
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

macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Unknown-1',"cell_type"] <- 'Astrocyte'

#Unknown-2,4,5 delete
macaque_RNA_seurat <- macaque_RNA_seurat[,!macaque_RNA_seurat$cell_type %in% c('Unknown-2','Unknown-4','Unknown-5')]
p1 <- DimPlot(macaque_RNA_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'umap')
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'umap')
p1+p2

meta_data <- macaque_RNA_seurat@meta.data
saveRDS(meta_data,file = './res/step_7_fig_210825/meta_data_raw_210902.rds')
express_matrix <- macaque_RNA_seurat@assays$RNA@counts
saveRDS(express_matrix,file = './res/step_7_fig_210825/express_matrix_210902.rds')

# final cluster and annotation --------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './res/step_7_fig_210825/express_matrix_210902.rds')
meta_data <- readRDS(file = './res/step_7_fig_210825/meta_data_raw_210902.rds')
meta_data <- meta_data[,-c(1:7)]
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',meta.data = meta_data,min.cells = 0,min.features = 0)

#SCTransform
macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)

macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 4000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

#cluster and annotation
ndims <- 40
macaque_RNA_seurat <- FindNeighbors(object = macaque_RNA_seurat,dims = 1:ndims)
macaque_RNA_seurat <- FindClusters(object = macaque_RNA_seurat,resolution = 1.5)
macaque_RNA_seurat <- RunUMAP(object = macaque_RNA_seurat,dims = 1:ndims)

DimPlot(macaque_RNA_seurat,group.by = 'clusters',label = TRUE,repel = TRUE)

#annotate In
#add sub_cell_type
macaque_RNA_seurat$sub_cell_type <- paste(macaque_RNA_seurat$cell_type,macaque_RNA_seurat$clusters,sep = '-')

dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
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
                             group.by = 'sub_cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
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

macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type %in% c('In-31','In-26','In-11','In-0'),"cell_type"] <- 'InMGE'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type %in% c('In-10','In-5'),"cell_type"] <- 'InCGE'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type %in% c('In-7','In-23','In-20','In-12'),"cell_type"] <- 'InPSB'

DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#annotate Ex
macaque_RNA_seurat$sub_cell_type <- paste(macaque_RNA_seurat$cell_type,macaque_RNA_seurat$clusters,sep = '-')
DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
temp <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type == 'Ex']
scibet::Confusion_heatmap(ori = temp$sub_cell_type,prd = temp$predicted.id)

macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type %in% c('Ex-4'),"cell_type"] <- 'Ex-1'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type %in% c('Ex-6'),"cell_type"] <- 'Ex-2'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type %in% c('Ex-1','Ex-15','Ex-18','Ex-19','Ex-2','Ex-22','Ex-3','Ex-30'),"cell_type"] <- 'Ex-3'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type %in% c('Ex-16','Ex-9'),"cell_type"] <- 'Ex-4'


DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)


#final show
macaque_RNA_seurat$sub_cell_type <- paste(macaque_RNA_seurat$cell_type,macaque_RNA_seurat$clusters,sep = '__')

dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
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
                             group.by = 'sub_cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
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

DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#save data
meta_data <- macaque_RNA_seurat@meta.data
colnames(meta_data)
meta_data <- meta_data[,c('percent.mt','clusters','species','donor','doublet_score','doublet','cell_type')]
macaque_RNA_seurat <- macaque_RNA_seurat@assays$RNA@counts
table(colnames(macaque_RNA_seurat) == rownames(meta_data))

macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque_210903',meta.data = meta_data,min.cells = 0,min.features = 0)
macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)

macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 4000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'

macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)

ndims <- 40
macaque_RNA_seurat <- FindNeighbors(object = macaque_RNA_seurat,dims = 1:ndims)
macaque_RNA_seurat <- FindClusters(object = macaque_RNA_seurat,resolution = 1.5)
macaque_RNA_seurat <- RunUMAP(object = macaque_RNA_seurat,dims = 1:ndims)

DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
saveRDS(macaque_RNA_seurat,file = './processed_data/macaque_RNA_210903_seurat_annotated_210903.rds')

# why marker express sparse -----------------------------------------------
#load data
macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
macaque_new_seurat <- readRDS(file = './processed_data/macaque_RNA_210903_seurat_annotated_210903.rds')
greenleaf_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
human_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')

#process
human_RNA_seurat <- my_process_seurat(object = human_RNA_seurat,assay = 'RNA',nfeatures = 2000,npcs = 50,preprocess = TRUE)
human_RNA_seurat <- FindNeighbors(object = human_RNA_seurat,dims = 1:30)
human_RNA_seurat <- RunUMAP(object = human_RNA_seurat,dims = 1:30)
DimPlot(human_RNA_seurat,group.by = 'Cluster',label = TRUE,repel = TRUE)

#QC plot of macaque_new_seurat

#annotation dimplot
pdf(file = './res/step_7_fig_210904/macaque_new_seurat_dimplot.pdf',width = 8,height = 8)
DimPlot(macaque_new_seurat,group.by = 'cell_type',reduction = 'umap',label = TRUE,repel = TRUE) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5)) + 
  labs('macaque new')
dev.off()

pdf(file = './res/step_7_fig_210904/macaque_new_seurat_dimplot_split_by_donor.pdf',width = 18,height = 6)
DimPlot(macaque_new_seurat,group.by = 'cell_type',reduction = 'umap',label = FALSE,split.by = 'donor') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,size = 1,colour = 'black')) + 
  labs('macaque new')
dev.off()

#doublet
pdf(file = './res/step_7_fig_210904/macaque_new_seurat_doublet_dimplot.pdf',width = 8,height = 8)
DimPlot(macaque_new_seurat,group.by = 'doublet',reduction = 'umap',label = FALSE,cols = c('lightgrey','black')) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        plot.title = element_text(hjust = 0.5))
dev.off()

#doublet distribution
macaque_new_seurat$sub_cell_type <- paste(macaque_new_seurat$cell_type,macaque_new_seurat$clusters,sep = '__')
temp <- My_Cell_Proportion(seu.obj = macaque_new_seurat,group.by = 'doublet',split.by = 'sub_cell_type')

pdf(file = './res/step_7_fig_210904/macaque_new_seurat_doublet_sub_cell_type_proportion.pdf',width = 12,height = 6)
ggplot(temp,aes(x=sub_cell_type,y=Proportion,fill=doublet)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  scale_fill_manual(values = c('lightgrey','black')) + 
  theme_classic() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()

#donor distribution
temp <- My_Cell_Proportion(seu.obj = macaque_new_seurat,group.by = 'donor',split.by = 'sub_cell_type')

pdf(file = './res/step_7_fig_210904/macaque_new_seurat_donor_sub_cell_type_proportion.pdf',width = 12,height = 6)
ggplot(temp,aes(x=sub_cell_type,y=Proportion,fill=donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.8) + 
  theme_classic() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()

#modify macaque_old_seurat
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type %in% c('oRG','vRG'),"cell_type"] <- 'RG'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('30','9','1','4','22','2','12'),'cell_type'] <- 'Ex-1'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('14','8'),'cell_type'] <- 'Ex-2'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('16','15','24','6','18'),'cell_type'] <- 'Ex-3'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('19','20'),'cell_type'] <- 'Ex-4'

macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'oRG',"sub_cell_type"] <- 'RG-1'
macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'vRG-1',"sub_cell_type"] <- 'RG-2'
macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'vRG-2',"sub_cell_type"] <- 'RG-3'

macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-1',"sub_cell_type"] <- 'Ex-1'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-2',"sub_cell_type"] <- 'Ex-2'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-3',"sub_cell_type"] <- 'Ex-3'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-4',"sub_cell_type"] <- 'Ex-4'

DimPlot(macaque_old_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
DimPlot(macaque_old_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#label transfer
macaque_new_seurat <- my_process_seurat(object = macaque_new_seurat,assay = 'RNA',nfeatures = 2000,npcs = 50,preprocess = TRUE)
macaque_old_seurat <- my_process_seurat(object = macaque_old_seurat,assay = 'RNA',nfeatures = 2000,npcs = 50,preprocess = TRUE)

anchor <- FindTransferAnchors(reference = macaque_old_seurat,query = macaque_new_seurat,dims = 1:30)
predictions <- TransferData(anchorset = anchor, refdata = macaque_old_seurat$sub_cell_type, dims = 1:30)
macaque_new_seurat <- AddMetaData(object = macaque_new_seurat,metadata = predictions)

DimPlot(macaque_new_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'umap')

pdf(file = './res/step_7_fig_210904/macaque_new_seurat_cell_type_predict.pdf',width = 12,height = 7)
scibet::Confusion_heatmap(ori = macaque_new_seurat$sub_cell_type,prd = macaque_new_seurat$predicted.id) + 
  theme(aspect.ratio = 18/32)
dev.off()

#marker express
#macaque_new_seurat
dotplot_matrix <- my_dotplot(macaque_new_seurat,assay = 'RNA', 
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
                             group.by = 'sub_cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

pdf(file = './res/step_7_fig_210904/macaque_new_seurat_marker_express.pdf',width = 18,height = 10)
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

#macaque_old_seurat
dotplot_matrix <- my_dotplot(macaque_old_seurat,assay = 'RNA', 
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
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

pdf(file = './res/step_7_fig_210904/macaque_old_seurat_marker_express.pdf',width = 18,height = 10)
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

#greenleaf_RNA_seurat
dotplot_matrix <- my_dotplot(greenleaf_RNA_seurat,assay = 'converted', 
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
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

pdf(file = './res/step_7_fig_210904/greenleaf_RNA_seurat_marker_express.pdf',width = 18,height = 10)
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

#human_RNA_seurat
dotplot_matrix <- my_dotplot(human_RNA_seurat,assay = 'RNA', 
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
                             group.by = 'Cluster', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

pdf(file = './res/step_7_fig_210904/human_RNA_seurat_marker_express.pdf',width = 18,height = 10)
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

#sequencing influence
macaque_new_seurat$project <- 'macaque_new'
macaque_old_seurat$project <- 'macaque_old'
greenleaf_RNA_seurat$project <- 'greenleaf'
human_RNA_seurat$project <- 'PDhuman'

macaque_new_meta <- macaque_new_seurat@meta.data[,c("nCount_RNA","nFeature_RNA",'project')]
colnames(macaque_new_meta) <- c('nCount','nFeature','project')
macaque_old_meta <- macaque_old_seurat@meta.data[,c("nCount_RNA","nFeature_RNA",'project')]
colnames(macaque_old_meta) <- c('nCount','nFeature','project')
greenleaf_meta <- greenleaf_RNA_seurat@meta.data[,c("nCount_originalexp","nFeature_originalexp",'project')]
colnames(greenleaf_meta) <- c('nCount','nFeature','project')
PDhuman_meta <- human_RNA_seurat@meta.data[,c("nCount_RNA","nFeature_RNA","project")]
colnames(PDhuman_meta) <- c('nCount','nFeature','project')

meta_data <- rbind(macaque_new_meta,macaque_old_meta,greenleaf_meta,PDhuman_meta)

stat_meta_data <- matrix(nrow = 4,ncol = 3)
rownames(stat_meta_data) <- unique(meta_data$project)
colnames(stat_meta_data) <- c('nCount','nFeature','project')
stat_meta_data <- as.data.frame(stat_meta_data)
for (i in rownames(stat_meta_data)) {
  stat_meta_data[i,3] <- i
  stat_meta_data[i,1] <- round(mean(meta_data[meta_data$project == i,"nCount"]))
  stat_meta_data[i,2] <- round(mean(meta_data[meta_data$project == i,"nFeature"]))
}

p1 <- ggplot(meta_data,aes(x=project,y=nCount,fill=project)) + 
  geom_violin() + 
  ylim(c(0,15000)) + 
  geom_point(data = stat_meta_data,aes(x=project,y=nCount),size = 1,colour = 'black') + 
  geom_text(data = stat_meta_data,aes(x=project,label=nCount,y=nCount+1000)) + 
  theme_classic() + 
  theme(aspect.ratio = 0.4,
        axis.text.x = element_blank(),
        legend.position = 'none') + xlab('')

p2 <- ggplot(meta_data,aes(x=project,y=nFeature,fill=project)) + 
  geom_violin() + 
  ylim(c(0,6000)) + 
  geom_point(data = stat_meta_data,aes(x=project,y=nFeature),size = 1,colour = 'black') + 
  geom_text(data = stat_meta_data,aes(x=project,label=nFeature,y=nFeature+400)) + 
  theme_classic() + 
  theme(aspect.ratio = 0.4,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = 'none') + xlab('')

pdf(file = './res/step_7_fig_210904/combined_QC_plot.pdf',width = 6,height = 6)
p1+p2+plot_layout(nrow = 2)
dev.off()

#marker gene featureplot
p1 <- geneSetAveragePlot(genes = 'HOPX',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new HOPX')

p2 <- geneSetAveragePlot(genes = 'HOPX',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old HOPX')

pdf(file = './res/step_7_fig_210904/macaque_HOPX_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'SOX9',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new SOX9')

p2 <- geneSetAveragePlot(genes = 'SOX9',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old SOX9')

pdf(file = './res/step_7_fig_210904/macaque_SOX9_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'EOMES',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new EOMES')

p2 <- geneSetAveragePlot(genes = 'EOMES',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old EOMES')

pdf(file = './res/step_7_fig_210904/macaque_EOMES_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'SOX10',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new SOX10')

p2 <- geneSetAveragePlot(genes = 'SOX10',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old SOX10')

pdf(file = './res/step_7_fig_210904/macaque_SOX10_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'OLIG2',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new OLIG2')

p2 <- geneSetAveragePlot(genes = 'OLIG2',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old OLIG2')

pdf(file = './res/step_7_fig_210904/macaque_OLIG2_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'FEZF2',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new FEZF2')

p2 <- geneSetAveragePlot(genes = 'FEZF2',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old FEZF2')

pdf(file = './res/step_7_fig_210904/macaque_FEZF2_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'NEUROD2',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new NEUROD2')

p2 <- geneSetAveragePlot(genes = 'NEUROD2',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old NEUROD2')

pdf(file = './res/step_7_fig_210904/macaque_NEUROD2_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'SOX5',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new SOX5')

p2 <- geneSetAveragePlot(genes = 'SOX5',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old SOX5')

pdf(file = './res/step_7_fig_210904/macaque_SOX5_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'LHX6',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new LHX6')

p2 <- geneSetAveragePlot(genes = 'LHX6',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old LHX6')

pdf(file = './res/step_7_fig_210904/macaque_LHX6_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- geneSetAveragePlot(genes = 'MEIS2',object = macaque_new_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque new MEIS2')

p2 <- geneSetAveragePlot(genes = 'MEIS2',object = macaque_old_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                         scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,plot.type = 'single',print = FALSE) + 
  labs(title = 'macaque old MEIS2')

pdf(file = './res/step_7_fig_210904/macaque_MEIS2_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#more sparse gene express in macaque new


# try merge ---------------------------------------------------------------
#load data
macaque_old_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
macaque_new_seurat <- readRDS(file = './processed_data/macaque_RNA_210903_seurat_annotated_210903.rds')

macaque_old_seurat@meta.data[macaque_old_seurat$cell_type %in% c('oRG','vRG'),"cell_type"] <- 'RG'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('30','9','1','4','22','2','12'),'cell_type'] <- 'Ex-1'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('14','8'),'cell_type'] <- 'Ex-2'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('16','15','24','6','18'),'cell_type'] <- 'Ex-3'
macaque_old_seurat@meta.data[macaque_old_seurat$SCT_snn_res.1.6 %in% c('19','20'),'cell_type'] <- 'Ex-4'

macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'oRG',"sub_cell_type"] <- 'RG-1'
macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'vRG-1',"sub_cell_type"] <- 'RG-2'
macaque_old_seurat@meta.data[macaque_old_seurat$sub_cell_type == 'vRG-2',"sub_cell_type"] <- 'RG-3'

macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-1',"sub_cell_type"] <- 'Ex-1'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-2',"sub_cell_type"] <- 'Ex-2'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-3',"sub_cell_type"] <- 'Ex-3'
macaque_old_seurat@meta.data[macaque_old_seurat$cell_type == 'Ex-4',"sub_cell_type"] <- 'Ex-4'

DimPlot(macaque_old_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
DimPlot(macaque_old_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#merge
gene_list <- dplyr::intersect(rownames(macaque_new_seurat@assays$RNA@counts),rownames(macaque_old_seurat@assays$RNA@counts))

macaque_RNA_seurat <- cbind(macaque_new_seurat@assays$RNA@counts[gene_list,],macaque_old_seurat@assays$RNA@counts[gene_list,])
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',min.cells = 0,min.features = 0)

#add annotation
macaque_RNA_seurat@meta.data[,'dataset'] <- NA
macaque_RNA_seurat@meta.data[colnames(macaque_new_seurat),'dataset'] <- 'macaque_new'
macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),'dataset'] <- 'macaque_old'

macaque_RNA_seurat@meta.data[,'cell_type'] <- NA
macaque_RNA_seurat@meta.data[colnames(macaque_new_seurat),"cell_type"] <- macaque_new_seurat$cell_type
macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"cell_type"] <- macaque_old_seurat$cell_type

macaque_RNA_seurat@meta.data[,'donor'] <- NA
macaque_RNA_seurat@meta.data[colnames(macaque_new_seurat),"donor"] <- macaque_new_seurat$donor
macaque_RNA_seurat@meta.data[colnames(macaque_old_seurat),"donor"] <- macaque_old_seurat$donor

#cluster and umap
rm(list = c('macaque_new_seurat','macaque_old_seurat','gene_list'))

macaque_RNA_seurat <- NormalizeData(object = macaque_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
macaque_RNA_seurat <- SCTransform(object = macaque_RNA_seurat,variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(macaque_RNA_seurat) <- 'SCT'
macaque_RNA_seurat <- RunPCA(object = macaque_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(macaque_RNA_seurat,ndims = 50)

ndim <- 42
macaque_RNA_seurat <- FindNeighbors(object = macaque_RNA_seurat,dims = 1:ndim)
macaque_RNA_seurat <- RunUMAP(macaque_RNA_seurat,dims = 1:ndim)

p1 <- DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'dataset',label = FALSE) + theme(aspect.ratio = 1)

pdf(file = './res/step_7_fig_210904/macaque_RNA_merge_dimplot.pdf',width = 15,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()