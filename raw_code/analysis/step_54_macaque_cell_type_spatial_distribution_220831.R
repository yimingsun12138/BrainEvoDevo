#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque cell type spatial distribution                          ##
## Data: 2022.08.31                                                                ##
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
library(DESeq2)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#get layer marker
layer_marker <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')

# convert gene id to gene symbol ------------------------------------------
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

my_temp_function <- function(layer_marker,layer){
  marker <- layer_marker[[layer]]
  marker <- unique(marker)
  print('marker in annotation list')
  print(table(marker %in% macaque_annotation$gene_id))
  marker <- marker[marker %in% macaque_annotation$gene_id]
  print('marker trans id to name')
  marker <- macaque_annotation$gene_name[base::match(marker,table = macaque_annotation$gene_id)]
  print('marker in express matrix')
  print(table(marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)))
  marker <- marker[marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)]
  print('duplicated markers')
  print(table(duplicated(marker)))
  marker <- unique(marker)
  return(marker)
}

# macaque bulk data -------------------------------------------------------
#load data
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')
colnames(macaque_express_matrix) <- sub(pattern = '-',replacement = '_',colnames(macaque_express_matrix))

macaque_coldata <- data.frame(species='macaque',sample=colnames(macaque_express_matrix))
rownames(macaque_coldata) <- colnames(macaque_express_matrix)
donor_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][1]
  return(temp)
})
donor_list <- unlist(donor_list)
macaque_coldata$donor <- donor_list
layer_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][2]
  return(temp)
})
layer_list <- unlist(layer_list)
macaque_coldata$layer <- layer_list

#create deseq2 data object
layer_list <- unique(layer_list)
temp <- do.call(cbind,base::lapply(layer_list,function(x){
  temp_list <- as.character(macaque_coldata$layer == x)
  temp_list[which(temp_list == 'TRUE')] <- 'YES'
  temp_list[which(temp_list == 'FALSE')] <- 'NO'
  return(temp_list)
}))
colnames(temp) <- layer_list
rownames(temp) <- rownames(macaque_coldata)
macaque_coldata <- cbind(macaque_coldata,temp)

dds <- DESeqDataSetFromMatrix(countData = macaque_express_matrix,colData = macaque_coldata,design = ~donor+layer)
rld <- rlog(dds, blind = FALSE)

#Sample distances analysis
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)

col_fun <- colorRamp2(breaks = 1:250,colors = colorRampPalette(rev(brewer.pal(9,"Blues")))(250))
top_anno <- HeatmapAnnotation(layer = unlist(base::lapply(X = colnames(sampleDistMatrix),FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][2])
})),donor = unlist(base::lapply(X = colnames(sampleDistMatrix),FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][1])
})),which = 'column')

pdf(file = './res/step_54_fig_220831/macaque_bulk_data_sample_distance_heatmap.pdf',width = 12,height = 10)
Heatmap(matrix = sampleDistMatrix,cluster_columns = TRUE,show_column_names = TRUE,cluster_rows = TRUE,show_row_names = TRUE,
        col = col_fun,top_annotation = top_anno,width = unit(7,'inches'),height = unit(7,'inches'),name = 'distance')
dev.off()

#PCAplot
pcaData <- plotPCA(rld,intgroup = c('donor','layer'),returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar

pdf(file = './res/step_54_fig_220831/macaque_bulk_data_PCA_plot.pdf',width = 6,height = 6)
ggplot(pcaData, aes(x = PC1, y = PC2, color = layer, shape = donor)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA plot") + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5))
dev.off()

#filter duplicated layer markers
gene_list <- unlist(layer_marker)
gene_list <- gene_list[duplicated(gene_list)]
gene_list <- unique(gene_list)

for (i in names(layer_marker)) {
  temp <- layer_marker[[i]]
  temp <- temp[!(temp %in% gene_list)]
  layer_marker[[i]] <- temp
}

table(duplicated(unlist(layer_marker)))
layer_marker_vector <- do.call(what = base::c,args = base::lapply(X = names(layer_marker),FUN = function(x){
  temp <- layer_marker[[x]]
  names(temp) <- rep(x,times = length(temp))
  return(temp)
}))

#Heatmap
macaque_express_matrix <- assay(rld)[layer_marker_vector,]
macaque_express_matrix <- t(scale(t(macaque_express_matrix)))
top_anno <- HeatmapAnnotation(layer = unlist(base::lapply(X = colnames(macaque_express_matrix),FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][2])
})),which = 'column')
row_anno <- HeatmapAnnotation(layer = names(layer_marker_vector),which = 'row')

pdf(file = './res/step_54_fig_220831/macaque_bulk_data_layer_marker_heatmap.pdf',width = 8,height = 10)
Heatmap(matrix = macaque_express_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,show_row_names = FALSE,
        top_annotation = top_anno,left_annotation = row_anno,row_split = factor(names(layer_marker_vector),levels = c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ')),border = TRUE,
        column_split = factor(rld$layer,levels = c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ')),width = unit(6,'inches'),height = unit(8,'inches'))
dev.off()

# macaque multiome layer marker module score ------------------------------
for (i in c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")){
  temp_marker <- my_temp_function(layer_marker = layer_marker,layer = i)
  macaque_multiome_Seurat <- My_add_module_score(seu.obj = macaque_multiome_Seurat,assay = 'RNA',features = temp_marker,meta_var = i,scale = FALSE,center = FALSE)
}

#heatmap
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")]
meta_data <- meta_data[!(macaque_multiome_Seurat$cell_type %in% c('End','Per','Mic','VLMC')),]

#create annotation
cell_type_annotation <- HeatmapAnnotation(cell_type=macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                                          show_annotation_name = TRUE,which = 'row',
                                          name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

#plot
col_fun <- colorRamp2(c(0,0.5,1),c("#440154FF","#21908CFF","#FDE725FF"))
pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_heatmap.pdf',width = 8,height = 10)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                           levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(8,units = 'inches'),width = unit(5,units = 'inches'),
        name = 'expression',col = col_fun,border = TRUE)
dev.off()

#umap plot
#with order
p1 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'MZ',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPo',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPi',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'SP',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'IZ',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p6 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'OSVZ',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p7 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'ISVZ',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p8 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'VZ',order = TRUE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_featureplot_with_order.pdf',width = 24,height = 11)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 4)
dev.off()

#without order
p1 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'MZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPo',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPi',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'SP',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'IZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p6 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'OSVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p7 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'ISVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p8 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'VZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_featureplot_without_order.pdf',width = 24,height = 11)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 4)
dev.off()

# only select most variable 100 genes -------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
DimPlot(object = macaque_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

## create DESeq2 object ----------------------------------------------------
#load data
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')
colnames(macaque_express_matrix) <- sub(pattern = '-',replacement = '_',colnames(macaque_express_matrix))

macaque_coldata <- data.frame(species='macaque',sample=colnames(macaque_express_matrix))
rownames(macaque_coldata) <- colnames(macaque_express_matrix)
donor_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][1]
  return(temp)
})
donor_list <- unlist(donor_list)
macaque_coldata$donor <- donor_list
layer_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][2]
  return(temp)
})
layer_list <- unlist(layer_list)
macaque_coldata$layer <- layer_list

#create deseq2 data object
layer_list <- unique(layer_list)
temp <- do.call(cbind,base::lapply(layer_list,function(x){
  temp_list <- as.character(macaque_coldata$layer == x)
  temp_list[which(temp_list == 'TRUE')] <- 'YES'
  temp_list[which(temp_list == 'FALSE')] <- 'NO'
  return(temp_list)
}))
colnames(temp) <- layer_list
rownames(temp) <- rownames(macaque_coldata)
macaque_coldata <- cbind(macaque_coldata,temp)

dds <- DESeqDataSetFromMatrix(countData = macaque_express_matrix,colData = macaque_coldata,design = ~donor+layer)
rld <- rlog(dds, blind = FALSE)

## calculate DEG ------------------------------------------------------------
layer_list <- unique(macaque_coldata$layer)
DE_list <- base::lapply(layer_list,function(x){
  formula_char <- paste0('~','donor','+',as.character(x))
  formula_char <- as.formula(formula_char)
  temp <- DESeqDataSetFromMatrix(countData = macaque_express_matrix,colData = macaque_coldata,design = formula_char)
  temp <- DESeq(temp)
  temp <- results(temp,contrast = c(as.character(x),'YES','NO'),lfcThreshold = 0,alpha = 0.05)
  temp <- temp[!is.na(temp$padj),]
  temp <- temp[!is.na(temp$log2FoldChange),]
  temp <- temp[(temp$log2FoldChange > 0) & (temp$padj < 0.05) & (temp$baseMean > 1),]
  temp <- temp[order(temp$padj,decreasing = FALSE),]
  temp <- rownames(temp)[1:100]
  return(temp)
})
names(DE_list) <- layer_list
my_send_sms('DE done!')
# saveRDS(DE_list,file = './res/step_54_fig_220831/macaque_bulk_layer_marker_top_100.rds')

#filter duplicated layer markers
layer_marker <- readRDS(file = './res/step_54_fig_220831/macaque_bulk_layer_marker_top_100.rds')
gene_list <- unlist(layer_marker)
gene_list <- gene_list[duplicated(gene_list)]
gene_list <- unique(gene_list)

for (i in names(layer_marker)) {
  temp <- layer_marker[[i]]
  temp <- temp[!(temp %in% gene_list)]
  layer_marker[[i]] <- temp
}

table(duplicated(unlist(layer_marker)))
layer_marker_vector <- do.call(what = base::c,args = base::lapply(X = names(layer_marker),FUN = function(x){
  temp <- layer_marker[[x]]
  names(temp) <- rep(x,times = length(temp))
  return(temp)
}))

#Heatmap
macaque_express_matrix <- assay(rld)[layer_marker_vector,]
macaque_express_matrix <- t(scale(t(macaque_express_matrix)))
top_anno <- HeatmapAnnotation(layer = unlist(base::lapply(X = colnames(macaque_express_matrix),FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][2])
})),which = 'column')
row_anno <- HeatmapAnnotation(layer = names(layer_marker_vector),which = 'row')

pdf(file = './res/step_54_fig_220831/macaque_bulk_data_layer_marker_heatmap_only_100_genes.pdf',width = 6,height = 8)
Heatmap(matrix = macaque_express_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,show_row_names = FALSE,
        top_annotation = top_anno,left_annotation = row_anno,row_split = factor(names(layer_marker_vector),levels = c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ')),border = TRUE,
        column_split = factor(rld$layer,levels = c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ')),width = unit(4,'inches'),height = unit(6,'inches'))
dev.off()

## convert gene id to symbol -----------------------------------------------
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

my_temp_function <- function(layer_marker,layer){
  marker <- layer_marker[[layer]]
  marker <- unique(marker)
  print('marker in annotation list')
  print(table(marker %in% macaque_annotation$gene_id))
  marker <- marker[marker %in% macaque_annotation$gene_id]
  print('marker trans id to name')
  marker <- macaque_annotation$gene_name[base::match(marker,table = macaque_annotation$gene_id)]
  print('marker in express matrix')
  print(table(marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)))
  marker <- marker[marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)]
  print('duplicated markers')
  print(table(duplicated(marker)))
  marker <- unique(marker)
  return(marker)
}

## macaque multiome layer marker module score ------------------------------
for (i in c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")){
  temp_marker <- my_temp_function(layer_marker = layer_marker,layer = i)
  macaque_multiome_Seurat <- My_add_module_score(seu.obj = macaque_multiome_Seurat,assay = 'RNA',features = temp_marker,meta_var = i,scale = FALSE,center = FALSE)
}

#heatmap
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")]
meta_data <- meta_data[!(macaque_multiome_Seurat$cell_type %in% c('End','Per','Mic','VLMC')),]

#create annotation
cell_type_annotation <- HeatmapAnnotation(cell_type=macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                                          show_annotation_name = TRUE,which = 'row',
                                          name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

#plot
col_fun <- colorRamp2(c(0,0.5,1),c("#440154FF","#21908CFF","#FDE725FF"))
pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_heatmap_only_100_genes.pdf',width = 8,height = 10)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                           levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(8,units = 'inches'),width = unit(5,units = 'inches'),
        name = 'expression',col = col_fun,border = TRUE)
dev.off()

#umap plot
#with order
p1 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'MZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPo',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPi',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'SP',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'IZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p6 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'OSVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p7 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'ISVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p8 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'VZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_featureplot_only_100_genes_without_order.pdf',width = 24,height = 11)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 4)
dev.off()

# Seurat add module score -------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
layer_marker <- readRDS(file = './res/step_54_fig_220831/macaque_bulk_layer_marker_top_100.rds')

#convert id to symbol
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

my_temp_function <- function(layer_marker,layer){
  marker <- layer_marker[[layer]]
  marker <- unique(marker)
  print('marker in annotation list')
  print(table(marker %in% macaque_annotation$gene_id))
  marker <- marker[marker %in% macaque_annotation$gene_id]
  print('marker trans id to name')
  marker <- macaque_annotation$gene_name[base::match(marker,table = macaque_annotation$gene_id)]
  print('marker in express matrix')
  print(table(marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)))
  marker <- marker[marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)]
  print('duplicated markers')
  print(table(duplicated(marker)))
  marker <- unique(marker)
  return(marker)
}

#add module score
for (i in c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")){
  temp_marker <- my_temp_function(layer_marker = layer_marker,layer = i)
  macaque_multiome_Seurat <- Seurat::AddModuleScore(object = macaque_multiome_Seurat,features = list(temp_marker),name = i)
}

#heatmap
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("MZ1","CPo1","CPi1","SP1","IZ1","OSVZ1","ISVZ1","VZ1")]
colnames(meta_data) <- sub(pattern = '1',replacement = '',x = colnames(meta_data))
meta_data <- meta_data[!(macaque_multiome_Seurat$cell_type %in% c('End','Per','Mic','VLMC')),]

#create annotation
cell_type_annotation <- HeatmapAnnotation(cell_type=macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                                          show_annotation_name = TRUE,which = 'row',
                                          name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

#plot
col_fun <- colorRamp2(c(-0.2,0,0.2),c("#440154FF","#21908CFF","#FDE725FF"))
pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_heatmap_only_100_genes_seurat_module_score.pdf',width = 8,height = 10)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                           levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(8,units = 'inches'),width = unit(5,units = 'inches'),
        name = 'expression',col = col_fun,border = TRUE)
dev.off()

#umap plot
#with order
colnames(macaque_multiome_Seurat@meta.data) <- sub(pattern = '1',replacement = '',x = colnames(macaque_multiome_Seurat@meta.data))
p1 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'MZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPo',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPi',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'SP',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'IZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p6 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'OSVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p7 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'ISVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p8 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'VZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_featureplot_only_100_genes_without_order_Seurat_modulescore.pdf',width = 24,height = 11)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 4)
dev.off()

# Seurat add module score all marker -------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
layer_marker <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')

#convert id to symbol
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

my_temp_function <- function(layer_marker,layer){
  marker <- layer_marker[[layer]]
  marker <- unique(marker)
  print('marker in annotation list')
  print(table(marker %in% macaque_annotation$gene_id))
  marker <- marker[marker %in% macaque_annotation$gene_id]
  print('marker trans id to name')
  marker <- macaque_annotation$gene_name[base::match(marker,table = macaque_annotation$gene_id)]
  print('marker in express matrix')
  print(table(marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)))
  marker <- marker[marker %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)]
  print('duplicated markers')
  print(table(duplicated(marker)))
  marker <- unique(marker)
  return(marker)
}

#add module score
for (i in c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")){
  temp_marker <- my_temp_function(layer_marker = layer_marker,layer = i)
  macaque_multiome_Seurat <- Seurat::AddModuleScore(object = macaque_multiome_Seurat,features = list(temp_marker),name = i)
}

#heatmap
meta_data <- macaque_multiome_Seurat@meta.data
meta_data <- meta_data[,c("MZ1","CPo1","CPi1","SP1","IZ1","OSVZ1","ISVZ1","VZ1")]
colnames(meta_data) <- sub(pattern = '1',replacement = '',x = colnames(meta_data))
meta_data <- meta_data[!(macaque_multiome_Seurat$cell_type %in% c('End','Per','Mic','VLMC')),]

#create annotation
cell_type_annotation <- HeatmapAnnotation(cell_type=macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                                          show_annotation_name = TRUE,which = 'row',
                                          name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

#plot
col_fun <- colorRamp2(c(-0.2,0,0.2),c("#440154FF","#21908CFF","#FDE725FF"))
pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_heatmap_genes_seurat_module_score.pdf',width = 8,height = 10)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(macaque_multiome_Seurat@meta.data[rownames(meta_data),"cell_type"],
                           levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(8,units = 'inches'),width = unit(5,units = 'inches'),
        name = 'expression',col = col_fun,border = TRUE)
dev.off()

#umap plot
#with order
colnames(macaque_multiome_Seurat@meta.data) <- sub(pattern = '1',replacement = '',x = colnames(macaque_multiome_Seurat@meta.data))
p1 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'MZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPo',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'CPi',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'SP',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'IZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p6 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'OSVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p7 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'ISVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p8 <- FeaturePlot(object = macaque_multiome_Seurat,features = 'VZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_54_fig_220831/macaque_multiome_Seurat_layer_marker_featureplot_genes_without_order_Seurat_modulescore.pdf',width = 24,height = 11)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 4)
dev.off()

# Cibersortx --------------------------------------------------------------
#load single cell data
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
table(macaque_multiome_Seurat$cell_type)

#comprehensive matrix -- down sample to 1000
cell_list <- c()
for (i in unique(macaque_multiome_Seurat$cell_type)) {
  if(i %in% c('End','Per','VLMC')){} else{
    temp <- colnames(macaque_multiome_Seurat)[macaque_multiome_Seurat$cell_type %in% c(i)]
    if(length(temp) <= 1000){
      cell_list <- append(cell_list,temp)
    } else{
      temp <- sample(x = temp,size = 1000,replace = FALSE)
      cell_list <- append(cell_list,temp)
    }
  }
}

express_matrix <- macaque_multiome_Seurat[,cell_list]
express_matrix <- express_matrix@assays$RNA@counts
express_matrix <- express_matrix[rowSums(express_matrix) > 0,]
express_matrix <- as.data.frame(express_matrix)
temp <- data.frame(cell = colnames(express_matrix),Gene = macaque_multiome_Seurat@meta.data[colnames(express_matrix),"cell_type"])
temp <- t(temp)
colnames(temp) <- colnames(express_matrix)
temp <- temp['Gene',,drop = FALSE]
express_matrix <- rbind(temp,express_matrix)
colnames(express_matrix) <- NULL

write.table(x = express_matrix,file = './res/step_54_fig_220831/macaque_multiome_RNA_expression_matrix.txt',quote = FALSE,sep = '\t',col.names = FALSE)

#bulk data
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')
colnames(macaque_express_matrix) <- sub(pattern = '-',replacement = '_',colnames(macaque_express_matrix))

macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_name','gene_id')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]
macaque_annotation <- cbind(macaque_annotation,macaque_annotation)

macaque_express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = macaque_express_matrix,anno = macaque_annotation,filter_anno = FALSE,future.globals.maxSize = 200*(1024^3),workers = 5)
macaque_express_matrix <- macaque_express_matrix[rowSums(macaque_express_matrix) > 0,]
macaque_express_matrix <- as.data.frame(macaque_express_matrix)
temp <- data.frame(Gene = colnames(macaque_express_matrix),layer = colnames(macaque_express_matrix))
temp <- t(temp)
colnames(temp) <- colnames(macaque_express_matrix)
temp <- temp[1,,drop = FALSE]
macaque_express_matrix <- rbind(temp,macaque_express_matrix)
colnames(macaque_express_matrix) <- NULL
write.table(x = macaque_express_matrix,file = './res/step_54_fig_220831/macaque_bulk_layer_RNA_expression_matrix.txt',quote = FALSE,sep = '\t',col.names = FALSE)

# boxplot asked by liuyt --------------------------------------------------
#load bulk data
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')
colnames(macaque_express_matrix) <- sub(pattern = '-',replacement = '_',colnames(macaque_express_matrix))

macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_name','gene_id')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]
macaque_annotation <- cbind(macaque_annotation,macaque_annotation)

macaque_express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = macaque_express_matrix,anno = macaque_annotation,filter_anno = FALSE,future.globals.maxSize = 200*(1024^3),workers = 5)
macaque_express_matrix <- macaque_express_matrix[rowSums(macaque_express_matrix) > 0,]

#create DESeq2 dataset
macaque_coldata <- data.frame(species='macaque',sample=colnames(macaque_express_matrix))
rownames(macaque_coldata) <- colnames(macaque_express_matrix)
donor_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][1]
  return(temp)
})
donor_list <- unlist(donor_list)
macaque_coldata$donor <- donor_list
layer_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][2]
  return(temp)
})
layer_list <- unlist(layer_list)
macaque_coldata$layer <- layer_list

#create deseq2 data object
layer_list <- unique(layer_list)
temp <- do.call(cbind,base::lapply(layer_list,function(x){
  temp_list <- as.character(macaque_coldata$layer == x)
  temp_list[which(temp_list == 'TRUE')] <- 'YES'
  temp_list[which(temp_list == 'FALSE')] <- 'NO'
  return(temp_list)
}))
colnames(temp) <- layer_list
rownames(temp) <- rownames(macaque_coldata)
macaque_coldata <- cbind(macaque_coldata,temp)

dds <- DESeqDataSetFromMatrix(countData = macaque_express_matrix,colData = macaque_coldata,design = ~donor+layer)
rld <- rlog(dds, blind = FALSE)

#EOMES
gene_list <- 'EOMES'
temp <- data.frame(express = as.numeric(assay(rld)[gene_list,]),sample = rld$sample,donor = rld$donor,layer = rld$layer)
temp$layer <- factor(temp$layer,levels = c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ'))

char <- paste('./res/step_54_fig_220831/macaque_bulk_layer',gene_list,'express_boxplot.pdf',sep = '_')
pdf(file = char,width = 7,height = 5)
ggplot(data = temp,aes(x = layer,y=express,fill = layer)) + 
  geom_boxplot(outlier.alpha = 0) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.6,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        axis.line = element_blank(),
        axis.title.x = element_blank()) + 
  scale_fill_manual(values = ggsci::pal_npg()(8)) + 
  labs(title = gene_list)
dev.off()

#create fig by for loop
for (i in c("EOMES","CORO1C","SMOC1","GRAMD1B","PPP1R17","RASGRP1","NEUROD4","NHLH1","TMEM158","PENK")) {
  gene_list <- i
  temp <- data.frame(express = as.numeric(assay(rld)[gene_list,]),sample = rld$sample,donor = rld$donor,layer = rld$layer)
  temp$layer <- factor(temp$layer,levels = c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ'))
  
  char <- paste('./res/step_54_fig_220831/macaque_bulk_layer',gene_list,'express_boxplot.pdf',sep = '_')
  p <- ggplot(data = temp,aes(x = layer,y=express,fill = layer)) + 
    geom_boxplot(outlier.alpha = 0) + 
    theme_cowplot() + 
    theme(aspect.ratio = 0.6,
          axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
          plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = NA,colour = 'black',size = 1),
          axis.line = element_blank(),
          axis.title.x = element_blank()) + 
    scale_fill_manual(values = ggsci::pal_npg()(8)) + 
    labs(title = gene_list)
  pdf(file = char,width = 7,height = 5)
  print(p)
  dev.off()
}

# human cell type layer distribution --------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/220718_summary/Greenleaf_RNA_HVG2kPC40_withPredictedAnno_Seurat_220721.rds')
DimPlot(object = Greenleaf_RNA_Seurat,group.by = 'ReAnno_celltype',label = TRUE,repel = TRUE)

#marker gene list
layer_marker <- list()
layer_marker[['MZ']] <- c('CALB2','B4GALNT4','CBX6','CCDC85C','CXCL14','IGF1','LRP3','MLLT6','WTIP')
layer_marker[['CPo']] <- c('MPC1','DACT1','GALNT3','MAN1C1','MYCN','RHOU','THG1L','TMOD1')
layer_marker[['CPi']] <- c('LMO4','CD1D','CYTIP','FOXP1','HS3ST2','PCDH19','PTPRK','SAMD5','TRAF5')
layer_marker[['SP']] <- c('ADTRP','CDH18','HAS3','HS3ST4','NXPH3','SLC35F3','TMEM178A','ZFPM2')
layer_marker[['IZ']] <- c('ADM','COL20A1','NEU4','NKX2-2','OLIG1','PDGFRA','PPP1R3G','SP5')
layer_marker[['OSVZ']] <- c('ACSBG1','BTBD17','PPP1R17','LYN','PLAC9','PLAGL1','PREX2','TNC','EOMES')
layer_marker[['ISVZ']] <- c('PAX6','CDK6','GAS1','JAG1','MIR99AHG','PTH2R','SETD7')
layer_marker[['VZ']] <- c('ST13','MGARP','CA3','CASP6','CD44','GLT8D2','PGAM2','POU4F1','SFRP2','GFAP')

#check marker 
table(layer_marker[['MZ']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))
table(layer_marker[['CPo']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))
table(layer_marker[['CPi']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))
table(layer_marker[['SP']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))
table(layer_marker[['IZ']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))
table(layer_marker[['OSVZ']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))
table(layer_marker[['ISVZ']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))
table(layer_marker[['VZ']] %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))

#add module score
for (i in names(layer_marker)) {
  temp_marker <- layer_marker[[i]]
  Greenleaf_RNA_Seurat <- Seurat::AddModuleScore(object = Greenleaf_RNA_Seurat,features = list(temp_marker),name = i)
}

#heatmap
meta_data <- Greenleaf_RNA_Seurat@meta.data
meta_data <- meta_data[,c("MZ1","CPo1","CPi1","SP1","IZ1","OSVZ1","ISVZ1","VZ1")]
colnames(meta_data) <- sub(pattern = '1',replacement = '',colnames(meta_data))
meta_data <- meta_data[!(Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('End','Per','Mic','VLMC')),]

#create annotation
cell_type_annotation <- HeatmapAnnotation(cell_type=Greenleaf_RNA_Seurat@meta.data[rownames(meta_data),"ReAnno_celltype"],
                                          show_annotation_name = TRUE,which = 'row',
                                          name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

#plot
col_fun <- colorRamp2(c(-0.3,0,0.3),c("#440154FF","#21908CFF","#FDE725FF"))
pdf(file = './res/step_54_fig_220831/Greenleaf_RNA_Seurat_layer_marker_heatmap.pdf',width = 8,height = 10)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(Greenleaf_RNA_Seurat@meta.data[rownames(meta_data),"ReAnno_celltype"],
                           levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(8,units = 'inches'),width = unit(5,units = 'inches'),
        name = 'expression',col = col_fun,border = TRUE)
dev.off()

#umap plot
colnames(Greenleaf_RNA_Seurat@meta.data) <- sub(pattern = '1$',replacement = '',x = colnames(Greenleaf_RNA_Seurat@meta.data),fixed = FALSE)
#without order
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'MZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'CPo',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p3 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'CPi',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p4 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'SP',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p5 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'IZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p6 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'OSVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p7 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'ISVZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)
p8 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'VZ',order = FALSE,max.cutoff = 'q99',cols = c('#e0ecf4','#8c6bb1','#4d004b')) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_54_fig_220831/Greenleaf_RNA_Seurat_layer_marker_featureplot_without_order.pdf',width = 24,height = 11)
p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol = 4)
dev.off()

# macaque bulk cibersortx analysis ----------------------------------------
#load data
Cibersortx_results <- read.csv(file = './res/step_54_fig_220831/macaque_bulk_CIBERSORTx_Results.csv',row.names = 1)
colnames(Cibersortx_results) <- sub(pattern = '.',replacement = '-',x = colnames(Cibersortx_results),fixed = TRUE)

#generate ggplot based data.frame
bar_table <- Cibersortx_results[,1:(length(Cibersortx_results) - 3)]
bar_table <- base::do.call(what = rbind,args = base::lapply(X = colnames(bar_table),FUN = function(x){
  temp <- data.frame(sample = rownames(bar_table),cell_type = x,proportion = bar_table[,x])
  return(temp)
}))
bar_table$layer <- unlist(base::lapply(X = bar_table$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][2])
}))
bar_table$donor <- unlist(base::lapply(X = bar_table$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][1])
}))

#re-create color param table
# col_param <- read.csv(file = './data/parameter/scRNAseq_color_param.csv',row.names = 1)
# col_param[col_param$term == 'Ex-3',"term"] <- 'Ex-4'
# col_param[col_param$term == 'Ex-2',"term"] <- 'Ex-3'
# col_param[col_param$term == 'Ependymal',"term"] <- 'VLMC'
# col_param <- rbind(col_param,data.frame(term = 'Ex-2',color = '#688EC1'))
# col_param <- rbind(col_param,data.frame(term = 'Cycling',color = '#DC494C'))
# write.csv(x = col_param,file = './data/parameter/scRNAseq_color_param.csv')
col_param <- read.csv(file = './data/parameter/scRNAseq_color_param.csv')
temp_col <- col_param$color
names(temp_col) <- col_param$term

#ggplot
integrated_bar_table <- bar_table
integrated_bar_table$proportion <- (integrated_bar_table$proportion)/4
integrated_bar_table$layer <- factor(integrated_bar_table$layer,levels = c('VZ','ISVZ','OSVZ','IZ','SP','CPi','CPo','MZ'))
integrated_bar_table$cell_type <- factor(integrated_bar_table$cell_type,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC','Mic'))

pdf(file = './res/step_54_fig_220831/macaque_bulk_CIBERSORTx_cell_type_layer_distribution.pdf',width = 8,height = 5)
ggplot(data = integrated_bar_table, mapping = aes(x = layer,y = proportion,fill = cell_type))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.6) + 
  scale_fill_manual(values = temp_col[c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC','Mic')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5)
dev.off()

#RG-2 ratio in all sample
ggplot(bar_table[bar_table$cell_type == 'RG-2',],aes(x = layer,y = proportion,fill = donor)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.6) + 
  theme_cowplot()
