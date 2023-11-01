#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Greenleaf scRNAseq process and reannotate                       ##
## Data: 2022.09.16                                                                ##
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
library(ggpubr)
library(ggtext)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
liuyt_human_data <- readRDS(file = './processed_data/220718_summary/Greenleaf_RNA_HVG2kPC40_withPredictedAnno_Seurat_220721.rds')

# re-process --------------------------------------------------------------
liuyt_human_data <- FindClusters(object = liuyt_human_data,resolution = seq(from = 0.3,to = 2,by = 0.1))

# display -----------------------------------------------------------------
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype

for (i in seq(from = 0.3,to = 2,by = 0.1)) {
  char <- paste('RNA_snn_res',as.character(i),sep = '.')
  p <- scibet::Confusion_heatmap(ori = liuyt_human_data@meta.data[,char],prd = liuyt_human_data$ReAnno_celltype)
  char <- paste0('./res/step_59_fig_220916/human_RNA_Seurat_cluster_res_',as.character(i),'_cell_type_confusion_heatmap.pdf')
  pdf(file = char,width = 8,height = 6)
  print(p)
  dev.off()
}

#cluster res 0.7
p1 <- DimPlot(object = liuyt_human_data,group.by = 'RNA_snn_res.0.7',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = liuyt_human_data,group.by = 'ReAnno_celltype',cols = temp_col,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#cluster res 0.8
p1 <- DimPlot(object = liuyt_human_data,group.by = 'RNA_snn_res.0.8',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = liuyt_human_data,group.by = 'ReAnno_celltype',cols = temp_col,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)

#cluster res 0.85
liuyt_human_data <- FindClusters(object = liuyt_human_data,resolution = 0.85)
scibet::Confusion_heatmap(ori = liuyt_human_data$RNA_snn_res.0.85,prd = liuyt_human_data$ReAnno_celltype)
p1 <- DimPlot(object = liuyt_human_data,group.by = 'RNA_snn_res.0.85',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = liuyt_human_data,group.by = 'ReAnno_celltype',cols = temp_col,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+plot_layout(ncol = 2)
#dertermined

# re-process and annotate -------------------------------------------------
#load data
liuyt_human_data <- readRDS(file = './processed_data/220718_summary/Greenleaf_RNA_HVG2kPC40_withPredictedAnno_Seurat_220721.rds')
liuyt_human_data <- FindClusters(object = liuyt_human_data,resolution = 0.85)
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype

#re-annotate
p1 <- DimPlot(object = liuyt_human_data,group.by = 'RNA_snn_res.0.85',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = liuyt_human_data,group.by = 'ReAnno_celltype',cols = temp_col,label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(object = liuyt_human_data,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p1+p2+p3+plot_layout(ncol = 3)
scibet::Confusion_heatmap(ori = liuyt_human_data$RNA_snn_res.0.85,prd = liuyt_human_data$ReAnno_celltype)
scibet::Confusion_heatmap(ori = liuyt_human_data$RNA_snn_res.0.85,prd = liuyt_human_data$cell_type)

liuyt_human_data$own_anno <- NA
#Cycling
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('12','17','27'),"own_anno"] <- 'Cycling'

#End
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('29'),"own_anno"] <- 'End'

#InCGE
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('4','14','21','25'),"own_anno"] <- 'InCGE'

#InMGE
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('5'),"own_anno"] <- 'InMGE'

#IP
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('7','11','20'),"own_anno"] <- 'IP'

#Mic
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('24','26'),"own_anno"] <- 'Mic'

#OPC
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('23'),"own_anno"] <- 'OPC'

#Per
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('28'),"own_anno"] <- 'per'

#tRG
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('22'),"own_anno"] <- 'tRG'

#oIPC
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('18'),"own_anno"] <- 'oIPC'

#Early RG
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('6','16','30'),"own_anno"] <- 'Early RG'

#Late RG
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('9'),"own_anno"] <- 'Late RG'

#Ex
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('3'),"own_anno"] <- 'Ex-1'
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('0'),"own_anno"] <- 'Ex-2'
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('2'),"own_anno"] <- 'Ex-3'
liuyt_human_data@meta.data[(liuyt_human_data$RNA_snn_res.0.85 %in% c('8') & liuyt_human_data@reductions$umap@cell.embeddings[,2] < 7),"own_anno"] <- 'Ex-4'
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('1'),"own_anno"] <- 'Ex-5'
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('10'),"own_anno"] <- 'Ex-6'
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('19'),"own_anno"] <- 'Ex-7'
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('15'),"own_anno"] <- 'Ex-8'
liuyt_human_data@meta.data[(liuyt_human_data$RNA_snn_res.0.85 %in% c('8') & liuyt_human_data@reductions$umap@cell.embeddings[,2] >= 7),"own_anno"] <- 'Ex-9'
liuyt_human_data@meta.data[liuyt_human_data$RNA_snn_res.0.85 %in% c('13'),"own_anno"] <- 'Ex-10'

# save data ---------------------------------------------------------------
saveRDS(object = liuyt_human_data,file = './processed_data/220802_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

# try cluster cell type ---------------------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

#aggregate expression matrix
express_matrix <- Seurat::AggregateExpression(object = Greenleaf_RNA_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'own_anno',slot = 'data',verbose = TRUE)
express_matrix <- express_matrix$RNA
col_name <- colnames(express_matrix)
row_name <- rownames(express_matrix)
express_matrix <- base::do.call(what = cbind,args = base::lapply(X = colnames(express_matrix),FUN = function(x){
  temp <- express_matrix[,x]
  temp <- temp/sum(temp)*1000000
  return(temp)
}))
colnames(express_matrix) <- col_name
rownames(express_matrix) <- row_name
express_matrix <- log1p(express_matrix)
express_matrix <- express_matrix[VariableFeatures(Greenleaf_RNA_Seurat),]

#dist
dist_matrix <- dist(t(express_matrix))
cl <- hclust(d = dist_matrix,method = 'ward.D')
cl <- as.dendrogram(cl)
cl <- reorder(x = cl,rep(1,times = 22))

pdf(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_anno_clustering.pdf',width = 8,height = 3)
plot(cl)
dev.off()

saveRDS(object = cl,file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_anno_clustering.rds')

# find markers of all cell type -------------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

#call marker
cell_type_list <- unique(Greenleaf_RNA_Seurat$own_anno)
marker_list <- do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  other_cell_list <- cell_type_list[cell_type_list != x]
  temp <- my_seurat_marker_wilcox_test(seu.obj = Greenleaf_RNA_Seurat,assay = 'RNA',ident.1 = x,ident.2 = other_cell_list,group.by = 'own_anno',min.express = 0,log2fc_thresholf = 0.5,pct.1 = 0.4,only.pos = TRUE,workers = 2,future.globals.maxSize = 200*(1024^3))
  temp$gene <- rownames(temp)
  temp$cluster <- x
  gc()
  print(paste(x,'done!',sep = ' '))
  return(temp)
}))

my_send_sms('call marker done!')

saveRDS(object = marker_list,file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')

# process microarray data -------------------------------------------------
#load data
pcw_1516 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_table.rds')
pcw_1516_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_sampleinfor_table.rds')
pcw_21 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_table.rds')
pcw_21_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_sampleinfor_table.rds')

#express_matrix
table(rownames(pcw_1516) == rownames(pcw_21))
express_matrix <- cbind(pcw_1516,pcw_21)
meta_data <- data.frame(region = c(pcw_1516_meta_data$region,pcw_21_meta_data$region),age = c(pcw_1516_meta_data$age,pcw_21_meta_data$age))

#average express matrix
layer_list <- c('MZ','CPo','CPi','SP','IZ','SZo','SZi','VZ')
layer_list %in% meta_data$region

express_matrix <- base::do.call(what = cbind,args = base::lapply(X = layer_list,FUN = function(x){
  temp <- express_matrix[,meta_data$region == x,drop = FALSE]
  temp <- rowMeans(temp)
  return(temp)
}))
colnames(express_matrix) <- layer_list
rownames(express_matrix) <- rownames(pcw_1516)

#marker gene
marker_list <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')
marker_list <- marker_list[marker_list$gene %in% rownames(express_matrix),]
table(marker_list$cluster)
marker_list$diff <- marker_list$mean1 - marker_list$mean2
marker_list <- marker_list[order(marker_list$diff,decreasing = TRUE),]

cell_type_list <- unique(marker_list$cluster)
marker_list <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(nrow(temp) <= 200){
    return(temp)
  }else{
    return(temp[1:200,])
  }
}))

#add module score
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
for (i in cell_type_list) {
  microarray_Seurat <- Seurat::AddModuleScore(object = microarray_Seurat,features = list((marker_list$gene)[marker_list$cluster == i]),assay = 'RNA',name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,grep(pattern = '1$',x = colnames(microarray_Seurat@meta.data),fixed = FALSE)]
colnames(meta_data) <- sub(pattern = '1$',replacement = '',x = colnames(meta_data),fixed = FALSE)
colnames(meta_data) <- sub(pattern = '.',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/average_all_microarray_cell_type_marker_seurat_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

#try my add module score
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
microarray_Seurat@assays$RNA@data <- log1p(express_matrix)
for (i in cell_type_list) {
  microarray_Seurat <- My_add_module_score(seu.obj = microarray_Seurat,assay = 'RNA',features = (marker_list$gene)[marker_list$cluster == i],meta_var = i,scale = TRUE,center = TRUE)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,cell_type_list]
colnames(meta_data) <- sub(pattern = ' ',replacement = '_',x = colnames(meta_data),fixed = TRUE)
colnames(meta_data) <- sub(pattern = '-',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/average_all_microarray_cell_type_marker_my_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

# use pcw16 microarray data -----------------------------------------------
#load data
pcw_1516 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_table.rds')
pcw_1516_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_sampleinfor_table.rds')

temp <- which(grepl(pattern = 'dorsolateral prefrontal cortex',x = pcw_1516_meta_data$structure_name,fixed = TRUE) & pcw_1516_meta_data$age == 'pcw15')
pcw_1516 <- pcw_1516[,temp]
pcw_1516_meta_data <- pcw_1516_meta_data[temp,]
colnames(pcw_1516) <- pcw_1516_meta_data$region
express_matrix <- pcw_1516

#marker gene
marker_list <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')
marker_list <- marker_list[marker_list$gene %in% rownames(express_matrix),]
table(marker_list$cluster)
marker_list$diff <- marker_list$mean1 - marker_list$mean2
marker_list <- marker_list[order(marker_list$diff,decreasing = TRUE),]

cell_type_list <- unique(marker_list$cluster)
marker_list <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(nrow(temp) <= 200){
    return(temp)
  }else{
    return(temp[1:200,])
  }
}))

#add module score
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
for (i in cell_type_list) {
  microarray_Seurat <- Seurat::AddModuleScore(object = microarray_Seurat,features = list((marker_list$gene)[marker_list$cluster == i]),assay = 'RNA',name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
layer_list <- c('CPo','CPi','SP','SZo','SZi','VZ')
meta_data <- microarray_Seurat@meta.data[,grep(pattern = '1$',x = colnames(microarray_Seurat@meta.data),fixed = FALSE)]
colnames(meta_data) <- sub(pattern = '1$',replacement = '',x = colnames(meta_data),fixed = FALSE)
colnames(meta_data) <- sub(pattern = '.',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/pcw16_microarray_cell_type_marker_seurat_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

#try my add module score
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
microarray_Seurat@assays$RNA@data <- as.matrix(log1p(express_matrix))
for (i in cell_type_list) {
  microarray_Seurat <- My_add_module_score(seu.obj = microarray_Seurat,assay = 'RNA',features = (marker_list$gene)[marker_list$cluster == i],meta_var = i,scale = TRUE,center = TRUE)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,cell_type_list]
colnames(meta_data) <- sub(pattern = ' ',replacement = '_',x = colnames(meta_data),fixed = TRUE)
colnames(meta_data) <- sub(pattern = '-',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/pcw16_microarray_cell_type_marker_my_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

# use pcw21 microarray data -----------------------------------------------
#load data
pcw_21 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_table.rds')
pcw_21_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_sampleinfor_table.rds')

temp <- which(grepl(pattern = 'dorsolateral prefrontal cortex',x = pcw_21_meta_data$structure_name,fixed = TRUE) & pcw_21_meta_data$age == 'pcw21_03')
pcw_21 <- pcw_21[,temp]
pcw_21_meta_data <- pcw_21_meta_data[temp,]
colnames(pcw_21) <- pcw_21_meta_data$region
express_matrix <- pcw_21

#marker gene
marker_list <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')
marker_list <- marker_list[marker_list$gene %in% rownames(express_matrix),]
table(marker_list$cluster)
marker_list$diff <- marker_list$mean1 - marker_list$mean2
marker_list <- marker_list[order(marker_list$diff,decreasing = TRUE),]

cell_type_list <- unique(marker_list$cluster)
marker_list <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(nrow(temp) <= 200){
    return(temp)
  }else{
    return(temp[1:200,])
  }
}))

#add module score
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
for (i in cell_type_list) {
  microarray_Seurat <- Seurat::AddModuleScore(object = microarray_Seurat,features = list((marker_list$gene)[marker_list$cluster == i]),assay = 'RNA',name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
layer_list <- c('MZ','CPo','CPi','SP','IZ','SZo','SZi','VZ')
meta_data <- microarray_Seurat@meta.data[,grep(pattern = '1$',x = colnames(microarray_Seurat@meta.data),fixed = FALSE)]
colnames(meta_data) <- sub(pattern = '1$',replacement = '',x = colnames(meta_data),fixed = FALSE)
colnames(meta_data) <- sub(pattern = '.',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/pcw21_microarray_cell_type_marker_seurat_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

#try my add module score
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
microarray_Seurat@assays$RNA@data <- as.matrix(log1p(express_matrix))
for (i in cell_type_list) {
  microarray_Seurat <- My_add_module_score(seu.obj = microarray_Seurat,assay = 'RNA',features = (marker_list$gene)[marker_list$cluster == i],meta_var = i,scale = TRUE,center = TRUE)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,cell_type_list]
colnames(meta_data) <- sub(pattern = ' ',replacement = '_',x = colnames(meta_data),fixed = TRUE)
colnames(meta_data) <- sub(pattern = '-',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/pcw21_microarray_cell_type_marker_my_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()


# try scale first then average --------------------------------------------
temp_add_module_score <- function(seu.obj,gene_list,name){
  temp <- seu.obj@assays$RNA@counts[gene_list,]
  temp <- t(base::scale(x = t(temp),center = TRUE,scale = TRUE))
  temp <- colMeans(x = temp)
  seu.obj@meta.data[,name] <- temp
  return(seu.obj)
}


## pcw 16 ------------------------------------------------------------------
#load data
pcw_1516 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_table.rds')
pcw_1516_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_sampleinfor_table.rds')

temp <- which(grepl(pattern = 'dorsolateral prefrontal cortex',x = pcw_1516_meta_data$structure_name,fixed = TRUE) & pcw_1516_meta_data$age == 'pcw15')
pcw_1516 <- pcw_1516[,temp]
pcw_1516_meta_data <- pcw_1516_meta_data[temp,]
colnames(pcw_1516) <- pcw_1516_meta_data$region
express_matrix <- pcw_1516

#marker gene
marker_list <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')
marker_list <- marker_list[marker_list$gene %in% rownames(express_matrix),]
table(marker_list$cluster)
marker_list$diff <- marker_list$mean1 - marker_list$mean2
marker_list <- marker_list[order(marker_list$diff,decreasing = TRUE),]

cell_type_list <- unique(marker_list$cluster)
marker_list <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(nrow(temp) <= 200){
    return(temp)
  }else{
    return(temp[1:200,])
  }
}))

#try my add module score
layer_list <- c('CPo','CPi','SP','SZo','SZi','VZ')
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
microarray_Seurat@assays$RNA@data <- as.matrix(log1p(express_matrix))
for (i in cell_type_list) {
  microarray_Seurat <- temp_add_module_score(seu.obj = microarray_Seurat,gene_list = (marker_list$gene)[marker_list$cluster == i],name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,cell_type_list]
colnames(meta_data) <- sub(pattern = ' ',replacement = '_',x = colnames(meta_data),fixed = TRUE)
colnames(meta_data) <- sub(pattern = '-',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/pcw16_microarray_cell_type_marker_temp_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

## pcw 21 ------------------------------------------------------------------
#load data
pcw_21 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_table.rds')
pcw_21_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_sampleinfor_table.rds')

temp <- which(grepl(pattern = 'dorsolateral prefrontal cortex',x = pcw_21_meta_data$structure_name,fixed = TRUE) & pcw_21_meta_data$age == 'pcw21_03')
pcw_21 <- pcw_21[,temp]
pcw_21_meta_data <- pcw_21_meta_data[temp,]
colnames(pcw_21) <- pcw_21_meta_data$region
express_matrix <- pcw_21

#marker gene
marker_list <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')
marker_list <- marker_list[marker_list$gene %in% rownames(express_matrix),]
table(marker_list$cluster)
marker_list$diff <- marker_list$mean1 - marker_list$mean2
marker_list <- marker_list[order(marker_list$diff,decreasing = TRUE),]

cell_type_list <- unique(marker_list$cluster)
marker_list <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(nrow(temp) <= 200){
    return(temp)
  }else{
    return(temp[1:200,])
  }
}))

#try my add module score
layer_list <- c('MZ','CPo','CPi','SP','IZ','SZo','SZi','VZ')
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
microarray_Seurat@assays$RNA@data <- as.matrix(log1p(express_matrix))
for (i in cell_type_list) {
  microarray_Seurat <- temp_add_module_score(seu.obj = microarray_Seurat,gene_list = (marker_list$gene)[marker_list$cluster == i],name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,cell_type_list]
colnames(meta_data) <- sub(pattern = ' ',replacement = '_',x = colnames(meta_data),fixed = TRUE)
colnames(meta_data) <- sub(pattern = '-',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/pcw21_microarray_cell_type_marker_temp_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

## average ------------------------------------------------------------------
#load data
pcw_1516 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_table.rds')
pcw_1516_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_sampleinfor_table.rds')
pcw_21 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_table.rds')
pcw_21_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_sampleinfor_table.rds')

#express_matrix
table(rownames(pcw_1516) == rownames(pcw_21))
express_matrix <- cbind(pcw_1516,pcw_21)
meta_data <- data.frame(region = c(pcw_1516_meta_data$region,pcw_21_meta_data$region),age = c(pcw_1516_meta_data$age,pcw_21_meta_data$age))

#average express matrix
layer_list <- c('MZ','CPo','CPi','SP','IZ','SZo','SZi','VZ')
layer_list %in% meta_data$region

express_matrix <- base::do.call(what = cbind,args = base::lapply(X = layer_list,FUN = function(x){
  temp <- express_matrix[,meta_data$region == x,drop = FALSE]
  temp <- rowMeans(temp)
  return(temp)
}))
colnames(express_matrix) <- layer_list
rownames(express_matrix) <- rownames(pcw_1516)

#marker gene
marker_list <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')
marker_list <- marker_list[marker_list$gene %in% rownames(express_matrix),]
table(marker_list$cluster)
marker_list$diff <- marker_list$mean1 - marker_list$mean2
marker_list <- marker_list[order(marker_list$diff,decreasing = TRUE),]

cell_type_list <- unique(marker_list$cluster)
marker_list <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(nrow(temp) <= 200){
    return(temp)
  }else{
    return(temp[1:200,])
  }
}))

#try my add module score
layer_list <- c('MZ','CPo','CPi','SP','IZ','SZo','SZi','VZ')
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
microarray_Seurat@assays$RNA@data <- as.matrix(log1p(express_matrix))
for (i in cell_type_list) {
  microarray_Seurat <- temp_add_module_score(seu.obj = microarray_Seurat,gene_list = (marker_list$gene)[marker_list$cluster == i],name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,cell_type_list]
colnames(meta_data) <- sub(pattern = ' ',replacement = '_',x = colnames(meta_data),fixed = TRUE)
colnames(meta_data) <- sub(pattern = '-',replacement = '_',x = colnames(meta_data),fixed = TRUE)
meta_data <- meta_data[layer_list,c('Early_RG','Late_RG','tRG','Cycling','IP','Ex_1','Ex_2','Ex_3','Ex_4','Ex_5','Ex_6','Ex_7','Ex_8','Ex_9','Ex_10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]

#plot
pdf(file = './res/step_59_fig_220916/all_average_microarray_cell_type_marker_temp_module_score_heatmap.pdf',width = 8,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE)
dev.off()

# use average and final process ---------------------------------------------
#my function
temp_add_module_score <- function(seu.obj,gene_list,name){
  temp <- seu.obj@assays$RNA@counts[gene_list,]
  temp <- t(base::scale(x = t(temp),center = TRUE,scale = TRUE))
  temp <- colMeans(x = temp)
  seu.obj@meta.data[,name] <- temp
  return(seu.obj)
}

#load data
pcw_1516 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_table.rds')
pcw_1516_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw15pcw16_majorregions_sampleinfor_table.rds')
pcw_21 <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_table.rds')
pcw_21_meta_data <- readRDS(file = './data/public/Transcriptional_landscape_of_the_prenatal_human_brain/microarray_from_liuyt/lmd_matrix_hg_pcw21_majorregions_sampleinfor_table.rds')

#express_matrix
table(rownames(pcw_1516) == rownames(pcw_21))
express_matrix <- cbind(pcw_1516,pcw_21)
meta_data <- data.frame(region = c(pcw_1516_meta_data$region,pcw_21_meta_data$region),age = c(pcw_1516_meta_data$age,pcw_21_meta_data$age))

#average express matrix
layer_list <- c('MZ','CPo','CPi','SP','IZ','SZo','SZi','VZ')
layer_list %in% meta_data$region

express_matrix <- base::do.call(what = cbind,args = base::lapply(X = layer_list,FUN = function(x){
  temp <- express_matrix[,meta_data$region == x,drop = FALSE]
  temp <- rowMeans(temp)
  return(temp)
}))
colnames(express_matrix) <- layer_list
rownames(express_matrix) <- rownames(pcw_1516)

#marker gene
marker_list <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_cell_type_marker.rds')
marker_list <- marker_list[marker_list$gene %in% rownames(express_matrix),]
table(marker_list$cluster)
marker_list$diff <- marker_list$mean1 - marker_list$mean2
marker_list <- marker_list[order(marker_list$diff,decreasing = TRUE),]

cell_type_list <- unique(marker_list$cluster)
marker_list <- base::do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- marker_list[marker_list$cluster == x,]
  if(nrow(temp) <= 200){
    return(temp)
  }else{
    return(temp[1:200,])
  }
}))

#try my add module score
microarray_Seurat <- CreateSeuratObject(counts = express_matrix,project = 'microarray',assay = 'RNA',min.cells = 0,min.features = 0)
microarray_Seurat$region <- colnames(microarray_Seurat)
microarray_Seurat@assays$RNA@data <- as.matrix(log1p(express_matrix))
for (i in cell_type_list) {
  microarray_Seurat <- temp_add_module_score(seu.obj = microarray_Seurat,gene_list = (marker_list$gene)[marker_list$cluster == i],name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- microarray_Seurat@meta.data[,cell_type_list]
meta_data <- meta_data[layer_list[layer_list != 'MZ'],c('Early RG','Late RG','tRG','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','Ex-5','Ex-6','Ex-7','Ex-8','Ex-9','Ex-10','InMGE','InCGE','oIPC','OPC','End','per','Mic')]
colnames(meta_data) <- c('Early RG','Late RG','tRG','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','Ex-5','Ex-6','Ex-7','Ex-8','Ex-9','Ex-10','InMGE','InCGE','oIPC','OPC','End','Per','Mic')

#col seq
cl <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_anno_clustering.rds')
cl <- as.hclust(cl)
cl <- cl$labels[cl$order]
cl[cl == 'per'] <- 'Per'
meta_data <- meta_data[,cl]

#change color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$species
col_fun <- colorRamp2(breaks = c(min(meta_data),max(meta_data)),colors = c('white',temp_col['human']))
col_fun <- colorRamp2(breaks = c(min(meta_data),min(meta_data)+6/10*(max(meta_data)-min(meta_data)),max(meta_data)),colors = c('white',col_fun(1/2*(min(meta_data)+max(meta_data))),temp_col['human']))
pdf(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_anno_layer_enrichment.pdf',width = 10,height = 4)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE,col = col_fun,
        width = unit(8,'inches'),height = unit(8/22*7,'inches'),name = 'enrichment')
dev.off()

