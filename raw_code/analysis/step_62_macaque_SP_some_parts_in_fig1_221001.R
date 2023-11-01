#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque SP some parts in fig1                                   ##
## Data: 2022.10.01                                                                ##
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
library(scales)
library(DESeq2)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#ggplot theme from liuyt
# theme blank: no axis, no legend, no title, no frame
# for plots saved as png
theme_blank <-
  NoAxes() +
  theme(
    plot.title = element_blank(),
    panel.border = element_blank(),
    aspect.ratio = 1,
    legend.position = 'none',
    plot.margin = margin(0,0,0,0, "cm")
  ) 


# theme with umap axis arrow and legend
# for plots not saved as png, just for recording information
theme_arrow <- 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5, size = 10),
        panel.grid = element_blank()) +
  theme(text=element_text(size=10,  family="sans"), line = element_line(size = 0.8)) + 
  guides(color = guide_legend(override.aes = list(size = 3, shape = 15)))

# unify three species hclust tree -----------------------------------------
## human -------------------------------------------------------------------
cl.human <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_anno_clustering.rds')
plot(cl.human)
cl.human <- as.hclust(cl.human)
cl.human <- my_order_hclust(cl.obj = cl.human,
                            order.list = c('End','per','Mic','OPC','tRG','oIPC','Late RG','Early RG','Cycling','IP',
                                           'InMGE','InCGE','Ex-1','Ex-5','Ex-4','Ex-2','Ex-3','Ex-6','Ex-7','Ex-8','Ex-9','Ex-10'))
nodePar <- list(pch = c(NA,15),cex = 2)

plot(cl.human,nodePar = nodePar)

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
cl <- cl.human
cl <- as.hclust(cl)
cl <- cl$labels[cl$order]
cl[cl == 'per'] <- 'Per'
meta_data <- meta_data[,cl]

#change color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
temp_col <- col_param$species
col_fun <- colorRamp2(breaks = c(min(meta_data),max(meta_data)),colors = c('white',temp_col['human']))
col_fun <- colorRamp2(breaks = c(min(meta_data),min(meta_data)+6/10*(max(meta_data)-min(meta_data)),max(meta_data)),colors = c('white',col_fun(1/2*(min(meta_data)+max(meta_data))),temp_col['human']))

#column anno
temp_col <- col_param$celltype.human
col_anno <- HeatmapAnnotation(celltype = colnames(meta_data),
                              col = list(celltype = temp_col[colnames(meta_data)]),
                              which = 'column',gp = gpar(color = 'black'),show_legend = FALSE,show_annotation_name = FALSE)

#cluster
cl <- cl.human
cl <- as.hclust(cl)
temp <- cl$labels
temp[temp == 'per'] <- 'Per'
cl$labels <- temp
plot(cl)
meta_data <- meta_data[,cl$labels]

#plot
pdf(file = './res/step_62_fig_221001/Greenleaf_RNA_Seurat_own_anno_layer_enrichment.pdf',width = 10,height = 5)
Heatmap(matrix = meta_data,cluster_columns = cl,cluster_rows = FALSE,col = col_fun,
        width = unit(8,'inches'),height = unit(8/22*7,'inches'),name = 'enrichment',
        top_annotation = col_anno)
dev.off()

## macaque -----------------------------------------------------------------
cl.macaque <- readRDS(file = './res/step_60_fig_220920/macaque_multiome_Seurat_cell_type_clustering.rds')
plot(cl.macaque)
cl.macaque <- as.hclust(cl.macaque)
cl.macaque <- my_order_hclust(cl.obj = cl.macaque,
                              order.list = c('End','Per','VLMC','Mic','OPC','RG-2','RG-1','Cycling','InMGE','InCGE','IP','Ex-1','Ex-2','Ex-3','Ex-4'))
nodePar <- list(pch = c(NA,15),cex = 2)

plot(cl.macaque,nodePar = nodePar)

#layer enrichment
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')
colnames(macaque_express_matrix) <- sub(pattern = '-',replacement = '_',colnames(macaque_express_matrix))

#merge
meta_data <- data.frame(sample = colnames(macaque_express_matrix),species = 'macaque')
meta_data$layer <- unlist(base::lapply(X = meta_data$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][2])
}))

meta_data$donor <- unlist(base::lapply(X = meta_data$sample,FUN = function(x){
  temp <- strsplit(x = x,split = '_',fixed = TRUE)
  return(temp[[1]][1])
}))

rownames(meta_data) <- meta_data$sample

#merge
layer_list <- c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ')
express_matrix <- base::do.call(what = cbind,args = base::lapply(X = layer_list,FUN = function(x){
  temp <- meta_data[meta_data$layer == x,]
  temp <- rownames(temp)
  temp <- macaque_express_matrix[,temp]
  temp <- rowSums(temp)
  return(temp)
}))
colnames(express_matrix) <- layer_list
rownames(express_matrix) <- rownames(macaque_express_matrix)

#convert id
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_name','gene_id')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]
macaque_annotation <- cbind(macaque_annotation,macaque_annotation)

express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = express_matrix,anno = macaque_annotation,filter_anno = FALSE,future.globals.maxSize = 200*(1024^3),workers = 6)

#normalize
meta_data <- data.frame(sample = colnames(express_matrix),layer = colnames(express_matrix))
temp <- DESeqDataSetFromMatrix(countData = express_matrix,colData = meta_data,design = ~ layer)
temp <- rlog(object = temp,blind = TRUE)
express_matrix <- assay(temp)

#marker list
marker_list <- readRDS(file = './res/step_57_fig_220913/cell_type_specific_marker_gene_lfc_0.3_pct_0.2_overlap_0.75.rds')
marker_list <- base::do.call(what = rbind,args = base::lapply(X = names(marker_list),FUN = function(x){
  temp <- marker_list[[x]]
  other_cell_type_list <- unique(macaque_multiome_Seurat$cell_type)
  other_cell_type_list <- other_cell_type_list[other_cell_type_list != x]
  temp <- my_seurat_marker_wilcox_test(seu.obj = macaque_multiome_Seurat,assay = 'RNA',ident.1 = x,ident.2 = other_cell_type_list,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0,pct.1 = 0,features = temp,only.pos = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
  temp$gene <- rownames(temp)
  temp$cluster <- x
  temp <- temp[temp$gene %in% rownames(express_matrix),]
  temp$diff <- temp$pct.1-temp$pct.2
  temp <- temp[order(temp$diff,decreasing = TRUE),]
  if(nrow(temp) <= 200){
    temp <- temp
  }else{
    temp <- temp[1:200,]
  }
  print(paste(x,'done!',sep = ' '))
  gc()
  return(temp)
}))
table(marker_list$cluster)
marker_list[1:10,]

#my function
temp_add_module_score <- function(seu.obj,gene_list,name){
  temp <- seu.obj@assays$RNA@counts[gene_list,]
  temp <- t(base::scale(x = t(temp),center = TRUE,scale = TRUE))
  temp <- colMeans(x = temp,na.rm = TRUE)
  seu.obj@meta.data[,name] <- temp
  return(seu.obj)
}

#try my add module score
cell_type_list <- readRDS(file = './res/step_57_fig_220913/cell_type_specific_marker_gene_lfc_0.3_pct_0.2_overlap_0.75.rds')
cell_type_list <- names(cell_type_list)
express_matrix <- CreateSeuratObject(counts = express_matrix,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)
express_matrix$region <- colnames(express_matrix)
for (i in cell_type_list) {
  express_matrix <- temp_add_module_score(seu.obj = express_matrix,gene_list = (marker_list$gene)[marker_list$cluster == i],name = i)
  print(paste(i,'done!',sep = ' '))
}

#get module score
meta_data <- express_matrix@meta.data[,cell_type_list]
meta_data <- meta_data[layer_list[layer_list != 'MZ'],cell_type_list]

#col seq
cl <- cl.macaque
cl <- as.hclust(cl)
meta_data <- meta_data[,cl$labels]

#change color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
temp_col <- col_param$species
col_fun <- colorRamp2(breaks = c(min(meta_data),max(meta_data)),colors = c('white',temp_col['macaque']))
col_fun <- colorRamp2(breaks = c(min(meta_data),min(meta_data)+6/10*(max(meta_data)-min(meta_data)),max(meta_data)),colors = c('white',col_fun(1/2*(min(meta_data)+max(meta_data))),temp_col['macaque']))

#top anno
temp_col <- col_param$celltype
col_anno <- HeatmapAnnotation(celltype = colnames(meta_data),
                              col = list(celltype = temp_col[colnames(meta_data)]),
                              which = 'column',gp = gpar(color = 'black'),show_legend = FALSE,show_annotation_name = FALSE)


#plot
pdf(file = './res/step_62_fig_221001/macaque_RNA_Seurat_cell_type_layer_enrichment.pdf',width = 10,height = 6)
Heatmap(matrix = meta_data,cluster_columns = cl,cluster_rows = FALSE,col = col_fun,
        width = unit(8,'inches'),height = unit(8/15*7,'inches'),name = 'enrichment',
        top_annotation = col_anno)
dev.off()

## mouse -------------------------------------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/mouse_multiome_Seurat_220922.rds')

#aggregate
express_matrix <- Seurat::AggregateExpression(object = mouse_multiome_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'data',verbose = TRUE)
express_matrix <- express_matrix$RNA

#normalize 
col_name <- colnames(express_matrix)
row_name <- rownames(express_matrix)
express_matrix <- base::do.call(what = cbind,args = base::lapply(X = col_name,FUN = function(x){
  temp <- express_matrix[,x]
  temp <- (temp/sum(temp))*1000000
  return(temp)
}))
colnames(express_matrix) <- col_name
rownames(express_matrix) <- row_name
express_matrix <- log1p(express_matrix)

#express_matrix
express_matrix <- express_matrix[VariableFeatures(mouse_multiome_Seurat),]
dist_matrix <- dist(t(express_matrix))

#clustering
cl.mouse <- hclust(d = dist_matrix,method = 'ward.D')
plot(cl.mouse)
cl.mouse <- my_order_hclust(cl.obj = cl.mouse,
                            order.list = c('End','Per','VLMC','Mic','RG','Cyc','Cyc-IP','IP','In','MGN','CPN','SCPN','CThPN'))
nodePar <- list(pch = c(NA,15),cex = 2)

plot(cl.mouse,nodePar = nodePar)

#cl
cl <- cl.mouse
cl <- as.hclust(cl)
plot(cl)

#plot
meta_data <- matrix(data = 1,nrow = 8,ncol = length(cl$labels))
rownames(meta_data) <- c('MZ','CPo','CPi','SP','IZ','OSVZ','ISVZ','VZ')
colnames(meta_data) <- cl$labels

#top_anno
temp_col <- col_param$celltype.mouse
col_anno <- HeatmapAnnotation(celltype = colnames(meta_data),
                              col = list(celltype = temp_col[colnames(meta_data)]),
                              which = 'column',gp = gpar(color = 'black'),show_legend = FALSE,show_annotation_name = FALSE)

#plot
pdf(file = './res/step_62_fig_221001/mouse_RNA_Seurat_cell_type_layer_enrichment.pdf',width = 10,height = 7)
Heatmap(matrix = meta_data,cluster_columns = cl,cluster_rows = FALSE,
        width = unit(8,'inches'),height = unit(8/13*8,'inches'),name = 'enrichment',
        top_annotation = col_anno)
dev.off()

# plot legend -------------------------------------------------------------
cl.human <- readRDS(file = './res/step_59_fig_220916/Greenleaf_RNA_Seurat_own_anno_clustering.rds')
plot(cl.human)
cl.human <- as.hclust(cl.human)
cl.human <- my_order_hclust(cl.obj = cl.human,
                            order.list = c('End','per','Mic','OPC','tRG','oIPC','Late RG','Early RG','Cycling','IP',
                                           'InMGE','InCGE','Ex-1','Ex-5','Ex-4','Ex-2','Ex-3','Ex-6','Ex-7','Ex-8','Ex-9','Ex-10'))
nodePar <- list(pch = c(NA,15),cex = 2)

plot(cl.human,nodePar = nodePar)

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
cl <- cl.human
cl <- as.hclust(cl)
cl <- cl$labels[cl$order]
cl[cl == 'per'] <- 'Per'
meta_data <- meta_data[,cl]

#change color
col_fun <- colorRamp2(breaks = c(min(meta_data),max(meta_data)),colors = c('white','black'))
col_fun <- colorRamp2(breaks = c(min(meta_data),min(meta_data)+6/10*(max(meta_data)-min(meta_data)),max(meta_data)),colors = c('white',col_fun(1/2*(min(meta_data)+max(meta_data))),'black'))

#column anno
temp_col <- col_param$celltype.human
col_anno <- HeatmapAnnotation(celltype = colnames(meta_data),
                              col = list(celltype = temp_col[colnames(meta_data)]),
                              which = 'column',gp = gpar(color = 'black'),show_legend = FALSE,show_annotation_name = FALSE)

#cluster
cl <- cl.human
cl <- as.hclust(cl)
temp <- cl$labels
temp[temp == 'per'] <- 'Per'
cl$labels <- temp
plot(cl)
meta_data <- meta_data[,cl$labels]

#plot
pdf(file = './res/step_62_fig_221001/heatmap_legend.pdf',width = 10,height = 5)
Heatmap(matrix = meta_data,cluster_columns = cl,cluster_rows = FALSE,col = col_fun,
        width = unit(8,'inches'),height = unit(8/22*7,'inches'),name = 'enrichment',
        top_annotation = col_anno)
dev.off()