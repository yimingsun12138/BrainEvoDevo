#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque multiome Seurat cell type layer distribution            ##
## Data: 2022.09.19                                                                ##
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
library(DESeq2)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
macaque_multiome_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_multiome_Seurat_220802.rds')
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')
colnames(macaque_express_matrix) <- sub(pattern = '-',replacement = '_',colnames(macaque_express_matrix))

# cluster macaque multiome Seurat -----------------------------------------
express_matrix <- Seurat::AggregateExpression(object = macaque_multiome_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'data',verbose = TRUE)
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
express_matrix <- express_matrix[VariableFeatures(macaque_multiome_Seurat),]

#dist
dist_matrix <- dist(t(express_matrix))
cl <- hclust(d = dist_matrix,method = 'ward.D')
cl <- as.dendrogram(cl)
cl <- reorder(x = cl,c(2^7,2^1,2^11,2^12,2^13,2^14,2^9,2^8,2^10,2^3,2^4,2^0,2^6,2^5,2^2))

pdf(file = './res/step_60_fig_220920/macaque_multiome_Seurat_cell_type_clustering.pdf',width = 8,height = 3)
plot(cl)
dev.off()

saveRDS(object = cl,file = './res/step_60_fig_220920/macaque_multiome_Seurat_cell_type_clustering.rds')

# merge bulk express ------------------------------------------------------
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

# convert gene id ---------------------------------------------------------
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_name','gene_id')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]
macaque_annotation <- cbind(macaque_annotation,macaque_annotation)

express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = express_matrix,anno = macaque_annotation,filter_anno = FALSE,future.globals.maxSize = 200*(1024^3),workers = 6)

# normalize ---------------------------------------------------------------
meta_data <- data.frame(sample = colnames(express_matrix),layer = colnames(express_matrix))
temp <- DESeqDataSetFromMatrix(countData = express_matrix,colData = meta_data,design = ~ layer)
temp <- rlog(object = temp,blind = TRUE)
express_matrix <- assay(temp)

# marker list -------------------------------------------------------------
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

# add module score --------------------------------------------------------
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
cl <- readRDS(file = './res/step_60_fig_220920/macaque_multiome_Seurat_cell_type_clustering.rds')
cl <- as.hclust(cl)
cl <- cl$labels[cl$order]
meta_data <- meta_data[,cl]

#change color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$species
col_fun <- colorRamp2(breaks = c(min(meta_data),max(meta_data)),colors = c('white',temp_col['macaque']))
col_fun <- colorRamp2(breaks = c(min(meta_data),min(meta_data)+6/10*(max(meta_data)-min(meta_data)),max(meta_data)),colors = c('white',col_fun(1/2*(min(meta_data)+max(meta_data))),temp_col['macaque']))
pdf(file = './res/step_60_fig_220920/macaque_RNA_Seurat_cell_type_layer_enrichment.pdf',width = 10,height = 6)
Heatmap(matrix = meta_data,cluster_columns = FALSE,cluster_rows = FALSE,col = col_fun,
        width = unit(8,'inches'),height = unit(8/15*7,'inches'),name = 'enrichment')
dev.off()