#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: using monocle3 to investigate Excite Neuron                     ##
## Data: 2021.09.07                                                                ##
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
library(monocle3)
library(factoextra)
library(circlize)
library(ggsci)
library(SeuratWrappers)
library(magrittr)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

# load data ---------------------------------------------------------------
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'

#modify the annotation
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('oRG','vRG'),"cell_type"] <- 'RG'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('30','9','1','4','22','2','12'),'cell_type'] <- 'Ex-1'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('14','8'),'cell_type'] <- 'Ex-2'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('16','15','24','6','18'),'cell_type'] <- 'Ex-3'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$SCT_snn_res.1.6 %in% c('19','20'),'cell_type'] <- 'Ex-4'

macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type == 'oRG',"sub_cell_type"] <- 'RG-1'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type == 'vRG-1',"sub_cell_type"] <- 'RG-2'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$sub_cell_type == 'vRG-2',"sub_cell_type"] <- 'RG-3'

macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-1',"sub_cell_type"] <- 'Ex-1'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-2',"sub_cell_type"] <- 'Ex-2'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-3',"sub_cell_type"] <- 'Ex-3'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'Ex-4',"sub_cell_type"] <- 'Ex-4'

DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)

#load and modify clolor paramater
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
color_paramater <- color_paramater[c(1:22,25,26),]
color_paramater[which(color_paramater$id == 'oRG'),'id'] <- 'RG'
color_paramater[which(color_paramater$id == 'vRG'),'id'] <- 'Ex-2'
color_paramater[which(color_paramater$id == 'Ex'),'id'] <- 'Ex-1'
color_paramater[which(color_paramater$id == 'Ex-U'),'id'] <- 'Ex-3'
color_paramater[which(color_paramater$id == 'Ex-SP'),'id'] <- 'Ex-4'

temp_col <- data.frame(id = c('RG-1','RG-2','RG-3'),col = colorRampPalette(c("#00A087FF","#3C5488FF"))(5)[2:4])
color_paramater <- rbind(color_paramater,temp_col)

temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,cols = temp_col)

# marker gene express -----------------------------------------------------
geneSetAveragePlot(genes = c('NRP1','NEUROD2','TUBB3','BCL11B','FEZF2','FOXP2','NR4A2','CRYM'),object = macaque_RNA_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,lims = NULL,print = FALSE)
#dead end

# layer marker express ----------------------------------------------------
macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]
layer_marker <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')

#layer marker
my_temp_function <- function(layer_marker,layer){
  marker <- layer_marker[[layer]]
  marker <- unique(marker)
  print('marker in annotation list')
  print(table(marker %in% macaque_annotation$gene_id))
  marker <- marker[marker %in% macaque_annotation$gene_id]
  print('marker trans id to name')
  marker <- macaque_annotation$gene_name[base::match(marker,table = macaque_annotation$gene_id)]
  print('marker in express matrix')
  print(table(marker %in% rownames(macaque_RNA_seurat)))
  marker <- marker[marker %in% rownames(macaque_RNA_seurat)]
  print('duplicated markers')
  print(table(duplicated(marker)))
  marker <- unique(marker)
  return(marker)
}

names(layer_marker)

for (i in c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")){
  temp_marker <- my_temp_function(layer_marker = layer_marker,layer = i)
  macaque_RNA_seurat <- My_add_module_score(seu.obj = macaque_RNA_seurat,assay = 'RNA',features = temp_marker,meta_var = i,scale = FALSE,center = FALSE)
}

macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG')]
meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")]

#annotate
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

cell_type_annotation <- HeatmapAnnotation(cell_type=macaque_RNA_seurat$cell_type,
                                          col = list(cell_type = temp_col[names(table(macaque_RNA_seurat$cell_type))]),
                                          show_annotation_name = TRUE,which = 'row',name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

col_fun <- colorRamp2(c(0,0.5,1),c("#440154FF","#21908CFF","#FDE725FF"))

pdf(file = './res/step_8_fig_210909/macaque_RNA_seurat_layer_marker_module_score_heatmap.pdf',width = 10,height = 8)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(macaque_RNA_seurat$cell_type,
                           levels = c('InPSB','InMGE','InCGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(6,units = 'inches'),width = unit(4,units = 'inches'),
        name = 'expression',border = TRUE,col = col_fun)
dev.off()
# monocle3 ----------------------------------------------------------------
## subset Ex ---------------------------------------------------------------
macaque_Ex_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('IP','Ex-1','Ex-2','Ex-3','Ex-4')]
macaque_Ex_seurat <- macaque_Ex_seurat[,macaque_Ex_seurat@reductions$umap@cell.embeddings[,'UMAP_1'] > -6 & macaque_Ex_seurat@reductions$umap@cell.embeddings[,'UMAP_2'] < 10]
DimPlot(macaque_Ex_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

## create cds object -------------------------------------------------------
cds <- SeuratWrappers::as.cell_data_set(macaque_Ex_seurat)
cds <- cluster_cells(cds)
plot_cells(cds,show_trajectory_graph = FALSE,color_cells_by = 'cell_type')
plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
cds <- learn_graph(cds,learn_graph_control=list(minimal_branch_len=10))

pdf(file = './res/step_8_fig_210909/macaque_Ex_seurat_cell_type_trajectory.pdf',width = 6,height = 6)
plot_cells(cds,show_trajectory_graph = TRUE,color_cells_by = 'cell_type',
           label_branch_points = FALSE,label_groups_by_cluster = FALSE,
           label_leaves = FALSE,group_label_size = 0) + 
  theme_classic() + 
  scale_color_manual(values = temp_col[c('IP','Ex-1','Ex-2','Ex-3','Ex-4')]) + 
  theme(aspect.ratio = 1)
dev.off()

#set root
cds <- order_cells(cds)
pdf(file = './res/step_8_fig_210909/macaque_Ex_seurat_pseudotime_trajectory.pdf',width = 6,height = 6)
plot_cells(cds,show_trajectory_graph = TRUE,color_cells_by = 'pseudotime',
           label_branch_points = FALSE,label_groups_by_cluster = FALSE,
           label_leaves = FALSE,group_label_size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 1)
dev.off()

p1 <- plot_cells(cds,show_trajectory_graph = TRUE,color_cells_by = 'cell_type',
                 label_branch_points = FALSE,label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,group_label_size = 0) + 
  theme_classic() + 
  scale_color_manual(values = temp_col[c('IP','Ex-1','Ex-2','Ex-3','Ex-4')]) + 
  theme(aspect.ratio = 1)

p2 <- plot_cells(cds,show_trajectory_graph = TRUE,color_cells_by = 'pseudotime',
                 label_branch_points = FALSE,label_groups_by_cluster = FALSE,
                 label_leaves = FALSE,group_label_size = 0) + 
  theme_classic() + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_8_fig_210909/macaque_Ex_seurat_combined_trajectory.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#save data
saveRDS(cds,file = './res/step_8_fig_210909/macaque_Ex_cds_210912.rds')
