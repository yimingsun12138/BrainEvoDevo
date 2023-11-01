#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: using monocle2 to investigate progenitor cell                   ##
## Data: 2021.06.29                                                                ##
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
library(monocle)
library(factoextra)
library(circlize)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
#subset
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$seurat_clusters %in% c('13','26','23','21','17','30','9','1','25')]
table(macaque_RNA_seurat$cell_type)
#create celldataset object
macaque_RNA_express <- macaque_RNA_seurat@assays$RNA@counts
macaque_RNA_metacol <- macaque_RNA_seurat@meta.data
macaque_RNA_metarow <- data.frame(gene_short_name=rownames(macaque_RNA_express),gene_id=rownames(macaque_RNA_express))
rownames(macaque_RNA_metarow) <- macaque_RNA_metarow$gene_short_name

pd <- new("AnnotatedDataFrame", data = macaque_RNA_metacol)
fd <- new("AnnotatedDataFrame", data = macaque_RNA_metarow)
macaque_RNA_monocle <- newCellDataSet(macaque_RNA_express,phenoData = pd,featureData = fd,lowerDetectionLimit = 0.1,expressionFamily = negbinomial.size())
rm(list = c('fd','pd','macaque_RNA_express','macaque_RNA_metacol','macaque_RNA_metarow'))

macaque_RNA_monocle <- estimateSizeFactors(macaque_RNA_monocle)
macaque_RNA_monocle <- estimateDispersions(macaque_RNA_monocle)
#QC
macaque_RNA_monocle <- detectGenes(macaque_RNA_monocle,min_expr = 3)
head(fData(macaque_RNA_monocle))
expressed_genes <- row.names(subset(fData(macaque_RNA_monocle),num_cells_expressed >= 50))
#trajectory
#find variable gene
diff_test_res <- differentialGeneTest(macaque_RNA_monocle[expressed_genes,],fullModelFormulaStr = "~sub_cell_type+donor+nCount_RNA",reducedModelFormulaStr = "~donor+nCount_RNA")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
macaque_RNA_monocle <- setOrderingFilter(macaque_RNA_monocle,ordering_genes = ordering_genes)
plot_ordering_genes(macaque_RNA_monocle)
#reduce dimension
macaque_RNA_monocle <- reduceDimension(macaque_RNA_monocle, max_components = 2, method = 'DDRTree',verbose = TRUE,norm_method = 'log')
#visulize
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id
macaque_RNA_monocle <- orderCells(macaque_RNA_monocle)

pdf(file = './res/step_4_fig_210629/progenitor_cell_trajectory.pdf',width = 4,height = 4)
plot_cell_trajectory(macaque_RNA_monocle, color_by = "sub_cell_type",show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_manual(values = temp_col[c('oRG','vRG-1','vRG-2','Cyc-S','Cyc-G2M','IP','OPC','Ex')]) + 
  theme_classic() + 
  theme(legend.position = 'none',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  guides(colour = guide_legend(override.aes = list(shape = 15,size = 5))) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = 'Cell type')
dev.off()

pdf(file = './res/step_4_fig_210629/progenitor_cell_trajectory_legend.pdf',width = 5,height = 5)
plot_cell_trajectory(macaque_RNA_monocle, color_by = "sub_cell_type",show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_manual(values = temp_col[c('oRG','vRG-1','vRG-2','Cyc-S','Cyc-G2M','IP','OPC','Ex')]) + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  guides(colour = guide_legend(override.aes = list(shape = 15,size = 5))) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = 'Cell type')
dev.off()

temp_col <- c("#4691B6",'grey','grey','grey','grey','grey','grey','grey')
names(temp_col) <- c('vRG-2','oRG','vRG-1','Cyc-S','Cyc-G2M','IP','OPC','Ex')
pdf(file = './res/step_4_fig_210629/progenitor_cell_trajectory_vRG_2.pdf',width = 4,height = 4)
plot_cell_trajectory(macaque_RNA_monocle, color_by = "sub_cell_type",show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_manual(values = temp_col[c('oRG','vRG-1','vRG-2','Cyc-S','Cyc-G2M','IP','OPC','Ex')]) + 
  theme_classic() + 
  theme(legend.position = 'none',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  guides(colour = guide_legend(override.aes = list(shape = 15,size = 5))) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = 'Cell type',title = 'vRG-2')
dev.off()



temp_col <- c("#49A6C5",'grey','grey','grey','grey','grey','grey','grey')
names(temp_col) <- c('vRG-1','oRG','vRG-2','Cyc-S','Cyc-G2M','IP','OPC','Ex')
pdf(file = './res/step_4_fig_210629/progenitor_cell_trajectory_vRG_1.pdf',width = 4,height = 4)
plot_cell_trajectory(macaque_RNA_monocle, color_by = "sub_cell_type",show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_manual(values = temp_col[c('oRG','vRG-1','vRG-2','Cyc-S','Cyc-G2M','IP','OPC','Ex')]) + 
  theme_classic() + 
  theme(legend.position = 'none',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  guides(colour = guide_legend(override.aes = list(shape = 15,size = 5))) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = 'Cell type',title = 'vRG-1')
dev.off()

#set root
plot_cell_trajectory(macaque_RNA_monocle,color_by = 'State')
macaque_RNA_monocle <- orderCells(macaque_RNA_monocle,root_state = '1')

pdf(file = './res/step_4_fig_210629/progenitor_cell_trajectory_pseudotime.pdf',width = 4,height = 4)
plot_cell_trajectory(macaque_RNA_monocle,color_by = 'Pseudotime',show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "#440154FF",high = "#FDE725FF") + 
  theme_classic() + 
  theme(legend.position = 'none',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = 'Pseudotime')
dev.off()

pdf(file = './res/step_4_fig_210629/progenitor_cell_trajectory_pseudotime_legend.pdf',width = 5,height = 5)
plot_cell_trajectory(macaque_RNA_monocle,color_by = 'Pseudotime',show_branch_points = FALSE,cell_size = 0.5) + 
  theme_classic() + 
  scale_colour_gradient(low = "#440154FF",high = "#FDE725FF") + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = 'Pseudotime')
dev.off()

#express change through pseudotime
BEAM_res <- BEAM(macaque_RNA_monocle,branch_point = 1,cores = 1,verbose = TRUE)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
saveRDS(BEAM_res,file = './res/step_4_fig_210629/progenitor_cell_monocle_BEAM_res.rds')
saveRDS(macaque_RNA_monocle,file = './res/step_4_fig_210629/progenitor_cell_monocle.rds')

plot_genes_branched_heatmap(macaque_RNA_monocle[row.names(subset(BEAM_res, qval < 1e-10 & num_cells_expressed > 100)),], 
                            branch_point = 1, num_clusters = 4, cores = 1, 
                            use_gene_short_name = T, show_rownames = T)

# #find marker gene
# Idents(macaque_RNA_seurat) <- 'sub_cell_type'
# DefaultAssay(macaque_RNA_seurat) <- 'RNA'
# progenitor_markers <- FindAllMarkers(macaque_RNA_seurat,assay = 'RNA',slot = 'data',
#                                      test.use = 'bimod',only.pos = TRUE,
#                                      logfc.threshold = 0.25,min.pct = 0.25)
# progenitor_markers <- progenitor_markers[progenitor_markers$p_val_adj < 0.01,]
# progenitor_markers <- unique(progenitor_markers$gene)
# progenitor_markers <- dplyr::intersect(row.names(subset(BEAM_res,qval < 1e-4)),progenitor_markers)

###############################
##reload data to plot heatmap##
###############################
macaque_RNA_monocle <- readRDS(file = './res/step_4_fig_210629/progenitor_cell_monocle.rds')
BEAM_res <- readRDS(file = './res/step_4_fig_210629/progenitor_cell_monocle_BEAM_res.rds')
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_monocle)]

progenitor_markers <- row.names(subset(BEAM_res,qval < 1e-20 & num_cells_expressed > 50))

#generate matrix
heatmap_matrix <- macaque_RNA_seurat@assays$RNA@data
heatmap_matrix <- heatmap_matrix[progenitor_markers,]
heatmap_matrix <- t(scale(t(heatmap_matrix),center = TRUE))
heatmap_matrix <- data.frame(heatmap_matrix)
colnames(heatmap_matrix) <- colnames(macaque_RNA_seurat)

#cell_order
cell_order <- data.frame(cell_id=colnames(macaque_RNA_seurat),cell_type=macaque_RNA_seurat$sub_cell_type)
rownames(cell_order) <- cell_order$cell_id
temp <- 1:8
names(temp) <- c('vRG-1','vRG-2','oRG','Cyc-S','Cyc-G2M','IP','Ex','OPC')
cell_order$cell_order <- as.numeric(temp[cell_order$cell_type])
cell_order <- cell_order[order(cell_order$cell_order,decreasing = FALSE),]

temp <- macaque_RNA_monocle$Pseudotime
names(temp) <- colnames(macaque_RNA_monocle)
temp <- temp[cell_order$cell_id]
cell_order$Pseudotime <- temp
cell_order <- cell_order[order(cell_order$Pseudotime,decreasing = FALSE),]

cell_order$state <- pData(macaque_RNA_monocle)[rownames(cell_order),'State']
cell_order <- cell_order[order(cell_order$state,decreasing = FALSE),]

heatmap_matrix <- heatmap_matrix[,cell_order$cell_id]
#feature cluster
# fviz_nbclust(heatmap_matrix,kmeans,method="wss",k.max = 10)
# k_means_cluster <- kmeans(heatmap_matrix,7)
# saveRDS(k_means_cluster,file = './res/step_4_fig_210629/progenitor_gene_kmeans.rds')
k_means_cluster <- readRDS(file = './res/step_4_fig_210629/progenitor_gene_kmeans.rds')
gene_order <- data.frame(gene=rownames(heatmap_matrix),cluster=k_means_cluster$cluster[rownames(heatmap_matrix)])
rownames(gene_order) <- gene_order$gene
temp <- letters[1:7]
names(temp) <- as.character(c(1,2,5,3,7,6,4))
gene_order$gene_order <- temp[as.character(gene_order$cluster)]
gene_order <- gene_order[order(gene_order$gene_order,decreasing = FALSE),]
heatmap_matrix <- heatmap_matrix[gene_order$gene,]
#create annotation
col_fun <- colorRamp2(c(0,25),c("#440154FF","#FDE725FF"))
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

cell_anno_heatmap <- HeatmapAnnotation(
  cell_type = as.character(cell_order$cell_type),
  Pseudotime = as.numeric(cell_order$Pseudotime),
  col = list(cell_type = temp_col[c('oRG','vRG-1','vRG-2','Cyc-S','Cyc-G2M','IP','OPC','Ex')],
             Pseudotime = col_fun)
)

gene_col <- ggsci::pal_igv()(7)
names(gene_col) <- unique(gene_order$cluster)
gene_anno_heatmap <- rowAnnotation(
  group = as.character(gene_order[rownames(heatmap_matrix),"cluster"]),
  col = list(group = gene_col)
)
#plot complexheatmap
heatmap_col <- colorRamp2(c(-4,0,4),c('#3B4992FF','white','#EE0000FF'))
Heatmap(heatmap_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,show_row_names = FALSE,
        height = unit(7,units = 'inches'),width = unit(5,units = 'inches'),top_annotation = cell_anno_heatmap,
        left_annotation = gene_anno_heatmap,
        row_split = factor(gene_order[rownames(heatmap_matrix),"cluster"],levels = as.character(c(1,2,5,3,7,6,4))),
        column_split = factor(cell_order[colnames(heatmap_matrix),"state"],levels = 1:3),col = heatmap_col,name = 'z score')

#plot TF
#TF list from http://humantfs.ccbr.utoronto.ca/download.php
TF <- c('GFAP','HOPX','TNC','EGFR','PCDH15')
TF %in% rownames(heatmap_matrix)

classic_TF <- rowAnnotation(TF = anno_mark(at = which(gene_order$gene %in% TF),
                                           labels = gene_order$gene[which(gene_order$gene %in% TF)]))
pdf(file = './res/step_4_fig_210629/progenitor_cell_gene_cascade.pdf',width = 9,height = 9)
Heatmap(heatmap_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,show_row_names = FALSE,
        height = unit(7,units = 'inches'),width = unit(5,units = 'inches'),top_annotation = cell_anno_heatmap,
        left_annotation = gene_anno_heatmap,right_annotation = classic_TF,
        row_split = factor(gene_order[rownames(heatmap_matrix),"cluster"],levels = as.character(c(1,2,5,3,7,6,4))),
        column_split = factor(cell_order[colnames(heatmap_matrix),"state"],levels = 1:3),col = heatmap_col,name = 'z score')
dev.off()


##################################
##reload data to investigate vRG##
##################################
macaque_RNA_monocle <- readRDS(file = './res/step_4_fig_210629/progenitor_cell_monocle.rds')
BEAM_res <- readRDS(file = './res/step_4_fig_210629/progenitor_cell_monocle_BEAM_res.rds')
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_monocle)]

#density plot
temp <- pData(macaque_RNA_monocle)
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

pdf(file = './res/step_4_fig_210629/progenitor_cell_pseudotime_density.pdf',width = 6,height = 5)
ggplot(temp,aes(x=Pseudotime,colour=sub_cell_type)) + geom_density() + 
  theme_classic() + 
  scale_colour_manual(values = temp_col[c('oRG','vRG-1','vRG-2','Cyc-S','Cyc-G2M','IP','OPC','Ex')]) + 
  theme(aspect.ratio = 1)
dev.off()

ggplot(temp[!(temp$cell_type %in% c('OPC','Ex')),],aes(x=Pseudotime,colour=sub_cell_type)) + geom_density() + 
  theme_classic() + 
  scale_colour_manual(values = temp_col[c('oRG','vRG-1','vRG-2','Cyc-S','Cyc-G2M','IP')]) + 
  theme(aspect.ratio = 1)

#featureplot
macaque_RNA_seurat@meta.data[colnames(macaque_RNA_monocle),'State'] <- macaque_RNA_monocle$State
macaque_RNA_seurat@meta.data[colnames(macaque_RNA_monocle),'Pseudotime'] <- macaque_RNA_monocle$Pseudotime

temp <- macaque_RNA_seurat[,macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_1'] > 0 & macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_2'] > -5]
p1 <- FeaturePlot(temp,features = 'Pseudotime',cols = c("#440154FF","#FDE725FF")) + 
  theme(aspect.ratio = 1)
p2 <- DimPlot(temp,group.by = 'State',label = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_4_fig_210629/progenitor_cell_state_pseudotime_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#gene express trajectory
macaque_RNA_monocle$PAX6 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['PAX6',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$TNC <- as.numeric(macaque_RNA_seurat@assays$RNA@data['TNC',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$VIM <- as.numeric(macaque_RNA_seurat@assays$RNA@data['VIM',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$EGFR <- as.numeric(macaque_RNA_seurat@assays$RNA@data['EGFR',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$OLIG2 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['OLIG2',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$SOX10 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['SOX10',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$PCDH15 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['PCDH15',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$EOMES <- as.numeric(macaque_RNA_seurat@assays$RNA@data['EOMES',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$BCL11B <- as.numeric(macaque_RNA_seurat@assays$RNA@data['BCL11B',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$NEUROD2 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['NEUROD2',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$NEUROD6 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['NEUROD6',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$SATB2 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['SATB2',colnames(macaque_RNA_monocle)])
macaque_RNA_monocle$TBR1 <- as.numeric(macaque_RNA_seurat@assays$RNA@data['TBR1',colnames(macaque_RNA_monocle)])
# gene_symbol <- 'PAX6'
# p1 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
#   scale_colour_gradient(low = "lightgrey",high = "blue") + 
#   theme_classic() + 
#   theme(legend.position = 'right',
#         aspect.ratio = 1,
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
#   xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'TNC'
p2 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

# gene_symbol <- 'VIM'
# p3 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
#   scale_colour_gradient(low = "lightgrey",high = "blue") + 
#   theme_classic() + 
#   theme(legend.position = 'right',
#         aspect.ratio = 1,
#         axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
#   xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'EGFR'
p4 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'OLIG2'
p5 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'SOX10'
p6 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'PCDH15'
p7 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'EOMES'
p8 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'BCL11B'
p9 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'NEUROD2'
p10 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'NEUROD6'
p11 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'SATB2'
p12 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

gene_symbol <- 'TBR1'
p13 <- plot_cell_trajectory(macaque_RNA_monocle,color_by = gene_symbol,show_branch_points = FALSE,cell_size = 0.5) + 
  scale_colour_gradient(low = "lightgrey",high = "blue") + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = gene_symbol)

pdf(file = './res/step_4_fig_210629/progenitor_cell_gene_express_trajectory.pdf',width = 8,height = 8)
p2+p4+p7+p6+p8+p13+plot_layout(nrow = 3)
dev.off()