#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque snRNA_seq quality control fig                           ##
## Data: 2022.09.13                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

# notice:
#donor sequence: A50A/A50B (E80), A82A/A82B (E90), A84B/A84C (E84), A68A/A68B (E92)

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

# vlnplt for basic QC -----------------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')

#MT ratio
# macaque_to_human_anno <- read.csv(file = './data/reference/BioMart_release_105/Mmul10_to_GRCh38.csv')
# macaque_to_human_anno <- macaque_to_human_anno[grep(pattern = '^MT-',x = macaque_to_human_anno[,4],fixed = FALSE),]
# macaque_to_human_anno <- macaque_to_human_anno[!(macaque_to_human_anno[,2] == ''),]
# MT_gene <- as.character(macaque_to_human_anno[,2])
# MT_gene <- unique(MT_gene)
# write(x = MT_gene,file = './data/reference/Macque_mitochondria_gene_list.txt')
MT_gene <- readLines(con = './data/reference/Macque_mitochondria_gene_list.txt')
macaque_RNA_Seurat[['percent.mt']] <- PercentageFeatureSet(object = macaque_RNA_Seurat,features = MT_gene,assay = 'RNA')

#col param
col_value <- col_param$sample

#plot
meta_data <- macaque_RNA_Seurat@meta.data
meta_data$donor <- factor(meta_data$donor,levels = c('A50A','A82B','A84B','A84C','A68A','A68B'))

p1 <- ggplot(meta_data,aes(x = donor,y = nCount_RNA,fill = donor,color = donor)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A82B','A84B','A84C','A68A','A68B')]) + 
  scale_color_manual(values = col_value[c('A50A','A82B','A84B','A84C','A68A','A68B')]) + ylim(c(0,15000)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1)) + 
  ylab('# Reads')

p2 <- ggplot(meta_data,aes(x = donor,y = nFeature_RNA,fill = donor,color = donor)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A82B','A84B','A84C','A68A','A68B')]) + 
  scale_color_manual(values = col_value[c('A50A','A82B','A84B','A84C','A68A','A68B')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1)) + 
  ylab('# Genes')

p3 <- ggplot(meta_data,aes(x = donor,y = percent.mt,fill = donor,color = donor)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  scale_fill_manual(values = col_value[c('A50A','A82B','A84B','A84C','A68A','A68B')]) + 
  scale_color_manual(values = col_value[c('A50A','A82B','A84B','A84C','A68A','A68B')]) + 
  ylim(c(0,1)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1),
        axis.ticks = element_blank()) + 
  ylab('# MT reads %')

pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_QC_vlnplot.pdf',width = 4,height = 6)
p1+p2+p3+plot_layout(ncol = 1)
dev.off()

# correlation plot --------------------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')
table(macaque_RNA_Seurat$batch)

#two batch
batch_1 <- macaque_RNA_Seurat[,macaque_RNA_Seurat$batch == '200919']
batch_1 <- batch_1@assays$RNA@counts
batch_1 <- log1p(rowSums(batch_1)/sum(rowSums(batch_1))*1000000)

batch_2 <- macaque_RNA_Seurat[,macaque_RNA_Seurat$batch == '210922']
batch_2 <- batch_2@assays$RNA@counts
batch_2 <- log1p(rowSums(batch_2)/sum(rowSums(batch_2))*1000000)

#plot
pdf(file = './res/step_57_fig_220913/snRNA_seq_batch_cor_scatter_plot.pdf',width = 4,height = 4)
ggplot(data = data.frame(batch_1 = batch_1,batch_2 = batch_2),aes(x = batch_1,y = batch_2)) + 
  geom_point(size = 0.1,alpha = 0.8,color = '#440154FF') + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1) + 
  stat_cor(method = 'pearson') + 
  xlim(c(0,7)) + ylim(c(0,7)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14,face = 'bold')) + 
  xlab('Batch 200919') + ylab('Batch 210922')
dev.off()

#correlation heatmap matrix
donor_list <- c('A50A','A82B','A84B','A84C','A68A','A68B')
express_matrix <- do.call(what = cbind,args = base::lapply(X = donor_list,FUN = function(x){
  temp <- macaque_RNA_Seurat[,macaque_RNA_Seurat$donor == x]
  temp <- temp@assays$RNA@counts
  temp <- log1p(rowSums(temp)/sum(rowSums(temp))*1000000)
  return(temp)
}))
colnames(express_matrix) <- donor_list
rownames(express_matrix) <- rownames(macaque_RNA_Seurat@assays$RNA@counts)
cor_matrix <- cor(x = express_matrix,method = 'pearson')

#color
mybreak <- seq(0.97,1, 0.001)
mycols <- colorRampPalette(c('white', 'red'))(length(mybreak))

#plot
r1 <- pheatmap::pheatmap(cor_matrix, display_numbers = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA)
r2 <- pheatmap::pheatmap(cor_matrix, display_numbers = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA)
r3 <- pheatmap::pheatmap(cor_matrix, display_numbers = TRUE, cluster_rows = TRUE, cluster_cols = TRUE, color = mycols, breaks = mybreak, border_color = NA)
r4 <- pheatmap::pheatmap(cor_matrix, display_numbers = FALSE, cluster_rows = TRUE, cluster_cols = TRUE, color = mycols, breaks = mybreak, border_color = NA)

pdf(file = './res/step_57_fig_220913/macaque_RNA_pseudobulk_by_sample_cor_heatmap_liuyt_style.pdf',width = 4.5,height = 4)
r1
plot.new()
r2
plot.new()
r3
plot.new()
r4
dev.off()

# correlation heatmap my own style ----------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')

#express matrix
donor_list <- c('A50A','A82B','A84B','A84C','A68A','A68B')
express_matrix <- do.call(what = cbind,args = base::lapply(X = donor_list,FUN = function(x){
  temp <- macaque_RNA_Seurat[,macaque_RNA_Seurat$donor == x]
  temp <- temp@assays$RNA@counts
  temp <- log1p(rowSums(temp)/sum(rowSums(temp))*1000000)
  return(temp)
}))
colnames(express_matrix) <- donor_list
rownames(express_matrix) <- rownames(macaque_RNA_Seurat@assays$RNA@counts)

#cor matrix
cor_matrix <- cor(x = express_matrix,method = 'pearson')

#plot
col_anno <- HeatmapAnnotation(batch = c('Batch 210922','Batch 210922','Batch 200919','Batch 200919','Batch 200919','Batch 200919'),
                              show_annotation_name = FALSE,which = 'column',
                              col = list(batch = c('Batch 200919' = '#F47E1F','Batch 210922' = '#FCBF6E')))
row_anno <- HeatmapAnnotation(batch = c('Batch 210922','Batch 210922','Batch 200919','Batch 200919','Batch 200919','Batch 200919'),
                              show_annotation_name = FALSE,show_legend = FALSE,which = 'row',
                              col = list(batch = c('Batch 200919' = '#F47E1F','Batch 210922' = '#FCBF6E')))
col_fun <- colorRamp2(breaks = c(0.95,1),colors = c('white','red'))

p1 <- Heatmap(matrix = cor_matrix,show_row_names = TRUE,show_column_names = TRUE,
              cluster_rows = cluster_within_group(mat = cor_matrix,factor = c('Batch 210922','Batch 210922','Batch 200919','Batch 200919','Batch 200919','Batch 200919')),
              cluster_columns = cluster_within_group(mat = cor_matrix,factor = c('Batch 210922','Batch 210922','Batch 200919','Batch 200919','Batch 200919','Batch 200919')),
              top_annotation = col_anno,left_annotation = row_anno,row_split = 2,column_split = 2,border = TRUE,
              width = unit(3,'inches'),height = unit(3,'inches'),col = col_fun,name = 'cor')
p2 <- Heatmap(matrix = cor_matrix,show_row_names = TRUE,show_column_names = TRUE,
              cluster_rows = cluster_within_group(mat = cor_matrix,factor = c('Batch 210922','Batch 210922','Batch 200919','Batch 200919','Batch 200919','Batch 200919')),
              cluster_columns = cluster_within_group(mat = cor_matrix,factor = c('Batch 210922','Batch 210922','Batch 200919','Batch 200919','Batch 200919','Batch 200919')),
              top_annotation = col_anno,left_annotation = row_anno,row_split = 2,column_split = 2,border = TRUE,
              width = unit(3,'inches'),height = unit(3,'inches'),col = col_fun,name = 'cor',cell_fun = function(j,i,x,y,width,height,fill){
                grid.text(sprintf("%.3f",cor_matrix[i,j]),x,y)
              })

pdf(file = './res/step_57_fig_220913/macaque_RNA_pseudobulk_by_sample_cor_heatmap_sunym_style.pdf',width = 8,height = 6)
print(p1)
print(p2)
dev.off()

# macaque snRNA-seq UMAP plot ---------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')

#dimplot
DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#load color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype

#legend sequence
donor_seq <- c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End','Per','VLMC','Mic')

#coordinate
x_lower <- min(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,1])
x_upper <- max(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,1])
y_lower <- min(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,2])
y_upper <- max(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,2])

#dimplot
p1 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

p2 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

p3 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'cell_type',label = FALSE) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_cell_type_dimplot.pdf',width = 8,height = 8)
p1
p2
p3
dev.off()

#show liuyt
my_dimplot(embedding = macaque_RNA_Seurat@reductions$umap@cell.embeddings,
           meta_data = as.data.frame(macaque_RNA_Seurat@meta.data),
           group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-1','RG-2','OPC','Mic','End','Per','VLMC')])

#by donor
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$sample

donor_seq <- c('A50A','A82B','A84B','A84C','A68A','A68B')

x_lower <- min(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,1])
x_upper <- max(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,1])
y_lower <- min(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,2])
y_upper <- max(macaque_RNA_Seurat@reductions$umap@cell.embeddings[,2])

p1 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'donor',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

p2 <- DimPlot(object = macaque_RNA_Seurat,group.by = 'donor',label = FALSE) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_donor_dimplot.pdf',width = 8,height = 8)
p1
p2
dev.off()

# marker gene dotplot -----------------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')
donor_seq <- c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End','Per','VLMC','Mic')

#dot matrix
# dotplot_matrix <- my_dotplot(macaque_RNA_Seurat,assay = 'RNA', 
#                              col.max = 2.5, col.min = -2.5, scale = TRUE, 
#                              features = list(VLMC = c('COL1A1','LUM'),
#                                              Per=c('PDGFRB'),
#                                              End=c('CLDN5','PECAM1'),
#                                              Mic=c('CX3CR1'),
#                                              OPC=c('SOX10','OLIG2'),
#                                              RG=c('SOX9','PAX6','AQP4','EGFR'),
#                                              Cyc=c('TOP2A','MKI67'),
#                                              IP=c('EOMES','PPP1R17'),
#                                              Ex=c('NEUROD2','NRP1','SATB2','STMN2','NEUROD6','LMO4'),
#                                              SP=c('NR4A2','CRYM','CDH18'),
#                                              In=c('DLX5','GAD2','SP8','NR2F2','SST','LHX6')),
#                              group.by = 'cell_type', cols = c('#3B4992FF','white','#EE0000FF'),
#                              return_data_plot = TRUE)
dotplot_matrix <- my_dotplot(macaque_RNA_Seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(RG=c('SLC1A3','SOX9','PAX6','PTPRZ1','EGFR'),
                                             Cyc=c('TOP2A','MKI67'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NRP1','SATB2','STMN2','NEUROD6','LMO4'),
                                             SP=c('NR4A2','CRYM','CDH18'),
                                             In=c('GAD1','DLX5','ADARB2','LHX6'),
                                             OPC=c('SOX10','OLIG2'),
                                             End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             VLMC = c('COL1A1','LUM'),
                                             Mic=c('CX3CR1')),
                             group.by = 'cell_type', cols = c('#3B4992FF','white','#EE0000FF'),
                             return_data_plot = TRUE)
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = donor_seq)

#plot
pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_cell_type_marker_dotplot.pdf',width = 13,height = 6)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom') + 
  xlab('') + ylab('')
dev.off()

# marker gene list heatmap ------------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')
donor_seq <- c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End','Per','VLMC','Mic')

#lfc 0.3 pct 0.2 overlap 0.75
marker_gene_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- my_seurat_find_specific_marker(seu.obj = macaque_RNA_Seurat,assay = 'RNA',ident.1 = x,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.3,pct.1 = 0.2,overlap_ratio = 0.75,workers = 4,features = NULL,future.globals.maxSize = 200*(1024^3))
  print(paste(x,'done!',sep = ' '))
  gc()
  return(temp)
})

saveRDS(object = marker_gene_list,file = './res/step_57_fig_220913/cell_type_specific_marker_gene_lfc_0.3_pct_0.2_overlap_0.75.rds')

#lfc 0.25 pct 0.1 overlap 0.75
marker_gene_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- my_seurat_find_specific_marker(seu.obj = macaque_RNA_Seurat,assay = 'RNA',ident.1 = x,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.25,pct.1 = 0.1,overlap_ratio = 0.75,workers = 4,features = NULL,future.globals.maxSize = 200*(1024^3))
  print(paste(x,'done!',sep = ' '))
  gc()
  return(temp)
})

saveRDS(object = marker_gene_list,file = './res/step_57_fig_220913/cell_type_specific_marker_gene_lfc_0.25_pct_0.1_overlap_0.75.rds')

#lfc 0.5 pct 0.4 overlap 0.8
marker_gene_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- my_seurat_find_specific_marker(seu.obj = macaque_RNA_Seurat,assay = 'RNA',ident.1 = x,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.5,pct.1 = 0.4,overlap_ratio = 0.8,workers = 4,features = NULL,future.globals.maxSize = 200*(1024^3))
  print(paste(x,'done!',sep = ' '))
  gc()
  return(temp)
})

saveRDS(object = marker_gene_list,file = './res/step_57_fig_220913/cell_type_specific_marker_gene_lfc_0.5_pct_0.4_overlap_0.8.rds')

#check marker
marker_gene_list <- readRDS(file = './res/step_57_fig_220913/cell_type_specific_marker_gene_lfc_0.3_pct_0.2_overlap_0.75.rds')
marker_table <- base::do.call(what = rbind,args = base::lapply(X = names(marker_gene_list),FUN = function(x){
  temp <- data.frame(gene = marker_gene_list[[x]],cluster = x)
  rownames(temp) <- temp$gene
  other_cell_type <- donor_seq[donor_seq != x]
  marker_test <- my_seurat_marker_wilcox_test(seu.obj = macaque_RNA_Seurat,assay = 'RNA',ident.1 = x,ident.2 = other_cell_type,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0,pct.1 = 0,features = temp$gene,only.pos = TRUE,workers = 1,future.globals.maxSize = 2200*(1024^3))
  marker_test <- marker_test[rownames(temp),]
  temp <- cbind(temp,marker_test)
  print(paste(x,'done!',sep = ' '))
  return(temp)
}))

#classic marker
classis_marker <- readRDS(file = './data/parameter/shared_param/Celltype-Canonical-gene-list.rds')
classis_marker <- unlist(classis_marker)
classis_marker <- classis_marker[!(classis_marker %in% c('ETV1','MEIS2'))]
classis_marker <- append(classis_marker,c('SLC1A3','SOX9','PAX6','PTPRZ1','EGFR','NEUROD2',
                                          'NRP1','SATB2','STMN2','NEUROD6','LMO4','NR4A2',
                                          'CRYM','CDH18','ADARB2','SOX10','OLIG2','COL1A1','LUM'))

#check duplicate
table(duplicated(marker_table$gene))
marker_table <- marker_table[order(marker_table$fdr,decreasing = FALSE),]
marker_table <- marker_table[!(duplicated(marker_table$gene)),]
table(duplicated(marker_table$gene))
table(marker_table$cluster)
rownames(marker_table) <- marker_table$gene

#subset marker table
marker_table_subset <- base::do.call(what = rbind,args = base::lapply(X = donor_seq,FUN = function(x){
  temp <- marker_table[marker_table$cluster == x,]
  temp$diff <- temp$mean1 - temp$mean2
  temp <- temp[order(temp$diff,decreasing = TRUE),]
  if(nrow(temp) <= 100){
    return(temp)
  }else{
    return(temp[1:100,])
  }
}))
rownames(marker_table_subset) <- marker_table_subset$gene
marker_table_subset <- marker_table_subset[order(marker_table_subset$cluster,decreasing = FALSE),]

#add marker
marker_table_subset <- rbind(marker_table_subset,data.frame(gene = c('NEUROD6','LMO4','OLIG2','SOX9','PAX6'),
                                                            cluster = c('Ex-3','Ex-3','OPC','RG-2','RG-1'),
                                                            log2FC = 100,pval = 0,fdr = 0,mean1 = 100,mean2 = 0,
                                                            pct.1 = 1,pct.2 = 0,diff = 100))
#subset cells
table(macaque_RNA_Seurat$cell_type)
cell_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- colnames(macaque_RNA_Seurat)[macaque_RNA_Seurat$cell_type == x]
  if(length(temp) < 200){
    temp <- sample(x = temp,size = 200,replace = TRUE)
  }else{
    temp <- sample(x = temp,size = 200,replace = FALSE)
  }
  return(temp)
})
cell_list <- unlist(cell_list)
names(cell_list) <- macaque_RNA_Seurat@meta.data[cell_list,"cell_type"]

#express matrix
express_matrix <- macaque_RNA_Seurat@assays$RNA@data[marker_table_subset$gene,cell_list]
express_matrix <- t(scale(x = t(express_matrix)))

#annotation
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype
col_anno <- HeatmapAnnotation(cell_type = factor(names(cell_list),levels = donor_seq),
                              col = list(cell_type = temp_col[c(donor_seq)]),
                              which = 'column')

temp <- which(rownames(express_matrix) %in% unlist(classis_marker))
row_anno <- HeatmapAnnotation(marker = anno_mark(at = temp,labels = rownames(express_matrix)[temp]),which = 'row')

col_fun <- colorRamp2(breaks = c(-2,-0.5,0.5,2),colors = c('#3B4992FF','white','white','#EE0000FF'))

#plot
pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_cell_type_marker_list_heatmap.pdf',width = 10,height = 10)
Heatmap(matrix = express_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,show_row_names = FALSE,
        top_annotation = col_anno,col = col_fun,use_raster = TRUE,column_split = factor(names(cell_list),levels = donor_seq),
        row_split = factor(marker_table_subset$cluster,levels = donor_seq),right_annotation = row_anno,name = 'z-score',
        height = unit(8,'inches'),width = unit(5.5,'inches'))
dev.off()

# marker feature plot -----------------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')
marker_gene_list <- readRDS(file = './res/step_57_fig_220913/cell_type_specific_marker_gene_lfc_0.3_pct_0.2_overlap_0.75.rds')

#feature plot legend
pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_marker_featureplot_legend.pdf',width = 8,height = 7)
FeaturePlot(features = 'PAX6',order = TRUE,object = macaque_RNA_Seurat,cols = c('#e0ecf4','#8c6bb1','#4d004b'),pt.size = 0.1,max.cutoff = 'q99',slot = 'data') + 
  theme(aspect.ratio = 1) + 
  NoAxes() + 
  theme(plot.title = element_text(face = 'italic'))
dev.off()

#batch plot
marker_gene_list <- c('SLC1A3','SOX9','PAX6','PTPRZ1','EGFR','TOP2A','MKI67','EOMES','PPP1R17',
                      'NEUROD2','NRP1','SATB2','STMN2','NEUROD6','LMO4','NR4A2','CRYM','CDH18',
                      'GAD1','DLX5','ADARB2','LHX6','SOX10','OLIG2','CLDN5','PECAM1','PDGFRB',
                      'COL1A1','LUM','CX3CR1')
temp <- readRDS(file = './data/parameter/shared_param/Celltype-Canonical-gene-list.rds')
temp <- unique(unlist(temp))
marker_gene_list <- append(marker_gene_list,temp)
marker_gene_list <- unique(marker_gene_list)
marker_gene_list <- marker_gene_list[marker_gene_list %in% rownames(macaque_RNA_Seurat@assays$RNA@data)]

#with order
for (i in marker_gene_list) {
  char <- paste0('./res/step_57_fig_220913/feature_plot/macaque_RNA_Seurat_',i,'_featureplot_ordered.pdf')
  p <- FeaturePlot(features = i,order = TRUE,object = macaque_RNA_Seurat,cols = c('#e0ecf4','#8c6bb1','#4d004b'),pt.size = 0.001,max.cutoff = 'q99',slot = 'data') + 
    NoAxes() + 
    theme(plot.title = element_text(face = 'bold',size = 10,hjust = 0.5),
          aspect.ratio = 1,legend.position = 'none',panel.background = element_rect(fill = NA,colour = 'black',size = 0.8))
  pdf(file = char,width = 3,height = 4)
  print(p)
  dev.off()
}

#with order
for (i in marker_gene_list) {
  char <- paste0('./res/step_57_fig_220913/feature_plot/macaque_RNA_Seurat_',i,'_featureplot_without_ordered.pdf')
  p <- FeaturePlot(features = i,order = FALSE,object = macaque_RNA_Seurat,cols = c('#e0ecf4','#8c6bb1','#4d004b'),pt.size = 0.001,max.cutoff = 'q99',slot = 'data') + 
    NoAxes() + 
    theme(plot.title = element_text(face = 'bold',size = 10,hjust = 0.5),
          aspect.ratio = 1,legend.position = 'none',panel.background = element_rect(fill = NA,colour = 'black',size = 0.8))
  pdf(file = char,width = 3,height = 4)
  print(p)
  dev.off()
}

# donor and age distribution ----------------------------------------------
#load data
macaque_RNA_Seurat <- readRDS(file = './processed_data/220802_summary/macaque_RNA_Seurat_220803.rds')
donor_proportion_matrix <- My_Cell_Proportion(meta_data = macaque_RNA_Seurat@meta.data,group.by = 'donor',split.by = 'cell_type')
donor_list <- c('A50A','A82B','A84B','A84C','A68A','A68B')
cell_type_list <- c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC','End','Per','VLMC','Mic')

donor_proportion_matrix$donor <- factor(donor_proportion_matrix$donor,levels = rev(donor_list))
donor_proportion_matrix$cell_type <- factor(donor_proportion_matrix$cell_type,levels = rev(cell_type_list))

#load color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$sample

p <- ggplot(donor_proportion_matrix, aes(fill = donor,y = Proportion, x = cell_type)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = temp_col[donor_list]) + 
  xlab('Cell type') + 
  scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0),name = 'Percentage') + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),aspect.ratio = 2.5) + 
  guides(fill=guide_legend(title="Donors")) + 
  coord_flip()

pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_cell_type_donor_proportion_barplot.pdf',width = 4,height = 4.5)
p
dev.off()

#age proportion A50A/A50B (E80), A82A/A82B (E90), A84B/A84C (E84), A68A/A68B (E92)
macaque_RNA_Seurat$Age <- NA
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A50A','A50B'),"Age"] <- 'E80'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A82A','A82B'),"Age"] <- 'E90'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A84B','A84C'),"Age"] <- 'E84'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A68A','A68B'),"Age"] <- 'E92'
age_proportion_matrix <- My_Cell_Proportion(meta_data = macaque_RNA_Seurat@meta.data,group.by = 'Age',split.by = 'cell_type')

age_proportion_matrix$cell_type <- factor(x = age_proportion_matrix$cell_type,levels = rev(cell_type_list))
age_proportion_matrix$Age <- factor(x = age_proportion_matrix$Age,levels = rev(c('E80','E84','E90','E92')))

#load color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$age.macaque

p <- ggplot(age_proportion_matrix, aes(fill = Age,y = Proportion, x = cell_type)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = temp_col[c('E80','E84','E90','E92')]) + 
  xlab('Age') + 
  scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0),name = 'Percentage') + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),aspect.ratio = 2.5) + 
  guides(fill=guide_legend(title="Age")) + 
  coord_flip()

pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_cell_type_Age_proportion_barplot.pdf',width = 4,height = 4.5)
p
dev.off()

#cell type proportion by Age
macaque_RNA_Seurat$Age <- NA
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A50A','A50B'),"Age"] <- 'E80'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A82A','A82B'),"Age"] <- 'E90'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A84B','A84C'),"Age"] <- 'E84'
macaque_RNA_Seurat@meta.data[macaque_RNA_Seurat$donor %in% c('A68A','A68B'),"Age"] <- 'E92'
age_proportion_matrix <- My_Cell_Proportion(meta_data = macaque_RNA_Seurat@meta.data,group.by = 'cell_type',split.by = 'Age')

age_proportion_matrix$cell_type <- factor(x = age_proportion_matrix$cell_type,levels = rev(cell_type_list))
age_proportion_matrix$Age <- factor(x = age_proportion_matrix$Age,levels = rev(c('E80','E84','E90','E92')))

#load color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220915.rds')
temp_col <- col_param$celltype

p <- ggplot(age_proportion_matrix, aes(fill = cell_type,y = Proportion, x = Age)) +
  geom_bar(position = "stack", stat = "identity") + 
  scale_fill_manual(values = temp_col[cell_type_list]) + 
  xlab('Age') + 
  scale_y_continuous(breaks = seq(0,1,0.1),expand = c(0,0),name = 'Percentage') + 
  theme_bw() + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),aspect.ratio = 0.4,
        legend.position = 'bottom') + 
  guides(fill=guide_legend(title="Cell Type")) + 
  coord_flip()

pdf(file = './res/step_57_fig_220913/macaque_RNA_Seurat_cell_type_proportion_by_Age_barplot.pdf',width = 5,height = 3.5)
p
dev.off()