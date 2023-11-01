#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: mouse multiome data quality control                             ##
## Data: 2022.10.14                                                                ##
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

# mouse multiome data specific marker -------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
donor_seq <- names(col_param$celltype.mouse)

#lfc 0.3 pct 0.2 overlap 0.75
marker_gene_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- my_seurat_find_specific_marker(seu.obj = mouse_multiome_Seurat,assay = 'RNA',ident.1 = x,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.3,pct.1 = 0.2,overlap_ratio = 0.75,workers = 8,features = NULL,future.globals.maxSize = 200*(1024^3))
  print(paste(x,'done!',sep = ' '))
  gc()
  return(temp)
})

saveRDS(object = marker_gene_list,file = './res/step_66_fig_221014/cell_type_specific_marker_gene_lfc_0.3_pct_0.2_overlap_0.75.rds')

#lfc 0.25 pct 0.1 overlap 0.75
marker_gene_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- my_seurat_find_specific_marker(seu.obj = mouse_multiome_Seurat,assay = 'RNA',ident.1 = x,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.25,pct.1 = 0.1,overlap_ratio = 0.75,workers = 8,features = NULL,future.globals.maxSize = 200*(1024^3))
  print(paste(x,'done!',sep = ' '))
  gc()
  return(temp)
})

saveRDS(object = marker_gene_list,file = './res/step_66_fig_221014/cell_type_specific_marker_gene_lfc_0.25_pct_0.1_overlap_0.75.rds')

#lfc 0.5 pct 0.4 overlap 0.8
marker_gene_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- my_seurat_find_specific_marker(seu.obj = mouse_multiome_Seurat,assay = 'RNA',ident.1 = x,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0.5,pct.1 = 0.4,overlap_ratio = 0.8,workers = 8,features = NULL,future.globals.maxSize = 200*(1024^3))
  print(paste(x,'done!',sep = ' '))
  gc()
  return(temp)
})

saveRDS(object = marker_gene_list,file = './res/step_66_fig_221014/cell_type_specific_marker_gene_lfc_0.5_pct_0.4_overlap_0.8.rds')

# vlnplot for basic QC ----------------------------------------------------

# #color palette
# col_list <- c('#a6dcde','#0c6c7f','#f4a684','#faa330','#f37460','#fbd8ca')
# show_col(col_list)
# names(col_list) <- c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2')
# col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
# col_param$sample.mouse <- col_list
# saveRDS(object = col_param,file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#MT ratio
MT_gene <- rownames(mouse_multiome_Seurat@assays$RNA@counts)[grep(pattern = '^mt-',x = rownames(mouse_multiome_Seurat@assays$RNA@counts),fixed = FALSE)]
mouse_multiome_Seurat[['percent.mt']] <- PercentageFeatureSet(object = mouse_multiome_Seurat,features = MT_gene,assay = 'RNA')

#col param
col_value <- col_param$sample.mouse

#plot
meta_data <- mouse_multiome_Seurat@meta.data
meta_data$donor <- factor(meta_data$donor,levels = c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2'))

p1 <- ggplot(meta_data,aes(x = donor,y = nCount_RNA,fill = donor,color = donor)) + 
  geom_violin() + 
  geom_boxplot(color = 'black',fill = 'white',width = 0.1,outlier.alpha = 0) + 
  ylim(c(0,20000)) + 
  scale_fill_manual(values = col_value[c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2')]) + 
  scale_color_manual(values = col_value[c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3) + 
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
  ylim(c(0,6000)) + 
  scale_fill_manual(values = col_value[c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2')]) + 
  scale_color_manual(values = col_value[c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2')]) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3) + 
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
  scale_fill_manual(values = col_value[c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2')]) + 
  scale_color_manual(values = col_value[c('E145_1','E145_2','E155_1','E155_2','E155_3','E165_2')]) + 
  ylim(c(0,2)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1/3) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title.y = element_text(size = 14,face = 'bold'),
        axis.line = element_line(color = 'black',size = 1),
        axis.ticks = element_blank()) + 
  ylab('# MT reads %')

pdf(file = './res/step_66_fig_221014/mouse_multiome_Seurat_QC_vlnplot.pdf',width = 7,height = 6)
p1+p2+p3+plot_layout(ncol = 1)
dev.off()

# correlation dotplot -----------------------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#add meta data
mouse_multiome_Seurat$batch <- NA
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_1','E155_1','E155_2'),"batch"] <- '220706'
mouse_multiome_Seurat@meta.data[mouse_multiome_Seurat$donor %in% c('E145_2','E155_3','E165_2'),"batch"] <- '220907'

#create pseudobulk
batch_220706 <- mouse_multiome_Seurat@assays$RNA@counts[,mouse_multiome_Seurat$batch == '220706']
batch_220706 <- log1p(rowSums(batch_220706)/sum(rowSums(batch_220706))*1000000)

batch_220907 <- mouse_multiome_Seurat@assays$RNA@counts[,mouse_multiome_Seurat$batch == '220907']
batch_220907 <- log1p(rowSums(batch_220907)/sum(rowSums(batch_220907))*1000000)

#plot
pdf(file = './res/step_66_fig_221014/mouse_multiome_data_cor_scatter_plot.pdf',width = 4,height = 4)
ggplot(data = data.frame(batch_1 = batch_220706,batch_2 = batch_220907),aes(x = batch_1,y = batch_2)) + 
  geom_point(size = 0.1,alpha = 0.8,color = '#440154FF') + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1) + 
  stat_cor(method = 'pearson') + 
  xlim(c(0,7)) + ylim(c(0,7)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 14,face = 'bold')) + 
  xlab('Batch 220706') + ylab('Batch 220907')
dev.off()

#save png
png(filename = './res/step_66_fig_221014/mouse_multiome_data_cor_scatter_plot.png',width = 350,height = 350)
ggplot(data = data.frame(batch_1 = batch_220706,batch_2 = batch_220907),aes(x = batch_1,y = batch_2)) + 
  geom_point(size = 0.1,alpha = 0.8,color = '#440154FF') + 
  xlim(c(0,7)) + ylim(c(0,7)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1) + 
  NoAxes()
dev.off()

# correlation heatmap -----------------------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
donor_seq <- names(col_param$sample.mouse)

#create pseudobulk
express_matrix <- mouse_multiome_Seurat@assays$RNA@counts
gene_list <- rownames(express_matrix)
express_matrix <- base::do.call(what = cbind,args = base::lapply(X = donor_seq,FUN = function(x){
  temp <- express_matrix[,mouse_multiome_Seurat$donor == x]
  temp <- log1p(rowSums(temp)/sum(rowSums(temp))*1000000)
  return(temp)
}))
rownames(express_matrix) <- gene_list
colnames(express_matrix) <- donor_seq

#cor matrix
cor_matrix <- cor(x = express_matrix,method = 'pearson')

#add annotation
col_anno <- HeatmapAnnotation(Batch = c('220706','220706','220706','220907','220907','220907'),
                              Age = c('E14','E14','E15','E15','E15','E16'),
                              col = list(Age = col_param$age.mouse,
                                         Batch = c('220706' = '#F47E1F','220907' = '#FCBF6E')),
                              show_annotation_name = FALSE,which = 'column')
row_anno <- HeatmapAnnotation(Batch = c('220706','220706','220706','220907','220907','220907'),
                              Age = c('E14','E14','E15','E15','E15','E16'),
                              col = list(Age = col_param$age.mouse,
                                         Batch = c('220706' = '#F47E1F','220907' = '#FCBF6E')),
                              show_annotation_name = FALSE,which = 'row')

col_fun <- colorRamp2(breaks = c(0.97,1),colors = c('white','red'))

#plot
p1 <- Heatmap(matrix = cor_matrix,show_row_names = TRUE,show_column_names = TRUE,
              cluster_rows = cluster_within_group(mat = cor_matrix,factor = c('220706','220706','220706','220907','220907','220907')),
              cluster_columns = cluster_within_group(mat = cor_matrix,factor = c('220706','220706','220706','220907','220907','220907')),
              top_annotation = col_anno,left_annotation = row_anno,row_split = 2,column_split = 2,border = TRUE,
              width = unit(5,'inches'),height = unit(5,'inches'),name = 'cor',col = col_fun)

p2 <- Heatmap(matrix = cor_matrix,show_row_names = TRUE,show_column_names = TRUE,
              cluster_rows = cluster_within_group(mat = cor_matrix,factor = c('220706','220706','220706','220907','220907','220907')),
              cluster_columns = cluster_within_group(mat = cor_matrix,factor = c('220706','220706','220706','220907','220907','220907')),
              top_annotation = col_anno,left_annotation = row_anno,row_split = 2,column_split = 2,border = TRUE,
              width = unit(5,'inches'),height = unit(5,'inches'),col = col_fun,name = 'cor',cell_fun = function(j,i,x,y,width,height,fill){
                grid.text(sprintf("%.3f",cor_matrix[i,j]),x,y)
              })

pdf(file = './res/step_66_fig_221014/mouse_multiome_data_pseudobulk_by_sample_cor_heatmap.pdf',width = 10,height = 10)
print(p1)
print(p2)
dev.off()

# mouse multiome Seurat umap plot -----------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#dimplot
DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#load color
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
temp_col <- col_param$celltype.mouse

#legend sequence
donor_seq <- names(temp_col)

#coordinate
x_lower <- min(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,1])
x_upper <- max(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,1])
y_lower <- min(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,2])
y_upper <- max(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,2])

#dimplot
p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

p3 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'cell_type',label = FALSE) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_blank()) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

pdf(file = './res/step_66_fig_221014/mouse_multiome_Seurat_cell_type_dimplot.pdf',width = 8,height = 8)
p1
p2
dev.off()

png(filename = './res/step_66_fig_221014/mouse_multiome_Seurat_cell_type_dimplot.png',width = 500,height = 500)
p3
dev.off()

#by donor
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
temp_col <- col_param$sample.mouse

donor_seq <- names(temp_col)

x_lower <- min(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,1])
x_upper <- max(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,1])
y_lower <- min(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,2])
y_upper <- max(mouse_multiome_Seurat@reductions$umap@cell.embeddings[,2])

p1 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'donor',label = FALSE) + 
  little_axis(x_range = c(x_lower,x_upper),y_range = c(y_lower,y_upper),
              ratio = 0.2,margin_value = 1,x_label = 'UMAP_1',y_label = 'UMAP_2') + 
  theme_cowplot() + 
  NoAxes() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(hjust = 0.5)) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

p2 <- DimPlot(object = mouse_multiome_Seurat,group.by = 'donor',label = FALSE) + 
  theme_cowplot() + 
  NoAxes() + NoLegend() + 
  theme(aspect.ratio = 1,
        plot.title = element_blank()) + 
  scale_color_manual(values = temp_col[donor_seq]) + 
  guides(color = guide_legend(override.aes = list(size = 3,shape = 15)))

pdf(file = './res/step_66_fig_221014/mouse_multiome_Seurat_donor_dimplot.pdf',width = 8,height = 8)
p1
dev.off()

png(filename = './res/step_66_fig_221014/mouse_multiome_Seurat_donor_dimplot.png',width = 500,height = 500)
p2
dev.off()

# marker gene dotplot -----------------------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
donor_seq <- names(col_param$celltype.mouse)

dotplot_matrix <- my_dotplot(mouse_multiome_Seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(RG=c('Sox2','Pax6','Hes5'),
                                             Cyc = c('Top2a','Mki67'),
                                             IP = c('Eomes','Btg2','Neurog2'),
                                             Ex = c('Neurod2','Neurod6','Nrp1','Satb2','Cux1','Ldb2','Fezf2','Bcl11b','Tle4','Foxp2'),
                                             In = c('Gad2','Dlx5','Lhx6','Adarb2'),
                                             End = c('Cldn5','Pecam1'),
                                             Per = c('Pdgfrb'),
                                             VLMC = c('Col1a1','Lum'),
                                             Mic = c('Cx3cr1','P2ry12','C1qb')),
                             group.by = 'cell_type', cols = c('#3B4992FF','white','#EE0000FF'),
                             return_data_plot = TRUE)
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = donor_seq)

#plot
pdf(file = './res/step_66_fig_221014/mouse_multiome_Seurat_cell_type_marker_dotplot.pdf',width = 13,height = 6)
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

# mouse multiome Seurat marker gene list heatmap --------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
donor_seq <- names(col_param$celltype.mouse)

#check marker
marker_gene_list <- readRDS(file = './data/parameter/specific_marker_list/221014/mouse_multiome/cell_type_specific_marker_gene_lfc_0.3_pct_0.2_overlap_0.75.rds')
marker_table <- base::do.call(what = rbind,args = base::lapply(X = names(marker_gene_list),FUN = function(x){
  temp <- data.frame(gene = marker_gene_list[[x]],cluster = x)
  rownames(temp) <- temp$gene
  other_cell_type <- donor_seq[donor_seq != x]
  marker_test <- my_seurat_marker_wilcox_test(seu.obj = mouse_multiome_Seurat,assay = 'RNA',ident.1 = x,ident.2 = other_cell_type,group.by = 'cell_type',min.express = 0,log2fc_thresholf = 0,pct.1 = 0,features = temp$gene,only.pos = TRUE,workers = 1,future.globals.maxSize = 2200*(1024^3))
  marker_test <- marker_test[rownames(temp),]
  temp <- cbind(temp,marker_test)
  print(paste(x,'done!',sep = ' '))
  return(temp)
}))

#classic marker
classis_marker <- readRDS(file = './data/parameter/shared_param/Celltype-Canonical-gene-list.rds')
classis_marker <- unlist(classis_marker)
classis_marker <- classis_marker[!(classis_marker %in% c('ETV1','MEIS2'))]

mouse_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_anno <- mouse_anno[,c(4,2)]
mouse_anno <- mouse_anno[!duplicated(mouse_anno[,1]),]
rownames(mouse_anno) <- mouse_anno[,1]
classis_marker <- mouse_anno[classis_marker,2]
classis_marker <- classis_marker[!is.na(classis_marker)]
classis_marker <- unique(classis_marker)

classis_marker <- append(classis_marker,c('Sox2','Pax6','Hes5','Top2a','Mki67','Eomes','Btg2','Neurog2',
                                          'Neurod2','Neurod6','Nrp1','Satb2','Cux1','Ldb2','Fezf2','Bcl11b',
                                          'Tle4','Foxp2','Gad2','Dlx5','Lhx6','Adarb2','Cldn5','Pecam1','Pdgfrb',
                                          'Col1a1','Lum','Cx3cr1','P2ry12','C1qb'))
classis_marker <- unique(classis_marker)

#filter gene
classis_marker <- classis_marker[!(classis_marker %in% c('Rbfox1','Igfbp7','Sparc','Nr2f2','Eya2'))]

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
marker_table_subset <- marker_table_subset[!(rownames(marker_table_subset) %in% c('Ldb2','Cux1','Satb2')),]

#add marker
marker_table_subset <- rbind(marker_table_subset,data.frame(gene = c('Ldb2','Cux1','Satb2','Lhx6'),
                                                            cluster = c('SCPN','CPN','CPN','In'),
                                                            log2FC = 100,pval = 0,fdr = 0,mean1 = 100,mean2 = 0,
                                                            pct.1 = 1,pct.2 = 0,diff = 100))
#subset cells
table(mouse_multiome_Seurat$cell_type)
cell_list <- base::lapply(X = donor_seq,FUN = function(x){
  temp <- colnames(mouse_multiome_Seurat)[mouse_multiome_Seurat$cell_type == x]
  if(length(temp) < 200){
    temp <- sample(x = temp,size = 200,replace = TRUE)
  }else{
    temp <- sample(x = temp,size = 200,replace = FALSE)
  }
  return(temp)
})
cell_list <- unlist(cell_list)
names(cell_list) <- mouse_multiome_Seurat@meta.data[cell_list,"cell_type"]

#express matrix
express_matrix <- mouse_multiome_Seurat@assays$RNA@data[marker_table_subset$gene,cell_list]
express_matrix <- t(scale(x = t(express_matrix)))

#annotation
col_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
temp_col <- col_param$celltype.mouse
col_anno <- HeatmapAnnotation(cell_type = factor(names(cell_list),levels = donor_seq),
                              col = list(cell_type = temp_col[c(donor_seq)]),
                              which = 'column')

temp <- which(rownames(express_matrix) %in% unlist(classis_marker))
row_anno <- HeatmapAnnotation(marker = anno_mark(at = temp,labels = rownames(express_matrix)[temp]),which = 'row')

col_fun <- colorRamp2(breaks = c(-2,-0.5,0.5,2),colors = c('#3B4992FF','white','white','#EE0000FF'))

#plot
pdf(file = './res/step_66_fig_221014/mouse_multiome_Seurat_cell_type_marker_list_heatmap.pdf',width = 10,height = 10)
Heatmap(matrix = express_matrix,cluster_rows = FALSE,cluster_columns = FALSE,show_column_names = FALSE,show_row_names = FALSE,
        top_annotation = col_anno,col = col_fun,use_raster = TRUE,column_split = factor(names(cell_list),levels = donor_seq),
        row_split = factor(marker_table_subset$cluster,levels = donor_seq),right_annotation = row_anno,name = 'z-score',
        height = unit(8,'inches'),width = unit(5.5,'inches'))
dev.off()

# mouse multiome data feature plot ----------------------------------------
#load data
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

#classic marker
classis_marker <- readRDS(file = './data/parameter/shared_param/Celltype-Canonical-gene-list.rds')
classis_marker <- unlist(classis_marker)
classis_marker <- classis_marker[!(classis_marker %in% c('ETV1','MEIS2'))]

mouse_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCm39_to_GRCh38.csv')
mouse_anno <- mouse_anno[,c(4,2)]
mouse_anno <- mouse_anno[!duplicated(mouse_anno[,1]),]
rownames(mouse_anno) <- mouse_anno[,1]
classis_marker <- mouse_anno[classis_marker,2]
classis_marker <- classis_marker[!is.na(classis_marker)]
classis_marker <- unique(classis_marker)

classis_marker <- append(classis_marker,c('Sox2','Pax6','Hes5','Top2a','Mki67','Eomes','Btg2','Neurog2',
                                          'Neurod2','Neurod6','Nrp1','Satb2','Cux1','Ldb2','Fezf2','Bcl11b',
                                          'Tle4','Foxp2','Gad2','Dlx5','Lhx6','Adarb2','Cldn5','Pecam1','Pdgfrb',
                                          'Col1a1','Lum','Cx3cr1','P2ry12','C1qb'))
classis_marker <- unique(classis_marker)

#feature plot legend
pdf(file = './res/step_66_fig_221014/mouse_multiome_Seurat_marker_featureplot_legend.pdf',width = 8,height = 7)
FeaturePlot(features = 'Pax6',order = TRUE,object = mouse_multiome_Seurat,cols = c('#e0ecf4','#8c6bb1','#4d004b'),pt.size = 0.1,max.cutoff = 'q99',slot = 'data') + 
  theme(aspect.ratio = 1) + 
  NoAxes() + 
  theme(plot.title = element_text(face = 'italic'))
dev.off()

#batch plot
marker_gene_list <- classis_marker
marker_gene_list <- marker_gene_list[marker_gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts)]

#with order
for (i in marker_gene_list) {
  char <- paste0('./res/step_66_fig_221014/feature_plot/mouse_multiome_Seurat_',i,'_featureplot_ordered.png')
  p <- FeaturePlot(features = i,order = TRUE,object = mouse_multiome_Seurat,cols = c('#e0ecf4','#8c6bb1','#4d004b'),pt.size = 0.001,max.cutoff = 'q99',slot = 'data') + 
    NoAxes() + 
    theme(plot.title = element_blank(),
          aspect.ratio = 1,legend.position = 'none')
  png(filename = char,width = 500,height = 500)
  print(p)
  dev.off()
}

#without order
for (i in marker_gene_list) {
  char <- paste0('./res/step_66_fig_221014/feature_plot/mouse_multiome_Seurat_',i,'_featureplot_without_ordered.png')
  p <- FeaturePlot(features = i,order = FALSE,object = mouse_multiome_Seurat,cols = c('#e0ecf4','#8c6bb1','#4d004b'),pt.size = 0.001,max.cutoff = 'q99',slot = 'data') + 
    NoAxes() + 
    theme(plot.title = element_blank(),
          aspect.ratio = 1,legend.position = 'none')
  png(filename = char,width = 500,height = 500)
  print(p)
  dev.off()
}
