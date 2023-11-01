#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque RNA seq fig 1                                           ##
## Data: 2021.06.24                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210623.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
#scatter plot
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
# ggsci::pal_npg()(10)
# show_col(ggsci::pal_npg()(10))
# show_col(ggsci::pal_aaas()(10))
# show_col(ggsci::pal_lancet()(10))
# 
# color_paramater <- data.frame(id = c('vRG','oRG','OPC','Cycling','Astrocyte','IP','Ex','Ex-U','Ex-SP','InCGE','InMGE','InPSB','Mic','Per','End'),
#                               col = c('#4DBBD5FF','#00A087FF','#F39B7FFF','#3C5488FF','#E64B35FF','#8491B4FF','#91D1C2FF','#008280FF','#00468BFF','#7E6148FF','#B09C85FF','#DC0000FF','#631879FF','#5F559BFF','#A20056FF'))
# write.csv(color_paramater,file = './data/parameter/color_paramater.csv')

color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id
pdf(file = './res/step_3_fig_210624/macaque_scRNAseq_scatter_plot.pdf',width = 5,height = 5)
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cols = temp_col[c('vRG','oRG','OPC','Cycling','Astrocyte','IP','Ex','Ex-U','Ex-SP','InCGE','InMGE','InPSB','Mic','Per','End')]) + 
  theme_classic() + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.position = 'none',
        aspect.ratio = 1) + 
  labs(title = '')
dev.off()

#marker dotplot
show_col(ggsci::pal_aaas()(10))
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
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('End','Per','Mic','Astrocyte','vRG','oRG','OPC','Cycling','IP','Ex','Ex-U','Ex-SP','InCGE','InMGE','InPSB'))

pdf(file = './res/step_3_fig_210624/dotplot_marker_macaque_RNAseq.pdf',width = 12,height = 6)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'none') + 
  xlab('') + ylab('')
dev.off()

pdf(file = './res/step_3_fig_210624/dotplot_marker_macaque_RNAseq_legend.pdf',width = 9.5,height = 5)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_blank(),
        legend.position = 'bottom',
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) + 
  xlab('') + ylab('')
dev.off()

#Proportion plot
cell_type_donor_proportion <- My_Cell_Proportion(macaque_RNA_seurat,group.by = 'donor',split.by = 'cell_type')
cell_type_donor_proportion$cell_type <- factor(cell_type_donor_proportion$cell_type,levels = c('End','Per','Mic','Astrocyte','vRG','oRG','OPC','Cycling','IP','Ex','Ex-U','Ex-SP','InCGE','InMGE','InPSB'))
show_col(ggsci::pal_d3(alpha = 0.75)(10))
# color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
# temp <- data.frame(id = c('A68A','A68B','A84B','A84C'),
#                    col = ggsci::pal_d3(alpha = 0.75)(10)[1:4])
# color_paramater <- rbind(color_paramater,temp)
# 
# write.csv(color_paramater,file = './data/parameter/color_paramater.csv')
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)

temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

pdf(file = './res/step_3_fig_210624/macaque_RNA_seurat_donor_proportion.pdf',width = 4,height = 5)
ggplot(data = cell_type_donor_proportion, mapping = aes(x = cell_type,y = Proportion,fill = donor))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.95) + 
  scale_fill_manual(values = temp_col[c('A68A','A68B','A84B','A84C')]) + 
  theme_classic() + 
  theme(aspect.ratio = 2.5,
        axis.text.x = element_text(size = 12,angle = 45,vjust = 1,hjust = 1),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5)) + 
  labs(fill = 'Donor') + xlab('') + ylab('Percent of Cells') + 
  coord_flip()
dev.off()

# QC plot
meta_data <- macaque_RNA_seurat@meta.data
meta_data$nCount_RNA <- as.numeric(meta_data$nCount_RNA)
meta_data$nFeature_RNA <- as.numeric(meta_data$nFeature_RNA)
meta_data$percent.mt <- as.numeric(meta_data$percent.mt)
meta_data$donor <- factor(meta_data$donor,levels = c('A68A','A68B','A84B','A84C'))
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

p1 <- ggplot(data = meta_data,aes(x=donor,y=nCount_RNA,fill=donor)) + 
  geom_violin(trim = TRUE) + 
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange", color = "black") + 
  scale_fill_manual(values = temp_col[c('A68A','A68B','A84B','A84C')]) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.position = 'none',
        aspect.ratio = 0.5) + 
  xlab('') + ylab('UMI') + labs(title = '')

p2 <- ggplot(data = meta_data,aes(x=donor,y=nFeature_RNA,fill=donor)) + 
  geom_violin(trim = TRUE) + 
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange", color = "black") + 
  scale_fill_manual(values = temp_col[c('A68A','A68B','A84B','A84C')]) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.position = 'none',
        aspect.ratio = 0.5) + 
  xlab('') + ylab('Genes detected') + labs(title = '')

p3 <- ggplot(data = meta_data,aes(x=donor,y=percent.mt,fill=donor)) + 
  geom_violin(trim = TRUE) + 
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="pointrange", color = "black") + 
  scale_fill_manual(values = temp_col[c('A68A','A68B','A84B','A84C')]) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.position = 'none',
        aspect.ratio = 0.5) + 
  xlab('') + ylab('Percent MT') + labs(title = '')

pdf(file = './res/step_3_fig_210624/macaque_RNA_seurat_QC.pdf',width = 5,height = 7.5)
p1 + p2 + p3 + plot_layout(ncol = 1)
dev.off()

# QC -- correlation between twins
temp <- data.frame(A68A=rowMeans(macaque_RNA_seurat[,macaque_RNA_seurat$donor == 'A68A']@assays$RNA@counts),
                   A68B=rowMeans(macaque_RNA_seurat[,macaque_RNA_seurat$donor == 'A68B']@assays$RNA@counts))
rownames(temp) <- rownames(macaque_RNA_seurat@assays$RNA@counts)
temp <- log1p(temp)

summary(lm(A68B ~ A68A,temp))

p1 <- ggplot(data = temp,aes(x=A68A,y=A68B)) + 
  geom_point(size = 0.1,alpha = 0.8,color = '#440154FF') + 
  xlim(c(0,2)) + ylim(c(0,2)) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        title = element_text(family="Times",face = 'italic',size = 12)) + 
  xlab('log(A68A+1)') + ylab('log(A68B+1)') + labs(title = TeX('R^{2} = 0.9867'))



temp <- data.frame(A84B=rowMeans(macaque_RNA_seurat[,macaque_RNA_seurat$donor == 'A84B']@assays$RNA@counts),
                   A84C=rowMeans(macaque_RNA_seurat[,macaque_RNA_seurat$donor == 'A84C']@assays$RNA@counts))
rownames(temp) <- rownames(macaque_RNA_seurat@assays$RNA@counts)
temp <- log1p(temp)

summary(lm(A84B ~ A84C,temp))

p2 <- ggplot(data = temp,aes(x=A84B,y=A84C)) + 
  geom_point(size = 0.1,alpha = 0.8,color = '#440154FF') + 
  xlim(c(0,2)) + ylim(c(0,2)) +
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        title = element_text(family="Times",face = 'italic',size = 12)) + 
  xlab('log(A84B+1)') + ylab('log(A84C+1)') + labs(title = TeX('R^{2} = 0.9908'))

pdf(file = './res/step_3_fig_210624/macaque_RNA_seurat_replicate_cor.pdf',width = 6,height = 3)
p1 + p2 + plot_layout(nrow = 1)
dev.off()

#layer marker express
macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

layer_marker <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')
show_col(ggsci::pal_futurama()(10))

# #CP
# CP_marker <- c(layer_marker[['CPo']],layer_marker[['CPi']])
# CP_marker <- unique(CP_marker)
# table(CP_marker %in% macaque_annotation$gene_id)
# CP_marker <- CP_marker[CP_marker %in% macaque_annotation$gene_id]
# CP_marker <- macaque_annotation$gene_name[base::match(CP_marker,table = macaque_annotation$gene_id)]
# table(CP_marker %in% rownames(macaque_RNA_seurat))
# CP_marker <- CP_marker[CP_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(CP_marker))
# CP_marker <- unique(CP_marker)
# p1 <- geneSetAveragePlot(genes = CP_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) + 
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts') + 
#   theme_classic() + 
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) + 
#   labs(title = 'CP')
# 
# #SP
# SP_marker <- layer_marker[['SP']]
# SP_marker <- unique(SP_marker)
# table(SP_marker %in% macaque_annotation$gene_id)
# SP_marker <- SP_marker[SP_marker %in% macaque_annotation$gene_id]
# SP_marker <- macaque_annotation$gene_name[base::match(SP_marker,table = macaque_annotation$gene_id)]
# table(SP_marker %in% rownames(macaque_RNA_seurat))
# SP_marker <- SP_marker[SP_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(SP_marker))
# SP_marker <- unique(SP_marker)
# p2 <- geneSetAveragePlot(genes = SP_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) + 
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts') + 
#   theme_classic() + 
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) + 
#   labs(title = 'SP')
# 
# #IZ
# IZ_marker <- layer_marker[['IZ']]
# IZ_marker <- unique(IZ_marker)
# table(IZ_marker %in% macaque_annotation$gene_id)
# IZ_marker <- IZ_marker[IZ_marker %in% macaque_annotation$gene_id]
# IZ_marker <- macaque_annotation$gene_name[base::match(IZ_marker,table = macaque_annotation$gene_id)]
# table(IZ_marker %in% rownames(macaque_RNA_seurat))
# IZ_marker <- IZ_marker[IZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(IZ_marker))
# IZ_marker <- unique(IZ_marker)
# p3 <- geneSetAveragePlot(genes = IZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) + 
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts') + 
#   theme_classic() + 
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) + 
#   labs(title = 'IZ')
# 
# #OSVZ
# OSVZ_marker <- layer_marker[['OSVZ']]
# OSVZ_marker <- unique(OSVZ_marker)
# table(OSVZ_marker %in% macaque_annotation$gene_id)
# OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% macaque_annotation$gene_id]
# OSVZ_marker <- macaque_annotation$gene_name[base::match(OSVZ_marker,table = macaque_annotation$gene_id)]
# table(OSVZ_marker %in% rownames(macaque_RNA_seurat))
# OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(OSVZ_marker))
# OSVZ_marker <- unique(OSVZ_marker)
# p4 <- geneSetAveragePlot(genes = OSVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) + 
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts') + 
#   theme_classic() + 
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) + 
#   labs(title = 'OSVZ')
# 
# #ISVZ
# ISVZ_marker <- layer_marker[['ISVZ']]
# ISVZ_marker <- unique(ISVZ_marker)
# table(ISVZ_marker %in% macaque_annotation$gene_id)
# ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% macaque_annotation$gene_id]
# ISVZ_marker <- macaque_annotation$gene_name[base::match(ISVZ_marker,table = macaque_annotation$gene_id)]
# table(ISVZ_marker %in% rownames(macaque_RNA_seurat))
# ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(ISVZ_marker))
# ISVZ_marker <- unique(ISVZ_marker)
# p5 <- geneSetAveragePlot(genes = ISVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) + 
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts') + 
#   theme_classic() + 
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) + 
#   labs(title = 'ISVZ')
# 
# #VZ
# VZ_marker <- layer_marker[['VZ']]
# VZ_marker <- unique(VZ_marker)
# table(VZ_marker %in% macaque_annotation$gene_id)
# VZ_marker <- VZ_marker[VZ_marker %in% macaque_annotation$gene_id]
# VZ_marker <- macaque_annotation$gene_name[base::match(VZ_marker,table = macaque_annotation$gene_id)]
# table(VZ_marker %in% rownames(macaque_RNA_seurat))
# VZ_marker <- VZ_marker[VZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(VZ_marker))
# VZ_marker <- unique(VZ_marker)
# p6 <- geneSetAveragePlot(genes = VZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) + 
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts') + 
#   theme_classic() + 
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) + 
#   labs(title = 'VZ')
# 
# pdf(file = './res/step_3_fig_210624/layer_marker_enrich_macaque_RNA_seurat_no_limit.pdf',width = 12,height = 5)
# p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
# dev.off()

# #CP
# CP_marker <- c(layer_marker[['CPo']],layer_marker[['CPi']])
# CP_marker <- unique(CP_marker)
# table(CP_marker %in% macaque_annotation$gene_id)
# CP_marker <- CP_marker[CP_marker %in% macaque_annotation$gene_id]
# CP_marker <- macaque_annotation$gene_name[base::match(CP_marker,table = macaque_annotation$gene_id)]
# table(CP_marker %in% rownames(macaque_RNA_seurat))
# CP_marker <- CP_marker[CP_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(CP_marker))
# CP_marker <- unique(CP_marker)
# p1 <- geneSetAveragePlot(genes = CP_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts',limit = c(-2,2)) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) +
#   labs(title = 'CP')
# 
# #SP
# SP_marker <- layer_marker[['SP']]
# SP_marker <- unique(SP_marker)
# table(SP_marker %in% macaque_annotation$gene_id)
# SP_marker <- SP_marker[SP_marker %in% macaque_annotation$gene_id]
# SP_marker <- macaque_annotation$gene_name[base::match(SP_marker,table = macaque_annotation$gene_id)]
# table(SP_marker %in% rownames(macaque_RNA_seurat))
# SP_marker <- SP_marker[SP_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(SP_marker))
# SP_marker <- unique(SP_marker)
# p2 <- geneSetAveragePlot(genes = SP_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts',limit = c(-2,2)) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) +
#   labs(title = 'SP')
# 
# #IZ
# IZ_marker <- layer_marker[['IZ']]
# IZ_marker <- unique(IZ_marker)
# table(IZ_marker %in% macaque_annotation$gene_id)
# IZ_marker <- IZ_marker[IZ_marker %in% macaque_annotation$gene_id]
# IZ_marker <- macaque_annotation$gene_name[base::match(IZ_marker,table = macaque_annotation$gene_id)]
# table(IZ_marker %in% rownames(macaque_RNA_seurat))
# IZ_marker <- IZ_marker[IZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(IZ_marker))
# IZ_marker <- unique(IZ_marker)
# p3 <- geneSetAveragePlot(genes = IZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts',limit = c(-2,2)) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) +
#   labs(title = 'IZ')
# 
# #OSVZ
# OSVZ_marker <- layer_marker[['OSVZ']]
# OSVZ_marker <- unique(OSVZ_marker)
# table(OSVZ_marker %in% macaque_annotation$gene_id)
# OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% macaque_annotation$gene_id]
# OSVZ_marker <- macaque_annotation$gene_name[base::match(OSVZ_marker,table = macaque_annotation$gene_id)]
# table(OSVZ_marker %in% rownames(macaque_RNA_seurat))
# OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(OSVZ_marker))
# OSVZ_marker <- unique(OSVZ_marker)
# p4 <- geneSetAveragePlot(genes = OSVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts',limit = c(-2,2)) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) +
#   labs(title = 'OSVZ')
# 
# #ISVZ
# ISVZ_marker <- layer_marker[['ISVZ']]
# ISVZ_marker <- unique(ISVZ_marker)
# table(ISVZ_marker %in% macaque_annotation$gene_id)
# ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% macaque_annotation$gene_id]
# ISVZ_marker <- macaque_annotation$gene_name[base::match(ISVZ_marker,table = macaque_annotation$gene_id)]
# table(ISVZ_marker %in% rownames(macaque_RNA_seurat))
# ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(ISVZ_marker))
# ISVZ_marker <- unique(ISVZ_marker)
# p5 <- geneSetAveragePlot(genes = ISVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts',limit = c(-2,2)) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) +
#   labs(title = 'ISVZ')
# 
# #VZ
# VZ_marker <- layer_marker[['VZ']]
# VZ_marker <- unique(VZ_marker)
# table(VZ_marker %in% macaque_annotation$gene_id)
# VZ_marker <- VZ_marker[VZ_marker %in% macaque_annotation$gene_id]
# VZ_marker <- macaque_annotation$gene_name[base::match(VZ_marker,table = macaque_annotation$gene_id)]
# table(VZ_marker %in% rownames(macaque_RNA_seurat))
# VZ_marker <- VZ_marker[VZ_marker %in% rownames(macaque_RNA_seurat)]
# table(duplicated(VZ_marker))
# VZ_marker <- unique(VZ_marker)
# p6 <- geneSetAveragePlot(genes = VZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
#                          object.class = 'matrix',plot.type = 'average',
#                          projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
#                          reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
#   scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'Average counts',limit = c(-2,2)) +
#   theme_classic() +
#   theme(axis.ticks = element_blank(),
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.line = element_blank(),
#         legend.position = 'right',
#         panel.border = element_rect(fill = NA,colour = 'black',size = 1.5),
#         aspect.ratio = 1,
#         title = element_text(family="Times",face = 'bold',size = 10)) +
#   labs(title = 'VZ')
# 
# pdf(file = './res/step_3_fig_210624/layer_marker_enrich_macaque_RNA_seurat_limit_2.pdf',width = 12,height = 5)
# p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
# dev.off()

#CP
CP_marker <- c(layer_marker[['CPo']],layer_marker[['CPi']])
CP_marker <- unique(CP_marker)
table(CP_marker %in% macaque_annotation$gene_id)
CP_marker <- CP_marker[CP_marker %in% macaque_annotation$gene_id]
CP_marker <- macaque_annotation$gene_name[base::match(CP_marker,table = macaque_annotation$gene_id)]
table(CP_marker %in% rownames(macaque_RNA_seurat))
CP_marker <- CP_marker[CP_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(CP_marker))
CP_marker <- unique(CP_marker)
p1 <- geneSetAveragePlot(genes = CP_marker,object = macaque_RNA_seurat@assays$RNA@counts,
                         object.class = 'matrix',plot.type = 'average',
                         projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                         reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
  scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'z score',limit = c(-1.5,1.5)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        aspect.ratio = 1,
        title = element_text(family="Times",face = 'bold',size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(face = 'bold',size = 10)) +
  labs(title = 'CP')

#SP
SP_marker <- layer_marker[['SP']]
SP_marker <- unique(SP_marker)
table(SP_marker %in% macaque_annotation$gene_id)
SP_marker <- SP_marker[SP_marker %in% macaque_annotation$gene_id]
SP_marker <- macaque_annotation$gene_name[base::match(SP_marker,table = macaque_annotation$gene_id)]
table(SP_marker %in% rownames(macaque_RNA_seurat))
SP_marker <- SP_marker[SP_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(SP_marker))
SP_marker <- unique(SP_marker)
p2 <- geneSetAveragePlot(genes = SP_marker,object = macaque_RNA_seurat@assays$RNA@counts,
                         object.class = 'matrix',plot.type = 'average',
                         projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                         reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
  scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'z score',limit = c(-1.5,1.5)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        aspect.ratio = 1,
        title = element_text(family="Times",face = 'bold',size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(face = 'bold',size = 10)) +
  labs(title = 'SP')

#IZ
IZ_marker <- layer_marker[['IZ']]
IZ_marker <- unique(IZ_marker)
table(IZ_marker %in% macaque_annotation$gene_id)
IZ_marker <- IZ_marker[IZ_marker %in% macaque_annotation$gene_id]
IZ_marker <- macaque_annotation$gene_name[base::match(IZ_marker,table = macaque_annotation$gene_id)]
table(IZ_marker %in% rownames(macaque_RNA_seurat))
IZ_marker <- IZ_marker[IZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(IZ_marker))
IZ_marker <- unique(IZ_marker)
p3 <- geneSetAveragePlot(genes = IZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
                         object.class = 'matrix',plot.type = 'average',
                         projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                         reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
  scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'z score',limit = c(-1.5,1.5)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        aspect.ratio = 1,
        title = element_text(family="Times",face = 'bold',size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(face = 'bold',size = 10)) +
  labs(title = 'IZ')

#OSVZ
OSVZ_marker <- layer_marker[['OSVZ']]
OSVZ_marker <- unique(OSVZ_marker)
table(OSVZ_marker %in% macaque_annotation$gene_id)
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% macaque_annotation$gene_id]
OSVZ_marker <- macaque_annotation$gene_name[base::match(OSVZ_marker,table = macaque_annotation$gene_id)]
table(OSVZ_marker %in% rownames(macaque_RNA_seurat))
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(OSVZ_marker))
OSVZ_marker <- unique(OSVZ_marker)
p4 <- geneSetAveragePlot(genes = OSVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
                         object.class = 'matrix',plot.type = 'average',
                         projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                         reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
  scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'z score',limit = c(-1.5,1.5)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        aspect.ratio = 1,
        title = element_text(family="Times",face = 'bold',size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(face = 'bold',size = 10)) +
  labs(title = 'OSVZ')

#ISVZ
ISVZ_marker <- layer_marker[['ISVZ']]
ISVZ_marker <- unique(ISVZ_marker)
table(ISVZ_marker %in% macaque_annotation$gene_id)
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% macaque_annotation$gene_id]
ISVZ_marker <- macaque_annotation$gene_name[base::match(ISVZ_marker,table = macaque_annotation$gene_id)]
table(ISVZ_marker %in% rownames(macaque_RNA_seurat))
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(ISVZ_marker))
ISVZ_marker <- unique(ISVZ_marker)
p5 <- geneSetAveragePlot(genes = ISVZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
                         object.class = 'matrix',plot.type = 'average',
                         projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                         reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
  scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'z score',limit = c(-1.5,1.5)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        aspect.ratio = 1,
        title = element_text(family="Times",face = 'bold',size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(face = 'bold',size = 10)) +
  labs(title = 'ISVZ')

#VZ
VZ_marker <- layer_marker[['VZ']]
VZ_marker <- unique(VZ_marker)
table(VZ_marker %in% macaque_annotation$gene_id)
VZ_marker <- VZ_marker[VZ_marker %in% macaque_annotation$gene_id]
VZ_marker <- macaque_annotation$gene_name[base::match(VZ_marker,table = macaque_annotation$gene_id)]
table(VZ_marker %in% rownames(macaque_RNA_seurat))
VZ_marker <- VZ_marker[VZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(VZ_marker))
VZ_marker <- unique(VZ_marker)
p6 <- geneSetAveragePlot(genes = VZ_marker,object = macaque_RNA_seurat@assays$RNA@counts,
                         object.class = 'matrix',plot.type = 'average',
                         projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                         reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
  scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'z score',limit = c(-1.5,1.5)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        aspect.ratio = 1,
        title = element_text(family="Times",face = 'bold',size = 10),
        legend.text = element_text(size = 8),
        legend.title = element_text(face = 'bold',size = 10)) +
  labs(title = 'VZ')

pdf(file = './res/step_3_fig_210624/layer_marker_enrich_macaque_RNA_seurat_limit_1.5.pdf',width = 7.5,height = 5)
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
dev.off()

pdf(file = './res/step_3_fig_210624/layer_marker_enrich_macaque_RNA_seurat_limit_1.5_legend.pdf',width = 7.5,height = 5)
geneSetAveragePlot(genes = SP_marker,object = macaque_RNA_seurat@assays$RNA@counts,
                   object.class = 'matrix',plot.type = 'average',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',color.palette = c("grey","#440154FF","#FDE725FF"),print = FALSE,scaled = TRUE) +
  scale_color_gradientn(colours = c("#008EA0FF","grey","#C71000FF"),name = 'z score',limit = c(-1.5,1.5)) +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'right',
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        aspect.ratio = 1,
        title = element_text(family="Times",face = 'bold',size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(face = 'bold',size = 12)) +
  labs(title = 'SP')
dev.off()

#sub cluster
#phase -- cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
s.genes[!(s.genes %in% rownames(macaque_RNA_seurat@assays$RNA@counts))]
s.genes <- c(s.genes[s.genes %in% rownames(macaque_RNA_seurat)],'ENSMMUG00000004522','CENPU')
table((s.genes %in% rownames(macaque_RNA_seurat)))
table((g2m.genes %in% rownames(macaque_RNA_seurat)))
g2m.genes[!(g2m.genes %in% rownames(macaque_RNA_seurat@assays$RNA@counts))]
g2m.genes <- c(g2m.genes[g2m.genes %in% rownames(macaque_RNA_seurat)],'PIMREG','JPT1')
table((g2m.genes %in% rownames(macaque_RNA_seurat@assays$RNA@counts)))
macaque_RNA_seurat <- CellCycleScoring(macaque_RNA_seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)

show_col(ggsci::pal_aaas()(10))
# color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
# temp <- data.frame(id = c('G1','G2M','S'),
#                    col = ggsci::pal_aaas()(10)[1:3])
# color_paramater <- rbind(color_paramater,temp)
# 
# write.csv(color_paramater,file = './data/parameter/color_paramater.csv')
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)

temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

pdf(file = './res/step_3_fig_210624/macaque_RNA_seurat_phase_dimplot.pdf',width = 5,height = 5)
DimPlot(macaque_RNA_seurat,group.by = 'Phase',cols = temp_col[c('G1','G2M','S')]) + 
  theme_classic() + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = c(0.1,0.9),
        aspect.ratio = 1,
        legend.text = element_text(size = 10),
        legend.title = element_text(face = 'bold',size = 12)) + 
  labs(title = '',colour = 'Phase')
dev.off()

#phase proportion
phase_proportion <- My_Cell_Proportion(macaque_RNA_seurat,group.by = 'Phase',split.by = 'cell_type')
phase_proportion$cell_type <- factor(phase_proportion$cell_type,levels = c('End','Per','Mic','Astrocyte','vRG','oRG','OPC','Cycling','IP','Ex','Ex-U','Ex-SP','InCGE','InMGE','InPSB'))
pdf(file = './res/step_3_fig_210624/phase_cell_proportion_barplot_all_cell_type.pdf',width = 3,height = 5)
ggplot(data = phase_proportion, mapping = aes(x = cell_type,y = Proportion,fill = Phase))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.95) + 
  scale_fill_manual(values = temp_col[c('G1','G2M','S')]) + 
  theme_classic() + 
  theme(aspect.ratio = 2.5,
        axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = 'none') + 
  labs(fill = 'Phase') + xlab('') + ylab('Percent of Cells') + 
  coord_flip()
dev.off()


#sub cluster divided by phase
macaque_RNA_seurat$sub_cell_type <- macaque_RNA_seurat$cell_type
macaque_RNA_seurat@meta.data[,"sub_cell_type"] <- as.character(macaque_RNA_seurat@meta.data[,"sub_cell_type"])
DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type == 'oRG',"sub_cell_type"] <- 'oRG'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters == '26',"sub_cell_type"] <- 'vRG-1'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$seurat_clusters == '23',"sub_cell_type"] <- 'vRG-2'

temp <- macaque_RNA_seurat
Idents(temp) <- 'cell_type'
temp <- FindSubCluster(temp,cluster = 'Cycling',graph.name = 'SCT_snn',resolution = 0.1)
DimPlot(temp,group.by = 'sub.cluster',label = TRUE,repel = TRUE)

macaque_RNA_seurat@meta.data[colnames(temp)[temp$sub.cluster == 'Cycling_0'],"sub_cell_type"] <- 'Cyc-G2M'
macaque_RNA_seurat@meta.data[colnames(temp)[temp$sub.cluster == 'Cycling_1'],"sub_cell_type"] <- 'Cyc-S'

#dimplot for the progenitor cell
temp <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('oRG','vRG','Cycling')]
temp <- temp[,(temp@reductions$umap@cell.embeddings[,"UMAP_1"] > 4 & temp@reductions$umap@cell.embeddings[,"UMAP_2"] > 3)]
DimPlot(temp,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)

show_col(c('#4DBBD5FF','#3C5488FF'))
show_col(colorRampPalette(colors = c('#4DBBD5FF','#3C5488FF'))(6))

# color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
# temp <- data.frame(id = c('vRG-1','vRG-2','Cyc-G2M','Cyc-S'),
#                    col = colorRampPalette(colors = c('#4DBBD5FF','#3C5488FF'))(6)[2:5])
# color_paramater <- rbind(color_paramater,temp)
# write.csv(color_paramater,file = './data/parameter/color_paramater.csv')

color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)

temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

pdf(file = './res/step_3_fig_210624/progenitor_UMAP_divided_by_Phase.pdf',width = 3,height = 3)
DimPlot(temp,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,cols = temp_col[c('oRG','vRG-1','vRG-2','Cyc-G2M','Cyc-S')]) + 
  theme_classic() + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none',
        aspect.ratio = 1) + 
  labs(title = '')
dev.off()

#heatmap for progenitor cell
show_col(ggsci::pal_gsea()(10))

temp <- MyPseudoBluk(macaque_RNA_seurat,meta.var = 'sub_cell_type')
g2m_signature <- cc.genes$g2m.genes[cc.genes$g2m.genes %in% rownames(macaque_RNA_seurat)]
s_signature <- cc.genes$s.genes[cc.genes$s.genes %in% rownames(macaque_RNA_seurat)]

temp <- temp[,c('oRG','vRG-1','vRG-2','Cyc-G2M','Cyc-S')]
temp <- temp[c('SOX9','HES1','FBXO32','CNN2','MOXD1','HOPX','FAM107A','MT3',g2m_signature,s_signature),]
temp <- rbind(temp,G2M_sig=colMeans(temp[g2m_signature,]))
temp <- rbind(temp,S_sig=colMeans(temp[s_signature,]))
temp <- temp[c('SOX9','HES1','FBXO32','CNN2','MOXD1','HOPX','FAM107A','MT3','S_sig','G2M_sig'),]
temp <- log1p(temp)
temp <- t(scale(t(temp),center = TRUE))

col_dend <- as.dendrogram(hclust(dist(t(temp)),method = 'ward.D'))
col_dend_reorder <- reorder(col_dend,1:5)

pdf(file = './res/step_3_fig_210624/progenitor_heatmap_divided_by_Phase.pdf',width = 3,height = 3)
ComplexHeatmap::Heatmap(temp,show_row_names = TRUE,show_column_names = TRUE,
                        cluster_rows = FALSE,cluster_columns = col_dend_reorder,
                        rect_gp = gpar(col = 'black'),height = unit(2,units = 'inches'),width = unit(1,units = 'inches'),
                        heatmap_legend_param = list(title = 'z score'),
                        row_names_gp = gpar(fontsize = 10),column_names_gp = gpar(fontsize = 10),
                        col = c('#4500ACFF','white','#FF5959FF'),
                        column_dend_height = unit(0.2,units = 'inches'))
dev.off()

#save macaque_RNA_seurat
saveRDS(macaque_RNA_seurat,file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')


#featureplot for some marker
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
macaque_RNA_seurat <- ScaleData(macaque_RNA_seurat,features = rownames(macaque_RNA_seurat))

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_1.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('TOP2A','MKI67','HES1','ATP1A2','MOXD1','HOPX','FAM107A'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_2.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('NPY','CD9','GPX3','TNC','CRYAB','NR4A1'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_3.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('SOX10','MBP','ENSMMUG00000056728','OLIG2','AQP4','APOE','AGT'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_4.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('PPP1R17','EOMES','PENK','NEUROG1','NEUROG2'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_5.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('PPP1R17','EOMES','PENK','NEUROG1','NEUROG2'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_6.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('TBR1','BCL11B','SATB2','SLC17A7','NEUROD2','NEUROD6'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_7.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('NR4A2','CRYM','ST18','CDH18','DLX5','GAD2'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_8.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('LHX6','SST','SP8','NR2F2','MEIS2','PAX6','ETV1'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_3_fig_210624/fig_s1_featureplot_9.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('CX3CR1','AIF1','CCL3','CLDN5','FOXC2','PDGFRB'),object = macaque_RNA_seurat@assays$RNA@data,
                   object.class = 'matrix',plot.type = 'panels',
                   projection = macaque_RNA_seurat@reductions$umap@cell.embeddings,
                   reduction_key = 'UMAP',print = FALSE,scaled = FALSE,num.panel.rows = 1) + 
  scale_color_gradientn(colours = c("grey","#440154FF","#FDE725FF"),name = 'express') + 
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        legend.position = 'none',
        panel.border = element_blank(),
        aspect.ratio = 1,
        strip.background = element_blank()) + 
  labs(title = '')
dev.off()