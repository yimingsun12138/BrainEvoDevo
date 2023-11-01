#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: integrate macaque all tech RNA_seq data                         ##
## Data: 2022.03.05                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# scRNA and multiome integration using harmony ----------------------------
#load data
cell_list <- readRDS(file = './ArchR/res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220306_1630.rds')
macaque_integration_Seurat <- readRDS(file = './processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_annotated_220302.rds')
macaque_RNA_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'scRNA']
macaque_multiome_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']

macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_seurat) %in% cell_list]
macaque_multiome_seurat <- macaque_multiome_seurat[,colnames(macaque_multiome_seurat) %in% cell_list]

# #separate analysis
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
# macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 25,resolution = 1.5,group.by = 'old_cell_type',label = TRUE)
# 
# macaque_multiome_seurat <- my_process_seurat(object = macaque_multiome_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
# macaque_multiome_seurat <- my_process_seurat(object = macaque_multiome_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 25,resolution = 1.5,group.by = 'old_cell_type',label = TRUE)

#harmony
DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

macaque_RNA_seurat <- FindVariableFeatures(object = macaque_RNA_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- VariableFeatures(macaque_RNA_seurat)
table(gene_list %in% rownames(macaque_multiome_seurat))

dim_used <- 28

macaque_integration_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_seurat),
                                                     assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                           multiome=c('nCount_RNA','donor')),
                                                     npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                     UMAP_dim = dim_used,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                                     yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_seurat),"cell_type"] <- macaque_multiome_seurat$old_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$old_cell_type

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#annotate
p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

macaque_integration_Seurat$old_cell_type <- macaque_integration_Seurat$cell_type
macaque_integration_Seurat$cell_type <- NA

#cluster 15 RG-1
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('15'),"cell_type"] <- 'RG-1'

#cluster 11 Cycling
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('11'),"cell_type"] <- 'Cycling'

#cluster 16 RG-2
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('16'),"cell_type"] <- 'RG-2'

#cluster 23 RG-3
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('23'),"cell_type"] <- 'RG-3'

#cluster 19 31 OPC
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('19','31'),"cell_type"] <- 'OPC'

#cluster 17 IP
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('17'),"cell_type"] <- 'IP'

#cluster 0 5 10 Ex-1
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('0','5','10'),"cell_type"] <- 'Ex-1'

#cluster 1 7 13 6 Ex-2
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('1','7','13','6'),"cell_type"] <- 'Ex-2'

#cluster 18 8 Ex-3
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('18','8'),"cell_type"] <- 'Ex-3'

#cluster 3 9 12 26 InCGE
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('3','9','12','26'),"cell_type"] <- 'InCGE'

#cluster 2 4 14 21 24 25 InMGE
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('2','4','14','21','24','25'),"cell_type"] <- 'InMGE'

#cluster 27 Ependymal
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('27'),"cell_type"] <- 'Ependymal'

#cluster 22 Mic
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('22'),"cell_type"] <- 'Mic'

#cluster 20 End
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('20'),"cell_type"] <- 'End'

#cluster 30 28 29 Per
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$seurat_clusters %in% c('28','29','30'),"cell_type"] <- 'Per'

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#cycling sub-cluster
Idents(macaque_integration_Seurat) <- 'cell_type'
macaque_integration_Seurat <- FindSubCluster(object = macaque_integration_Seurat,cluster = 'Cycling',resolution = 0.1,subcluster.name = 'temp',graph.name = 'integration_snn')
p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'temp',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

macaque_integration_Seurat@meta.data[macaque_integration_Seurat$temp %in% c('Cycling_0'),"cell_type"] <- 'Cyc-S'
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$temp %in% c('Cycling_1'),"cell_type"] <- 'Cyc-G2M'

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#validate
#old cell type
scibet::Confusion_heatmap(ori = macaque_integration_Seurat$old_cell_type,prd = macaque_integration_Seurat$cell_type)

#marker gene
macaque_integration_Seurat <- NormalizeData(object = macaque_integration_Seurat)
dotplot_matrix <- my_dotplot(macaque_integration_Seurat,assay = 'integration', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('AQP4','FOXJ1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))
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

meta_data <- macaque_integration_Seurat@meta.data
meta_data <- meta_data[,c('cell_id','cell_type')]

#save data
saveRDS(meta_data,file = './res/step_17_fig_220305/macaque_integration_Seurat_meta_data_220306_2326.rds')
saveRDS(macaque_integration_Seurat@reductions$pca,file = './res/step_17_fig_220305/macaque_integration_Seurat_harmony_pca_220306_2326.rds')
saveRDS(macaque_integration_Seurat@reductions$umap,file = './res/step_17_fig_220305/macaque_integration_Seurat_harmony_umap_220306_2326.rds')

# scRNA and multiome merge ------------------------------------------------
macaque_integration_Seurat <- readRDS(file = './processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_annotated_220302.rds')
meta_data <- readRDS(file = './res/step_17_fig_220305/macaque_integration_Seurat_meta_data_220306_2326.rds')
temp <- macaque_integration_Seurat@meta.data

macaque_integration_Seurat <- macaque_integration_Seurat@assays$RNA@counts
macaque_integration_Seurat <- macaque_integration_Seurat[,rownames(meta_data)]
macaque_integration_Seurat <- CreateSeuratObject(counts = macaque_integration_Seurat,project = 'macaque',assay = 'RNA',min.cells = 0,min.features = 0)

macaque_integration_Seurat$batch <- temp[colnames(macaque_integration_Seurat),'batch']
macaque_integration_Seurat$donor <- temp[colnames(macaque_integration_Seurat),'donor']
macaque_integration_Seurat$tech <- temp[colnames(macaque_integration_Seurat),'tech']
macaque_integration_Seurat$cell_type <- meta_data[colnames(macaque_integration_Seurat),'cell_type']

#process
macaque_integration_Seurat <- my_process_seurat(object = macaque_integration_Seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','batch','donor','tech'),npcs = 50,preprocess = TRUE)
macaque_integration_Seurat <- my_process_seurat(object = macaque_integration_Seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 28,resolution = c(0.5,1,1.5),group.by = 'cell_type',label = TRUE)

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'tech',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#add dimreduction object
temp <- readRDS(file = './res/step_17_fig_220305/macaque_integration_Seurat_harmony_pca_220306_2326.rds')
macaque_integration_Seurat@reductions$harmonyPCA <- temp
temp <- readRDS(file = './res/step_17_fig_220305/macaque_integration_Seurat_harmony_umap_220306_2326.rds')
macaque_integration_Seurat@reductions$harmonyUMAP <- temp

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'tech',label = FALSE,reduction = 'harmonyUMAP')
p1+p2+p3+plot_layout(ncol = 3)

#dimplot
temp <- read.csv(file = './data/parameter/color_paramater.csv')
p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP',cols = temp$col[1:16]) + theme(aspect.ratio = 1)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'tech',label = FALSE,reduction = 'harmonyUMAP') + theme(aspect.ratio = 1)

pdf(file = './res/step_17_fig_220305/macaque_integration_Seurat_cell_type_dimplot.pdf',width = 16,height = 7)
p1+p2+plot_layout(ncol = 2)
dev.off()

#marker dotplot
dotplot_matrix <- my_dotplot(macaque_integration_Seurat,assay = 'RNA', 
                             col.max = 2, col.min = -2, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Epe=c('AQP4','FOXJ1'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))

pdf(file = './res/step_17_fig_220305/macaque_integration_Seurat_cell_type_marker_gene_dotplot.pdf',width = 13,height = 6)
my_dotplot(data_plot = dotplot_matrix,col.max = 2, col.min = -2,
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

#save data
saveRDS(macaque_integration_Seurat,file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
