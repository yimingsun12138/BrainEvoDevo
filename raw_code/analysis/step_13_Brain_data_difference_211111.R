#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Brain data difference                                           ##
## Data: 2021.11.11                                                                ##
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
library(SummarizedExperiment)
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
library(scales)
library(ggsci)
library(VennDiagram)
library(harmony)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# mouse unknown_Ex --------------------------------------------------------
#load data
mouse_RNA_seurat <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/scRNA_raw_data/mouse_RNA_seurat_E14_to_E18_211111.rds')
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
Brain_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_mouse_RNA_harmony_MNN_seurat_211111.rds')

DimPlot(mouse_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'umap')

#umap
pdf(file = './res/step_13_fig_211111/macaque_RNA_200919_210922_mouse_RNA_harmony_mnn_umap.pdf',width = 16,height = 8)
DimPlot(object = Brain_RNA_seurat,group.by = 'predict_label',label = TRUE,repel = TRUE,split.by = 'dataset')
dev.off()

#marker dotpot
mouse_RNA_seurat <- my_process_seurat(object = mouse_RNA_seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

dotplot_matrix <- my_dotplot(mouse_RNA_seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('Cldn5','Pecam1'),
                                             Per=c('Pdgfrb'),
                                             Mic=c('Cx3cr1'),
                                             Ast=c('Aqp4','Apoe'),
                                             RG=c('Sox2','Pax6','Hes5'),
                                             OPC=c('Sox10','Olig2','Egfr'),
                                             Cyc=c('Top2a','Mki67','Clspn','Aurka'),
                                             IP=c('Eomes','Ppp1r17','Neurog2'),
                                             Ex=c('Neurod2','Neurod6','Tbr1','Satb2','Slc17a7','Fezf2'),
                                             In=c('Dlx5','Gad2','Gad1','Dlx2'),
                                             MGE=c('Lhx6','Sst'),
                                             CGE=c('Sp8','Cxcl14','Htr3a'),
                                             PSB=c('Meis2','Etv1')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('low_quality','InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_13_fig_211111/mouse_RNA_seurat_mnn_label_marker_dotplot.pdf',width = 16,height = 8)
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
dev.off()

#Identify the unknown-Ex
pdf(file = './res/step_13_fig_211111/mouse_RNA_seurat_mnn_label_marker_vlnplot.pdf',width = 24,height = 14)
VlnPlot(mouse_RNA_seurat,features = c('Fezf2','Bcl11b','Tle4','Foxp2','Sox5','Ldb2','Crym','Pcp4','Nr4a2','Tshz2',
                                      'Satb2','Ptn','Cux1','Cux2','Rorb','Lmo4','Lpl','Ldb2','Reln','Lhx5'),group.by = 'cell_type',pt.size = 0)
dev.off()

#marker of unknown_Ex
unknown_Ex_marker <- FindMarkers(object = mouse_RNA_seurat,ident.1 = 'low_quality',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE)

pdf(file = './res/step_13_fig_211111/mouse_RNA_seurat_mnn_label_unknown_Ex_marker_vlnplot.pdf',width = 18,height = 9)
VlnPlot(mouse_RNA_seurat,features = rownames(unknown_Ex_marker)[1:12],pt.size = 0,group.by = 'cell_type',assay = 'RNA',slot = 'data')
dev.off()

# cell type proportion ----------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
greenleaf_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
mouse_RNA_seurat <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/scRNA_raw_data/mouse_RNA_seurat_E14_to_E18_211111.rds')

#modify annotation
table(macaque_RNA_seurat$cell_type)
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('End','Per','Astrocyte','InPSB'))]
macaque_RNA_seurat$basic_cell_type <- as.character(macaque_RNA_seurat$cell_type)
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('RG-1','RG-2','RG-3'),"basic_cell_type"] <- 'RG'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('Ex-1','Ex-2'),"basic_cell_type"] <- 'ExN'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('Ex-3'),"basic_cell_type"] <- 'ExUp'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('Ex-4'),"basic_cell_type"] <- 'ExDp'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$cell_type %in% c('Cyc-G2M','Cyc-S'),"basic_cell_type"] <- 'Cyc'
DimPlot(macaque_RNA_seurat,group.by = 'basic_cell_type',label = TRUE,repel = TRUE)

temp <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_greenleaf_RNA_harmony_MNN_seurat_211102.rds')
temp <- temp[,temp$dataset == 'human']
greenleaf_RNA_seurat <- greenleaf_RNA_seurat[,colnames(temp)]
greenleaf_RNA_seurat$basic_cell_type <- as.character(temp$predict_label)
greenleaf_RNA_seurat <- greenleaf_RNA_seurat[,!(greenleaf_RNA_seurat$basic_cell_type %in% c('End','Per','Astrocyte','InPSB','Unknown'))]
greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$basic_cell_type %in% c('RG-1','RG-2','RG-3'),"basic_cell_type"] <- 'RG'
greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$basic_cell_type %in% c('Ex-1','Ex-2'),"basic_cell_type"] <- 'ExN'
greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$basic_cell_type %in% c('Ex-3'),"basic_cell_type"] <- 'ExUp'
greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$basic_cell_type %in% c('Ex-4'),"basic_cell_type"] <- 'ExDp'
greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$basic_cell_type %in% c('Cyc-G2M','Cyc-S'),"basic_cell_type"] <- 'Cyc'
DimPlot(greenleaf_RNA_seurat,group.by = 'basic_cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP')
rm(temp)
gc()

temp <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_PDhuman_RNA_harmony_MNN_seurat_211029.rds')
temp <- temp[,temp$dataset == 'human']
PDhuman_RNA_seurat <- PDhuman_RNA_seurat[,colnames(temp)]
PDhuman_RNA_seurat$basic_cell_type <- as.character(temp$predict_label)
PDhuman_RNA_seurat <- PDhuman_RNA_seurat[,!(PDhuman_RNA_seurat$basic_cell_type %in% c('End','Per','Astrocyte','InPSB','Unknown'))]
PDhuman_RNA_seurat@meta.data[PDhuman_RNA_seurat$basic_cell_type %in% c('RG-1','RG-2','RG-3'),"basic_cell_type"] <- 'RG'
PDhuman_RNA_seurat@meta.data[PDhuman_RNA_seurat$basic_cell_type %in% c('Ex-1','Ex-2'),"basic_cell_type"] <- 'ExN'
PDhuman_RNA_seurat@meta.data[PDhuman_RNA_seurat$basic_cell_type %in% c('Ex-3'),"basic_cell_type"] <- 'ExUp'
PDhuman_RNA_seurat@meta.data[PDhuman_RNA_seurat$basic_cell_type %in% c('Ex-4'),"basic_cell_type"] <- 'ExDp'
PDhuman_RNA_seurat@meta.data[PDhuman_RNA_seurat$basic_cell_type %in% c('Cyc-G2M','Cyc-S'),"basic_cell_type"] <- 'Cyc'
PDhuman_RNA_seurat@reductions$umap <- temp@reductions$umap
DimPlot(PDhuman_RNA_seurat,group.by = 'basic_cell_type',label = TRUE,repel = TRUE,reduction = 'umap')
rm(temp)
gc()

temp <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_mouse_RNA_harmony_MNN_seurat_211111.rds')
temp <- temp[,temp$dataset == 'mouse']
mouse_RNA_seurat <- mouse_RNA_seurat[,colnames(temp)]
mouse_RNA_seurat$basic_cell_type <- as.character(temp$predict_label)
mouse_RNA_seurat <- mouse_RNA_seurat[,!(mouse_RNA_seurat$basic_cell_type %in% c('End','Per','Astrocyte','InPSB','Unknown','low_quality'))]
mouse_RNA_seurat@meta.data[mouse_RNA_seurat$basic_cell_type %in% c('RG-1','RG-2','RG-3'),"basic_cell_type"] <- 'RG'
mouse_RNA_seurat@meta.data[mouse_RNA_seurat$basic_cell_type %in% c('Ex-1','Ex-2'),"basic_cell_type"] <- 'ExN'
mouse_RNA_seurat@meta.data[mouse_RNA_seurat$basic_cell_type %in% c('Ex-3'),"basic_cell_type"] <- 'ExUp'
mouse_RNA_seurat@meta.data[mouse_RNA_seurat$basic_cell_type %in% c('Ex-4'),"basic_cell_type"] <- 'ExDp'
mouse_RNA_seurat@meta.data[mouse_RNA_seurat$basic_cell_type %in% c('Cyc-G2M','Cyc-S'),"basic_cell_type"] <- 'Cyc'
mouse_RNA_seurat@reductions$umap <- temp@reductions$umap
DimPlot(mouse_RNA_seurat,group.by = 'basic_cell_type',label = TRUE,repel = TRUE,reduction = 'umap')
rm(temp)
gc()

#calculate cell type proportion
macaque_RNA_seurat$dataset <- 'macaque'
greenleaf_RNA_seurat$dataset <- 'greenleaf'
PDhuman_RNA_seurat$dataset <- 'Neuron'
mouse_RNA_seurat$dataset <- 'mouse'
PDhuman_RNA_seurat$Age <- paste(PDhuman_RNA_seurat$Gestation_week,'week',sep = '_')
macaque_RNA_seurat$Age <- NA
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$donor %in% c('A50A','A50B'),"Age"] <- 'E80'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$donor %in% c('A68A','A68B'),"Age"] <- 'E92'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$donor %in% c('A82A','A82B'),"Age"] <- 'E90'
macaque_RNA_seurat@meta.data[macaque_RNA_seurat$donor %in% c('A84B','A84C'),"Age"] <- 'E84'

temp <- My_Cell_Proportion(seu.obj = greenleaf_RNA_seurat,group.by = 'basic_cell_type',split.by = 'Age')
temp$dataset <- 'greenleaf'
cell_proportion_dataframe <- temp

temp <- My_Cell_Proportion(seu.obj = PDhuman_RNA_seurat,group.by = 'basic_cell_type',split.by = 'Age')
temp$dataset <- 'Neuron'
cell_proportion_dataframe <- rbind(cell_proportion_dataframe,temp)

temp <- My_Cell_Proportion(seu.obj = macaque_RNA_seurat,group.by = 'basic_cell_type',split.by = 'Age')
temp$dataset <- 'macaque'
cell_proportion_dataframe <- rbind(cell_proportion_dataframe,temp)

temp <- My_Cell_Proportion(seu.obj = mouse_RNA_seurat,group.by = 'basic_cell_type',split.by = 'Age')
temp$dataset <- 'mouse'
cell_proportion_dataframe <- rbind(cell_proportion_dataframe,temp)

cell_proportion_dataframe$dataset <- factor(cell_proportion_dataframe$dataset,levels = c('greenleaf','Neuron','macaque','mouse'))
cell_proportion_dataframe$basic_cell_type <- factor(cell_proportion_dataframe$basic_cell_type,levels = c('Mic','InCGE','InMGE','unknown_Ex','ExDp','ExUp','ExN','IP','Cyc','OPC','RG'))

pdf(file = './res/step_13_fig_211111/all_species_cell_type_Proportion_dotplot.pdf',width = 6,height = 5)
ggplot(data = cell_proportion_dataframe,aes(x=Age,y=basic_cell_type)) + 
  geom_point(aes(size = Proportion,color=basic_cell_type)) + 
  facet_grid(~dataset,scales = "free_x",space = "free_x",switch = "y") + 
  scale_color_manual(values = c(ggsci::pal_aaas()(10),'grey')) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        panel.grid = element_line(color="grey",size = 0.1))
dev.off()

#normalize the proportion
for (i in unique(cell_proportion_dataframe$basic_cell_type)) {
  temp <- cell_proportion_dataframe[cell_proportion_dataframe$basic_cell_type == i,]
  temp <- temp$Proportion/max(temp$Proportion)
  cell_proportion_dataframe[cell_proportion_dataframe$basic_cell_type == i,"Proportion"] <- temp
}

pdf(file = './res/step_13_fig_211111/all_species_cell_type_Proportion_normalized_dotplot.pdf',width = 6,height = 5)
ggplot(data = cell_proportion_dataframe,aes(x=Age,y=basic_cell_type)) + 
  geom_point(aes(size = Proportion,color=basic_cell_type)) + 
  facet_grid(~dataset,scales = "free_x",space = "free_x",switch = "y") + 
  scale_color_manual(values = c(ggsci::pal_aaas()(10),'grey')) + 
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 10),
        strip.text = element_text(size = 8),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        panel.grid = element_line(color="grey",size = 0.1))
dev.off()

# marker gene intersect ---------------------------------------------------
mouse_to_human_anno <- read.csv(file = './data/reference/GRCm39_to_GRCh38.csv')
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')

temp <- mouse_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
mouse_RNA_seurat[['converted']] <- temp

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp

#preprocess
greenleaf_RNA_seurat <- my_process_seurat(object = greenleaf_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
PDhuman_RNA_seurat <- my_process_seurat(object = PDhuman_RNA_seurat,assay = 'RNA',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
mouse_RNA_seurat <- my_process_seurat(object = mouse_RNA_seurat,assay = 'converted',reduction.name = 'PCA',variable.feature = NULL,nfeatures = 3000,vars.to.regress = NULL,npcs = 50,TRUE)

#find all markers
Idents(greenleaf_RNA_seurat) <- 'basic_cell_type'
greenleaf_marker <- FindAllMarkers(object = greenleaf_RNA_seurat,assay = 'converted',test.use = 'bimod',slot = 'data',only.pos = TRUE)

Idents(PDhuman_RNA_seurat) <- 'basic_cell_type'
PDhuman_marker <- FindAllMarkers(object = PDhuman_RNA_seurat,assay = 'RNA',test.use = 'bimod',slot = 'data',only.pos = TRUE)

Idents(macaque_RNA_seurat) <- 'basic_cell_type'
macaque_marker <- FindAllMarkers(object = macaque_RNA_seurat,assay = 'converted',test.use = 'bimod',slot = 'data',only.pos = TRUE)

Idents(mouse_RNA_seurat) <- 'basic_cell_type'
mouse_marker <- FindAllMarkers(object = mouse_RNA_seurat,assay = 'converted',test.use = 'bimod',slot = 'data',only.pos = TRUE)

my_send_sms(body = 'find all markers done!')

write.csv(greenleaf_marker,file = './res/step_13_fig_211111/greenleaf_marker_211111.csv')
write.csv(PDhuman_marker,file = './res/step_13_fig_211111/PDhuman_marker_211111.csv')
write.csv(macaque_marker,file = './res/step_13_fig_211111/macaque_marker_211111.csv')
write.csv(mouse_marker,file = './res/step_13_fig_211111/mouse_marker_211111.csv')

#marker gene overlap
venn.plot <- venn.diagram(
  x = list(greenleaf=unique(greenleaf_marker$gene),
           Neuron=unique(PDhuman_marker$gene),
           macaque=unique(macaque_marker$gene),
           mouse=unique(mouse_marker$gene)),
  filename = NULL,
  col = "white",
  lty = 1,
  lwd = 4,
  fill = c("#ffb2b2","#b2e7cb","#b2d4ec","#d3c0e2"),
  alpha = 0.50,
  cex = 2.5,
  fontfamily = "serif",
  fontface = "bold",
  cat.col = c("#ffb2b2","#b2e7cb","#b2d4ec","#d3c0e2"),
  cat.cex = 2.5,
  cat.fontfamily = "serif"
)

pdf(file = './res/step_13_fig_211111/all_species_all_marker_intersect_vennplot.pdf',width = 12,height = 10)
grid.draw(venn.plot)
dev.off()
grid.newpage()

#different cell type marker overlap to show the mainly difference in Pg or mature neuron
cell_list <- c('Mic','InCGE','InMGE','ExDp','ExUp','ExN','IP','Cyc','OPC','RG')

overlap_heatmap <- matrix(nrow = length(cell_list),ncol = length(cell_list))
rownames(overlap_heatmap) <- cell_list
colnames(overlap_heatmap) <- cell_list
overlap_heatmap <- as.data.frame(overlap_heatmap)

my_find_intersect <- function(x,y){
  temp_greenleaf <- dplyr::intersect(greenleaf_marker[greenleaf_marker$cluster == x,"gene"],greenleaf_marker[greenleaf_marker$cluster == y,"gene"])
  temp_PDhuman <- dplyr::intersect(PDhuman_marker[PDhuman_marker$cluster == x,"gene"],PDhuman_marker[PDhuman_marker$cluster == y,"gene"])
  temp_macaque <- dplyr::intersect(macaque_marker[macaque_marker$cluster == x,"gene"],macaque_marker[macaque_marker$cluster == y,"gene"])
  temp_mouse <- dplyr::intersect(mouse_marker[mouse_marker$cluster == x,"gene"],mouse_marker[mouse_marker$cluster == y,"gene"])
  temp <- dplyr::intersect(temp_greenleaf,temp_PDhuman)
  temp <- dplyr::intersect(temp,temp_macaque)
  temp <- dplyr::intersect(temp,temp_mouse)
  return(length(temp))
}

my_find_union <- function(x,y){
  temp_greenleaf <- dplyr::union(greenleaf_marker[greenleaf_marker$cluster == x,"gene"],greenleaf_marker[greenleaf_marker$cluster == y,"gene"])
  temp_PDhuman <- dplyr::union(PDhuman_marker[PDhuman_marker$cluster == x,"gene"],PDhuman_marker[PDhuman_marker$cluster == y,"gene"])
  temp_macaque <- dplyr::union(macaque_marker[macaque_marker$cluster == x,"gene"],macaque_marker[macaque_marker$cluster == y,"gene"])
  temp_mouse <- dplyr::union(mouse_marker[mouse_marker$cluster == x,"gene"],mouse_marker[mouse_marker$cluster == y,"gene"])
  temp <- dplyr::union(temp_greenleaf,temp_PDhuman)
  temp <- dplyr::union(temp,temp_macaque)
  temp <- dplyr::union(temp,temp_mouse)
  return(length(temp))
}

for (i in rownames(overlap_heatmap)) {
  for (j in colnames(overlap_heatmap)) {
    overlap_heatmap[i,j] <- my_find_intersect(i,j)/my_find_union(i,j)
  }
}

row_cluster <- hclust(dist(overlap_heatmap),method = 'complete')
col_fun = colorRamp2(c(0,0.075,0.15), c('#3B4992FF','white','#EE0000FF'))

pdf(file = './res/step_13_fig_211111/all_species_all_marker_intersect_proportion_heatmap.pdf',width = 6,height = 5)
Heatmap(overlap_heatmap,show_column_names = TRUE,show_row_names = TRUE,cluster_rows = row_cluster,cluster_columns = row_cluster,
        col = col_fun,height = unit(3,'inches'),width = unit(3,'inches'),name = 'overlap proportion')
dev.off()

#cluster using express
gene_list <- dplyr::intersect(greenleaf_marker$gene,PDhuman_marker$gene)
gene_list <- dplyr::intersect(gene_list,macaque_marker$gene)
gene_list <- dplyr::intersect(gene_list,mouse_marker$gene)

greenleaf_bulk <- AverageExpression(object = greenleaf_RNA_seurat,assays = 'converted',features = gene_list,return.seurat = FALSE,group.by = 'basic_cell_type',slot = 'data')
greenleaf_bulk <- greenleaf_bulk$converted
colnames(greenleaf_bulk) <- paste('g',colnames(greenleaf_bulk),sep = '_')

PDhuman_bulk <- AverageExpression(object = PDhuman_RNA_seurat,assays = 'RNA',features = gene_list,return.seurat = FALSE,group.by = 'basic_cell_type',slot = 'data')
PDhuman_bulk <- PDhuman_bulk$RNA
colnames(PDhuman_bulk) <- paste('n',colnames(PDhuman_bulk),sep = '_')

macaque_bulk <- AverageExpression(object = macaque_RNA_seurat,assays = 'converted',features = gene_list,return.seurat = FALSE,group.by = 'basic_cell_type',slot = 'data')
macaque_bulk <- macaque_bulk$converted
colnames(macaque_bulk) <- paste('M',colnames(macaque_bulk),sep = '_')

mouse_bulk <- AverageExpression(object = mouse_RNA_seurat,assays = 'converted',features = gene_list,return.seurat = FALSE,group.by = 'basic_cell_type',slot = 'data')
mouse_bulk <- mouse_bulk$converted
colnames(mouse_bulk) <- paste('m',colnames(mouse_bulk),sep = '_')

temp <- cbind(greenleaf_bulk,PDhuman_bulk,macaque_bulk,mouse_bulk)
temp <- CreateSeuratObject(counts = temp,assay = 'RNA',min.cells = 0,min.features = 0,project = 'temp')
temp <- NormalizeData(object = temp,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)
temp <- temp@assays$RNA@data
temp <- t(scale(t(temp)))

temp_col <- ggsci::pal_npg()(4)
names(temp_col) <- c('greenleaf','Neuron','macaque','mouse')

col_anno <- HeatmapAnnotation(dataset = rep(c('greenleaf','Neuron','macaque','mouse'),each=10),
                              col = list(dataset = temp_col))

pdf(file = './res/step_13_fig_211111/all_species_all_cell_type_scale_together_heatmap.pdf',width = 10,height = 6)
Heatmap(as.matrix(temp),show_row_names = FALSE,show_column_names = TRUE,cluster_rows = TRUE,cluster_columns = TRUE,
        row_km = 10,clustering_method_columns = 'ward.D',name = 'z score',top_annotation = col_anno,
        height = unit(4,'inches'),width = unit(8,'inches'))
dev.off()

#scale separately
temp <- CreateSeuratObject(counts = greenleaf_bulk)
temp <- NormalizeData(temp)
greenleaf_bulk <- t(scale(t(temp@assays$RNA@data)))

temp <- CreateSeuratObject(counts = PDhuman_bulk)
temp <- NormalizeData(temp)
PDhuman_bulk <- t(scale(t(temp@assays$RNA@data)))

temp <- CreateSeuratObject(counts = macaque_bulk)
temp <- NormalizeData(temp)
macaque_bulk <- t(scale(t(temp@assays$RNA@data)))

temp <- CreateSeuratObject(counts = mouse_bulk)
temp <- NormalizeData(temp)
mouse_bulk <- t(scale(t(temp@assays$RNA@data)))

temp <- cbind(greenleaf_bulk,PDhuman_bulk,macaque_bulk,mouse_bulk)

pdf(file = './res/step_13_fig_211111/all_species_all_cell_type_scale_separately_heatmap.pdf',width = 10,height = 6)
Heatmap(as.matrix(temp),show_row_names = FALSE,show_column_names = TRUE,cluster_rows = TRUE,cluster_columns = TRUE,
        row_km = 10,clustering_method_columns = 'ward.D',name = 'z score',top_annotation = col_anno,
        height = unit(4,'inches'),width = unit(8,'inches'))
dev.off()