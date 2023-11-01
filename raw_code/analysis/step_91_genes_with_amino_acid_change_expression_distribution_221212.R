#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: genes with amino acid change expression distribution            ##
## Data: 2022.12.12                                                                ##
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
setwd('/content/data/sunym/project/Brain/')
.libPaths('/content/data/sunym/software/R_lib/R_4.1.3/')
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
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(UpSetR)
library(ggbreak)
library(ggvenn)
library(EnrichedHeatmap)
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(universalmotif)
library(topGO)
library(future.apply)
library(transPlotR)
library(aplot)
library(Biostrings)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# generate species combined dotplot ---------------------------------------

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
chimp_RNA_Seurat <- readRDS(file = './data/public/Organoid_single_cell_genomic_atlas_uncovers_human_specific_features_of_brain_development/processed_data/Chimp_filted_Seurat_object.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

human_to_macaque_anno <- readRDS(file = './res/step_91_fig_221205/macaque_homology_gene_list.rds')
human_to_mouse_anno <- readRDS(file = './res/step_91_fig_221205/mouse_homology_gene_list.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_221212.rds')

#gene list
Progenitor_gene <- c('KNL1','LMNB2','LAMA1','SPAG5','FOXM1','POLD3','CKAP5','KIF18A','CHEK1','CDK5RAP2','SCARA3',
                     'SCRN1','CALD1','PLEKHG1','MAST4','FAM114A1','KIF15','MCM2','GCA','PGM1','AKR7A2','COL11A1',
                     'VCAM1','PTBP1','TCF3','TOP2A','LAMP1','AHNAK','IFT74')
Ex_gene <- c('CHD3','SLC6A15','TBC1D30','TENM4','PCLO','SRRM3','EPHB6','NLN','GPRIN3','SLC8A1','PTPN4','PLXNA2','NFASC',
             'CSMD2','ERC1')
In_gene <- c('SLC32A1','C11orf96')
Glial_gene <- c('COL20A1','C3','NFKBID','MFSD12','PLD4','SH2B3','CTSC','SLCO2B1','DOCK8','ADAM28','TREM2','STAB1',
                'PTPRC','IL10RA')
Global_gene <- c('COMMD3','ATRX','IGBP1','FMR1','URI1','SAMD1','DHPS','TPGS2','SPAG9','MYEF2','SYNE2','NOVA1',
                 'CCND2','GOLM1','PHPT1','RB1CC1','KMT2C','ZNF292','DST','BRD2','VCAN','BOD1L1','MIF','RPL3',
                 'SPG7','SETD5','PRDM2')

#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_RNA_Seurat$ReAnno_celltype & cell_type_list %in% macaque_multiome_Seurat$cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

#filter cells
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$ReAnno_celltype %in% cell_type_list]
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% cell_type_list]
mouse_multiome_Seurat <- mouse_multiome_Seurat[,mouse_multiome_Seurat$macaque_cell_type %in% cell_type_list]

#add missing gene expression in human
Greenleaf_RNA_Seurat$ReAnno_celltype <- factor(Greenleaf_RNA_Seurat$ReAnno_celltype,levels = cell_type_list)

#add missing gene expression in chimp
gene_list <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene,Global_gene)
table(gene_list %in% rownames(chimp_RNA_Seurat@assays$RNA@counts))

chimp_RNA_Seurat$cell_type <- NA
chimp_RNA_Seurat@meta.data[chimp_RNA_Seurat$PredCellType == 'RG',"cell_type"] <- 'RG'
chimp_RNA_Seurat@meta.data[chimp_RNA_Seurat$PredCellType == 'IPC',"cell_type"] <- 'IP'
chimp_RNA_Seurat@meta.data[chimp_RNA_Seurat$PredCellType == 'EN',"cell_type"] <- 'Ex'
chimp_RNA_Seurat@meta.data[chimp_RNA_Seurat$PredCellType == 'IN',"cell_type"] <- 'In'
chimp_RNA_Seurat@meta.data[chimp_RNA_Seurat$PredCellType == 'Microglia',"cell_type"] <- 'Mic'
chimp_RNA_Seurat$cell_type <- factor(chimp_RNA_Seurat$cell_type,levels = c('RG','IP','Ex','In','Mic'))

#add missing gene expression in macaque
gene_list <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene,Global_gene)
gene_list <- gene_list[!(gene_list %in% names(human_to_macaque_anno))]
names(gene_list) <- gene_list
human_to_macaque_anno <- c(human_to_macaque_anno,gene_list)

gene_list <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene,Global_gene)
gene_list <- human_to_macaque_anno[gene_list]
gene_list <- gene_list[!(gene_list %in% rownames(macaque_multiome_Seurat@assays$RNA@counts))]
temp <- matrix(data = 0,nrow = length(gene_list),ncol = ncol(macaque_multiome_Seurat@assays$RNA@counts))
rownames(temp) <- gene_list
colnames(temp) <- colnames(macaque_multiome_Seurat@assays$RNA@counts)
temp <- as(temp,'dgCMatrix')

macaque_express_matrix <- rbind(macaque_multiome_Seurat@assays$RNA@counts,temp)
macaque_express_matrix <- CreateSeuratObject(counts = macaque_express_matrix,project = 'temp',assay = 'RNA',meta.data = macaque_multiome_Seurat@meta.data,min.cells = 0,min.features = 0)
macaque_express_matrix <- NormalizeData(object = macaque_express_matrix,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
macaque_express_matrix$cell_type <- factor(macaque_express_matrix$cell_type,levels = cell_type_list)

#add missing gene expression in mouse
gene_list <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene,Global_gene)
gene_list <- gene_list[!(gene_list %in% names(human_to_mouse_anno))]
i <- gene_list
gene_list <- str_to_title(gene_list)
names(gene_list) <- i
human_to_mouse_anno <- c(human_to_mouse_anno,gene_list)

gene_list <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene,Global_gene)
gene_list <- human_to_mouse_anno[gene_list]
gene_list <- gene_list[!(gene_list %in% rownames(mouse_multiome_Seurat@assays$RNA@counts))]
temp <- matrix(data = 0,nrow = length(gene_list),ncol = ncol(mouse_multiome_Seurat@assays$RNA@counts))
rownames(temp) <- gene_list
colnames(temp) <- colnames(mouse_multiome_Seurat@assays$RNA@counts)
temp <- as(temp,'dgCMatrix')

mouse_express_matrix <- rbind(mouse_multiome_Seurat@assays$RNA@counts,temp)
mouse_express_matrix <- CreateSeuratObject(counts = mouse_express_matrix,project = 'temp',assay = 'RNA',meta.data = mouse_multiome_Seurat@meta.data,min.cells = 0,min.features = 0)
mouse_express_matrix <- NormalizeData(object = mouse_express_matrix,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
mouse_express_matrix$macaque_cell_type <- factor(mouse_express_matrix$macaque_cell_type,levels = cell_type_list)

## progenitor gene ---------------------------------------------------------

#global setting
mk_gene <- Progenitor_gene

#get species marker gene
mk_gene_macaque <- human_to_macaque_anno[mk_gene]
names(mk_gene_macaque) <- NULL
mk_gene_mouse <- human_to_mouse_anno[mk_gene]
names(mk_gene_mouse) <- NULL

#get col.max and col.min
temp_human <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)
temp_chimp <- my_dotplot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))
col.min <- min(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#dot plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 3) + 
  theme(axis.title = element_blank())

temp_chimp <- DotPlot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['chimp']),
                      col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 3/5*12) + 
  theme(axis.title = element_blank())

temp_macaque <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,cols = c('lightgrey',color_param$species['macaque']),
                        col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 3) + 
  theme(axis.title = element_blank())

temp_mouse <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,cols = c('lightgrey',color_param$species['mouse']),
                      col.min = col.min,col.max = col.max,group.by = 'macaque_cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 3/10*12) + 
  theme(axis.title = element_blank())

pdf(file = './res/step_91_fig_221212/species_combined_dotplot_Progenitor_gene.pdf',width = 12,height = 10)
temp_human+temp_chimp+temp_macaque+temp_mouse+plot_layout(ncol = 4)
dev.off()

## Ex gene ---------------------------------------------------------

#global setting
mk_gene <- Ex_gene

#get species marker gene
mk_gene_macaque <- human_to_macaque_anno[mk_gene]
names(mk_gene_macaque) <- NULL
mk_gene_mouse <- human_to_mouse_anno[mk_gene]
names(mk_gene_mouse) <- NULL

#get col.max and col.min
temp_human <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)
temp_chimp <- my_dotplot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))
col.min <- min(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#dot plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 2) + 
  theme(axis.title = element_blank())

temp_chimp <- DotPlot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['chimp']),
                      col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 2/5*12) + 
  theme(axis.title = element_blank())

temp_macaque <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,cols = c('lightgrey',color_param$species['macaque']),
                        col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 2) + 
  theme(axis.title = element_blank())

temp_mouse <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,cols = c('lightgrey',color_param$species['mouse']),
                      col.min = col.min,col.max = col.max,group.by = 'macaque_cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 2/10*12) + 
  theme(axis.title = element_blank())

pdf(file = './res/step_91_fig_221212/species_combined_dotplot_Ex_gene.pdf',width = 12,height = 8)
temp_human+temp_chimp+temp_macaque+temp_mouse+plot_layout(ncol = 4)
dev.off()

## In gene ---------------------------------------------------------

#global setting
mk_gene <- In_gene

#get species marker gene
mk_gene_macaque <- human_to_macaque_anno[mk_gene]
names(mk_gene_macaque) <- NULL
mk_gene_mouse <- human_to_mouse_anno[mk_gene]
names(mk_gene_mouse) <- NULL

#get col.max and col.min
temp_human <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)
temp_chimp <- my_dotplot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))
col.min <- min(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#dot plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 0.25) + 
  theme(axis.title = element_blank())

temp_chimp <- DotPlot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['chimp']),
                      col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 0.25/5*12) + 
  theme(axis.title = element_blank())

temp_macaque <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,cols = c('lightgrey',color_param$species['macaque']),
                        col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 0.25) + 
  theme(axis.title = element_blank())

temp_mouse <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,cols = c('lightgrey',color_param$species['mouse']),
                      col.min = col.min,col.max = col.max,group.by = 'macaque_cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 0.25/10*12) + 
  theme(axis.title = element_blank())

pdf(file = './res/step_91_fig_221212/species_combined_dotplot_In_gene.pdf',width = 12,height = 8)
temp_human+temp_chimp+temp_macaque+temp_mouse+plot_layout(ncol = 4)
dev.off()

## Glial gene ---------------------------------------------------------

#global setting
mk_gene <- Glial_gene

#get species marker gene
mk_gene_macaque <- human_to_macaque_anno[mk_gene]
names(mk_gene_macaque) <- NULL
mk_gene_mouse <- human_to_mouse_anno[mk_gene]
names(mk_gene_mouse) <- NULL

#get col.max and col.min
temp_human <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)
temp_chimp <- my_dotplot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))
col.min <- min(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#dot plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 2) + 
  theme(axis.title = element_blank())

temp_chimp <- DotPlot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['chimp']),
                      col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 2/5*12) + 
  theme(axis.title = element_blank())

temp_macaque <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,cols = c('lightgrey',color_param$species['macaque']),
                        col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 2) + 
  theme(axis.title = element_blank())

temp_mouse <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,cols = c('lightgrey',color_param$species['mouse']),
                      col.min = col.min,col.max = col.max,group.by = 'macaque_cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 2/10*12) + 
  theme(axis.title = element_blank())

pdf(file = './res/step_91_fig_221212/species_combined_dotplot_Glial_gene.pdf',width = 12,height = 8)
temp_human+temp_chimp+temp_macaque+temp_mouse+plot_layout(ncol = 4)
dev.off()

## Global gene ---------------------------------------------------------

#global setting
mk_gene <- Global_gene

#get species marker gene
mk_gene_macaque <- human_to_macaque_anno[mk_gene]
names(mk_gene_macaque) <- NULL
mk_gene_mouse <- human_to_mouse_anno[mk_gene]
names(mk_gene_mouse) <- NULL

#get col.max and col.min
temp_human <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)
temp_chimp <- my_dotplot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))
col.min <- min(na.omit(temp_human$avg.exp.scaled,temp_chimp$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled))

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#dot plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 3) + 
  theme(axis.title = element_blank())

temp_chimp <- DotPlot(object = chimp_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['chimp']),
                      col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 3/5*12) + 
  theme(axis.title = element_blank())

temp_macaque <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,cols = c('lightgrey',color_param$species['macaque']),
                        col.min = col.min,col.max = col.max,group.by = 'cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 3) + 
  theme(axis.title = element_blank())

temp_mouse <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,cols = c('lightgrey',color_param$species['mouse']),
                      col.min = col.min,col.max = col.max,group.by = 'macaque_cell_type',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  theme(axis.ticks.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank()) + 
  theme(aspect.ratio = 3/10*12) + 
  theme(axis.title = element_blank())

pdf(file = './res/step_91_fig_221212/species_combined_dotplot_Global_gene.pdf',width = 12,height = 10)
temp_human+temp_chimp+temp_macaque+temp_mouse+plot_layout(ncol = 4)
dev.off()