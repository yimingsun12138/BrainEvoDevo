#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: investigate OPC                                                 ##
## Data: 2021.07.27                                                                ##
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
library(rtracklayer)
library(reshape2)
library(patchwork)
library(scales)
library(DESeq2)
library(dendextend)
library(ggsci)
library(circlize)
library(topGO)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

###############
##confirm vRG##
###############

# load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'

# subset
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('oRG','vRG','OPC','Cycling','IP')]
DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
# marker express
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id
marker_list <- c('PAX6','GFAP','VIM','NES','HOPX','TNC','FAM107A','EGFR','OLIG2','CSPG4','SOX10','PCDH15')
table(marker_list %in% rownames(macaque_RNA_seurat))

temp <- macaque_RNA_seurat[,macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_1'] > 0 & macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_2'] > -2.5]

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_OPC_dimplot.pdf',width = 4,height = 4)
DimPlot(temp,group.by = 'sub_cell_type',label = TRUE,repel = TRUE,cols = temp_col[c('oRG','vRG-1','vRG-2','OPC','Cyc-G2M','Cyc-S','IP')]) + 
  theme_classic() + 
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.position = 'none',
        aspect.ratio = 1) + 
  labs(title = '')
dev.off()

pdf(file = './res/step_5_fig_210727/OPC_marker_featureplot.pdf',width = 10,height = 7.5)
geneSetAveragePlot(genes = marker_list,object = temp,object.class = 'seurat',embedding = 'umap',reduction_key = 'UMAP',
                   assay = 'RNA',plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,
                   trim = NULL,print = TRUE,num.panel.rows = 3)
dev.off()

pdf(file = './res/step_5_fig_210727/OPC_marker_vlnplot.pdf',width = 12,height = 8)
VlnPlot(object = temp,features = marker_list,group.by = 'sub_cell_type')
dev.off()
# validate vRG
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
table(PDhuman_RNA_seurat$Cluster)
vRG_marker <- FindMarkers(object = PDhuman_RNA_seurat,ident.1 = 'vRG',group.by = 'Cluster',assay = 'RNA',slot = 'data',
                          logfc.threshold = 0.25,min.pct = 0.25,test.use = 'bimod',only.pos = TRUE)
vRG_marker <- vRG_marker[order(vRG_marker$pct.1,decreasing = TRUE),]
vRG_marker <- vRG_marker[order(vRG_marker$pct.2,decreasing = FALSE),]
DefaultAssay(macaque_RNA_seurat) <- 'RNA'

pdf(file = './res/step_5_fig_210727/PDhuman_vRG_marker_vlnplot.pdf',width = 12,height = 10)
VlnPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',features = rownames(vRG_marker)[1:20],assay = 'RNA',ncol = 4)
dev.off()

SZo_marker <- c('ACSBG1','BTBD17','PPP1R17','LYN','PLAC9','PLAGL1','PREX2','TNC')
SZo_marker %in% rownames(macaque_RNA_seurat)

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_OSVZ_marker_featureplot.pdf',width = 10,height = 6)
geneSetAveragePlot(genes = SZo_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,
                   print = TRUE)
dev.off()

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_OSVZ_marker_vlnplot.pdf',width = 10,height = 6)
VlnPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',features = SZo_marker,ncol = 4)
dev.off()

SZi_marker <- c('EOMES','PAX6','CDK6','GAS1','JAG1','PTH2R','SETD7')
SZi_marker %in% rownames(macaque_RNA_seurat)
geneSetAveragePlot(genes = SZi_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,
                   print = TRUE)
VlnPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',features = SZi_marker)


VZ_marker <- c('ST13','MGARP','CA3','CASP6','CD44','GLT8D2','PGAM2','POU4F1','SFRP2','GFAP')
VZ_marker %in% rownames(macaque_RNA_seurat)

pdf(file = './res/step_5_fig_210727/VZ_marker_express.pdf',width = 10,height = 7)
geneSetAveragePlot(genes = VZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,
                   print = TRUE)
dev.off()

pdf(file = './res/step_5_fig_210727/VZ_marker_express_vlnplot.pdf',width = 12,height = 7)
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
VlnPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',features = VZ_marker)
dev.off()

RG_marker <- c('SOX9','HES1','ATP1A2')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
table(RG_marker %in% rownames(macaque_RNA_seurat))

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_RG_marker_vlnplot.pdf',width = 8,height = 5)
VlnPlot(macaque_RNA_seurat,features = RG_marker,assay = 'RNA',group.by = 'sub_cell_type') + 
  theme(aspect.ratio = 1)
dev.off()
# label transfer
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')

table(macaque_RNA_seurat$cell_type)
table(PDhuman_RNA_seurat$Cluster)
PDhuman_RNA_seurat <- PDhuman_RNA_seurat[,PDhuman_RNA_seurat$Cluster %in% c('PgG2M','PgS','IP','OPC','oRG','vRG')]

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp
DefaultAssay(macaque_RNA_seurat) <- 'converted'

temp <- PDhuman_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno,filter_anno = FALSE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
PDhuman_RNA_seurat[['converted']] <- temp
DefaultAssay(PDhuman_RNA_seurat) <- 'converted'

macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',nfeatures = 2000,npcs = 50,preprocess = TRUE)
PDhuman_RNA_seurat <- my_process_seurat(object = PDhuman_RNA_seurat,assay = 'converted',nfeatures = 2000,npcs = 50,preprocess = TRUE)

OPC_anchor <- FindTransferAnchors(reference = PDhuman_RNA_seurat,query = macaque_RNA_seurat,dims = 1:20)
predictions <- TransferData(anchorset = OPC_anchor, refdata = PDhuman_RNA_seurat$Cluster, dims = 1:20)
macaque_RNA_seurat <- AddMetaData(object = macaque_RNA_seurat,metadata = predictions)

pdf(file = './res/step_5_fig_210727/PDhuman_annotation_label_transfer.pdf',width = 5,height = 4)
scibet::Confusion_heatmap(ori = macaque_RNA_seurat$sub_cell_type,prd = macaque_RNA_seurat$predicted.id) + 
  theme(aspect.ratio = (6/7),
        axis.title = element_text()) + 
  xlab('Original annotation') + ylab('Predicted annotation')
dev.off()

# macaque_RNA_seurat vRG marker
temp <- macaque_RNA_seurat[,macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_1'] > 0 & macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_2'] > -2.5]
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
vRG_1_marker <- FindMarkers(object = macaque_RNA_seurat,ident.1 = 'vRG-1',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',
                            logfc.threshold = 0.25,min.pct = 0.25,test.use = 'bimod',only.pos = TRUE)

vRG_1_marker <- vRG_1_marker[order(vRG_1_marker$pct.1,decreasing = TRUE),]
vRG_1_marker <- vRG_1_marker[order(vRG_1_marker$pct.2,decreasing = FALSE),]

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_vRG_1_marker_vlnplot.pdf',width = 15,height = 8)
VlnPlot(macaque_RNA_seurat,features = rownames(vRG_1_marker)[1:12],assay = 'RNA',group.by = 'sub_cell_type')
dev.off()

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_vRG_1_marker_featureplot.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = rownames(vRG_1_marker)[1:12],object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,print = TRUE,num.panel.rows = 3,plot.type = 'panels')
dev.off()


vRG_2_marker <- FindMarkers(object = macaque_RNA_seurat,ident.1 = 'vRG-2',group.by = 'sub_cell_type',assay = 'RNA',slot = 'data',
                            logfc.threshold = 0.25,min.pct = 0.25,test.use = 'bimod',only.pos = TRUE)
vRG_2_marker <- vRG_2_marker[order(vRG_2_marker$pct.1,decreasing = TRUE),]
vRG_2_marker <- vRG_2_marker[order(vRG_2_marker$pct.2,decreasing = FALSE),]

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_vRG_2_marker_vlnplot.pdf',width = 15,height = 8)
VlnPlot(macaque_RNA_seurat,features = rownames(vRG_2_marker)[1:12],assay = 'RNA',group.by = 'sub_cell_type')
dev.off()

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_vRG_2_marker_featureplot.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = rownames(vRG_2_marker)[1:12],object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,print = TRUE,num.panel.rows = 3,plot.type = 'panels')
dev.off()

#try label tranfer from greenleaf
# greenleaf_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat.rds')
# DimPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'UMAP')
# 
# #re annotate
# temp <- greenleaf_RNA_seurat@assays$originalexp@counts
# temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno,filter_anno = FALSE)
# temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
# greenleaf_RNA_seurat[['converted']] <- temp
# DefaultAssay(greenleaf_RNA_seurat) <- 'converted'
# greenleaf_RNA_seurat <- my_process_seurat(object = greenleaf_RNA_seurat,assay = 'converted',nfeatures = 2000,npcs = 50,preprocess = TRUE)
# greenleaf_RNA_seurat$cell_type <- greenleaf_RNA_seurat$seurat_clusters
# greenleaf_RNA_seurat@meta.data[,"cell_type"] <- as.character(greenleaf_RNA_seurat@meta.data[,"cell_type"])
# #Cyc
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('TOP2A','MKI67','CLSPN','AURKA'))
# #c8
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c8'),"cell_type"] <- 'Cyc'
# 
# #early RG
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('SOX9','HES1','ATP1A2','FBXO32','CCN2','CCN1','MOXD1','HOPX','FAM107A','MT3'))
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('NPY','FGFR3'))
# #c6
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c6'),"cell_type"] <- 'Early RG'
# 
# #late RG
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('CD9','GPX3','TNC'))
# #c10
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c10'),"cell_type"] <- 'Late RG'
# 
# #tRG
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('CRYAB','NR4A1','FOXJ1'))
# #c18
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c18'),"cell_type"] <- 'tRG'
# 
# #OPC
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('SOX10','NKX2-2','MBP'))
# #c17
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c17'),"cell_type"] <- 'OPC'
# 
# #oIPC
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('ASCL1','OLIG2','PDGFRA','EGFR'))
# #c11
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c11'),"cell_type"] <- 'oIPC'
# 
# #nIPC
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('EOMES','PPP1R17','PENK','NEUROG1','NEUROG2'))
# #c14
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c14'),"cell_type"] <- 'nIPC'
# 
# #SP
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('NR4A2','CRYM','ST18','CDH18'))
# #c13
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c13'),"cell_type"] <- 'SP-GluN'
# 
# #InMGE
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('LHX6','SST'))
# #c3 c21
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c3','c21'),"cell_type"] <- 'InMGE'
# 
# #InCGE
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('SP8','NR2F2'))
# #c1
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c1'),"cell_type"] <- 'InCGE'
# 
# #GluN
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('NEUROD2','TBR1','BCL11B','SATB2','SLC17A7'))
# #c2,c0,c9,c5,c7,c4,c12,c15
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c2'),"cell_type"] <- 'GluN-1'
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c0'),"cell_type"] <- 'GluN-2'
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c9'),"cell_type"] <- 'GluN-3'
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c5'),"cell_type"] <- 'GluN-4'
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c7'),"cell_type"] <- 'GluN-5'
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c4'),"cell_type"] <- 'GluN-6'
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c12'),"cell_type"] <- 'GluN-7'
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c15'),"cell_type"] <- 'GluN-8'
# 
# #Mic
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('AIF1','CCL3'))
# #c16
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c16'),"cell_type"] <- 'Mic'
# 
# #End
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('CLDN5','PECAM1'))
# #c20
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c20'),"cell_type"] <- 'End'
# 
# #Per
# VlnPlot(greenleaf_RNA_seurat,group.by = 'seurat_clusters',features = c('FOXC2','PDGFRB'))
# #c19,c22
# greenleaf_RNA_seurat@meta.data[greenleaf_RNA_seurat$seurat_clusters %in% c('c19','c22'),"cell_type"] <- 'Per'
# 
# #dimplot
# DimPlot(greenleaf_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP')
# table(greenleaf_RNA_seurat$cell_type)
# saveRDS(greenleaf_RNA_seurat,file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')

#try label transfer from greenleaf data
greenleaf_human_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
DefaultAssay(greenleaf_human_seurat) <- 'converted'
greenleaf_human_seurat <- greenleaf_human_seurat[,greenleaf_human_seurat$cell_type %in% c('Cyc','Early RG','Late RG','tRG','oIPC','nIPC','OPC')]

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp
DefaultAssay(macaque_RNA_seurat) <- 'converted'

macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',nfeatures = 2000,npcs = 50,preprocess = TRUE)
greenleaf_human_seurat <- my_process_seurat(object = greenleaf_human_seurat,assay = 'converted',nfeatures = 2000,npcs = 50,preprocess = TRUE)

OPC_anchor <- FindTransferAnchors(reference = greenleaf_human_seurat,query = macaque_RNA_seurat,dims = 1:30)
predictions <- TransferData(anchorset = OPC_anchor, refdata = greenleaf_human_seurat$cell_type, dims = 1:30)
macaque_RNA_seurat <- AddMetaData(object = macaque_RNA_seurat,metadata = predictions)

pdf(file = './res/step_5_fig_210727/greenleaf_human_seurat_label_transfer.pdf',width = 5,height = 5)
scibet::Confusion_heatmap(ori = macaque_RNA_seurat$sub_cell_type,prd = macaque_RNA_seurat$predicted.id) + 
  theme(aspect.ratio = 1,
        axis.title = element_text()) + 
  xlab('Original annotation') + ylab('Predicted annotation')
dev.off()

#vRG in greenleaf data
greenleaf_human_seurat <- greenleaf_human_seurat[,greenleaf_human_seurat$Age == 'pcw16']
greenleaf_human_seurat <- greenleaf_human_seurat[,greenleaf_human_seurat$cell_type == 'Early RG']
greenleaf_human_seurat <- my_process_seurat(object = greenleaf_human_seurat,assay = 'originalexp',nfeatures = 2000,npcs = 50,preprocess = TRUE)
greenleaf_human_seurat <- my_process_seurat(object = greenleaf_human_seurat,assay = 'originalexp',preprocess = FALSE,group.by = 'seurat_clusters',dim_to_use = 30,resolution = 0.4,label = TRUE)
DimPlot(greenleaf_human_seurat,group.by = 'cell_type',reduction = 'umap',label = TRUE,repel = TRUE)
DimPlot(greenleaf_human_seurat,group.by = 'seurat_clusters',reduction = 'umap',label = TRUE,repel = TRUE)
DimPlot(greenleaf_human_seurat,group.by = 'Sample.ID',reduction = 'umap',label = TRUE,repel = TRUE)
DefaultAssay(greenleaf_human_seurat) <- 'converted'

geneSetAveragePlot(genes = c('SOX9','HES1','ATP1A2','MOXD1','HOPX','FAM107A','MT3'),object = greenleaf_human_seurat,object.class = 'seurat',
                   assay = 'converted',embedding = 'umap',reduction_key = 'umap',scaled = FALSE,plot.type = 'panels',
                   color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,print = TRUE)

geneSetAveragePlot(genes = c('SOX9','HES1','ATP1A2','FBXO32','CCN1','CCN2'),object = greenleaf_human_seurat,object.class = 'seurat',
                   assay = 'converted',embedding = 'umap',reduction_key = 'umap',scaled = FALSE,plot.type = 'panels',
                   color.palette = c('lightgrey','blue'),aspectratio = 1,trim = NULL,print = TRUE)

greenleaf_human_seurat$sub_cell_type <- greenleaf_human_seurat$cell_type
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('0'),"sub_cell_type"] <- 'vRG-1'
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('1'),"sub_cell_type"] <- 'vRG-2'
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('3'),"sub_cell_type"] <- 'oRG-1'
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('4'),"sub_cell_type"] <- 'oRG-2'
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('5'),"sub_cell_type"] <- 'oRG-3'
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('6'),"sub_cell_type"] <- 'oRG-4'
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('2'),"sub_cell_type"] <- 'oRG-5'

pdf(file = './res/step_5_fig_210727/greenleaf_human_seurat_pcw16_Early_RG_sub_cell_type_dimplot.pdf',width = 6,height = 5)
DimPlot(greenleaf_human_seurat,group.by = 'sub_cell_type',reduction = 'umap',label = TRUE,repel = TRUE)
dev.off()

pdf(file = './res/step_5_fig_210727/greenleaf_human_seurat_pcw16_Early_RG_sub_cell_type_vlnplot.pdf',width = 12,height = 8)
VlnPlot(greenleaf_human_seurat,features = c('SOX9','HES1','ATP1A2','MOXD1','HOPX','FAM107A','MT3','FBXO32','CCN1','CCN2'),group.by = 'sub_cell_type')
dev.off()

greenleaf_human_seurat$cell_type <- 'oRG'
greenleaf_human_seurat@meta.data[greenleaf_human_seurat$seurat_clusters %in% c('0','1'),"cell_type"] <- 'vRG'

#label transfer
greenleaf_pcw16_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
greenleaf_pcw16_seurat <- greenleaf_pcw16_seurat[,greenleaf_pcw16_seurat$Age == 'pcw16']
greenleaf_pcw16_seurat <- greenleaf_pcw16_seurat[,greenleaf_pcw16_seurat$cell_type %in% c('Cyc','Early RG','oIPC','nIPC','OPC')]

DimPlot(greenleaf_pcw16_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP')

greenleaf_pcw16_seurat@meta.data[colnames(greenleaf_human_seurat[,greenleaf_human_seurat$cell_type == 'oRG']),"cell_type"] <- 'oRG'
greenleaf_pcw16_seurat@meta.data[colnames(greenleaf_human_seurat[,greenleaf_human_seurat$cell_type == 'vRG']),"cell_type"] <- 'vRG'
greenleaf_pcw16_seurat@meta.data[greenleaf_pcw16_seurat$cell_type %in% c('OPC','oIPC'),"cell_type"] <- 'oIPC/OPC'

macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type != 'IP']
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',nfeatures = 2000,npcs = 50,preprocess = TRUE)
greenleaf_pcw16_seurat <- my_process_seurat(object = greenleaf_pcw16_seurat,assay = 'converted',nfeatures = 2000,npcs = 50,preprocess = TRUE)

OPC_anchor <- FindTransferAnchors(reference = greenleaf_pcw16_seurat,query = macaque_RNA_seurat,dims = 1:30)
predictions <- TransferData(anchorset = OPC_anchor, refdata = greenleaf_pcw16_seurat$cell_type, dims = 1:30)
macaque_RNA_seurat <- AddMetaData(object = macaque_RNA_seurat,metadata = predictions)

pdf(file = './res/step_5_fig_210727/greenleaf_human_seurat_pcw16_sub_cell_type_confusion_heatmap.pdf',width = 5,height = 4)
scibet::Confusion_heatmap(ori = macaque_RNA_seurat$sub_cell_type,prd = macaque_RNA_seurat$predicted.id) + 
  theme(aspect.ratio = (4/6),
        axis.title = element_text()) + 
  xlab('Original annotation') + ylab('Predicted annotation')
dev.off()


# try macaque bulk RNA-seq data
macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

layer_marker <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')
show_col(ggsci::pal_futurama()(10))

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

macaque_RNA_seurat <- My_add_module_score(seu.obj = macaque_RNA_seurat,assay = 'RNA',features = OSVZ_marker,meta_var = 'OSVZ',scale = FALSE,center = TRUE)

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_OSVZ_module_score_boxplot.pdf',width = 4,height = 4)
ggplot(data = macaque_RNA_seurat@meta.data,aes(x=sub_cell_type,y=OSVZ,fill=sub_cell_type)) + 
  geom_boxplot(position = 'identity',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('vRG-1','vRG-2','oRG','OPC','Cyc-S','Cyc-G2M','IP')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  CenterTitle() + xlab('cell type') + ylab('OSVZ signature') + labs(title = 'centered log1p(mean counts)') + 
  geom_hline(yintercept = 0,linetype = 'dashed',size = 1,color = 'brown')
dev.off()
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

macaque_RNA_seurat <- My_add_module_score(seu.obj = macaque_RNA_seurat,assay = 'RNA',features = ISVZ_marker,meta_var = 'ISVZ',scale = FALSE,center = TRUE)

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_ISVZ_module_score_boxplot.pdf',width = 4,height = 4)
ggplot(data = macaque_RNA_seurat@meta.data,aes(x=sub_cell_type,y=ISVZ,fill=sub_cell_type)) + 
  geom_boxplot(position = 'identity',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('vRG-1','vRG-2','oRG','OPC','Cyc-S','Cyc-G2M','IP')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  CenterTitle() + xlab('cell type') + ylab('ISVZ signature') + labs(title = 'centered log1p(mean counts)') + 
  geom_hline(yintercept = 0,linetype = 'dashed',size = 1,color = 'brown')
dev.off()
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

macaque_RNA_seurat <- My_add_module_score(seu.obj = macaque_RNA_seurat,assay = 'RNA',features = VZ_marker,meta_var = 'VZ',scale = FALSE,center = TRUE)

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_VZ_module_score_boxplot.pdf',width = 4,height = 4)
ggplot(data = macaque_RNA_seurat@meta.data,aes(x=sub_cell_type,y=VZ,fill=sub_cell_type)) + 
  geom_boxplot(position = 'identity',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('vRG-1','vRG-2','oRG','OPC','Cyc-S','Cyc-G2M','IP')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  CenterTitle() + xlab('cell type') + ylab('VZ signature') + labs(title = 'centered log1p(mean counts)') + 
  geom_hline(yintercept = 0,linetype = 'dashed',size = 1,color = 'brown')
dev.off()

# using more strict threshold to find marker

#load data
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')
colnames(macaque_express_matrix) <- sub(pattern = '-',replacement = '_',colnames(macaque_express_matrix))

macaque_coldata <- data.frame(species='macaque',sample=colnames(macaque_express_matrix))
rownames(macaque_coldata) <- colnames(macaque_express_matrix)
donor_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][1]
  return(temp)
})
donor_list <- unlist(donor_list)
macaque_coldata$donor <- donor_list
layer_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][2]
  return(temp)
})
layer_list <- unlist(layer_list)
macaque_coldata$layer <- layer_list
#coldata for DE analysis
layer_list <- unique(layer_list)
temp <- do.call(cbind,base::lapply(layer_list,function(x){
  temp_list <- as.character(macaque_coldata$layer == x)
  temp_list[which(temp_list == 'TRUE')] <- 'YES'
  temp_list[which(temp_list == 'FALSE')] <- 'NO'
  return(temp_list)
}))
colnames(temp) <- layer_list
rownames(temp) <- rownames(macaque_coldata)
macaque_coldata <- cbind(macaque_coldata,temp)
#find DE
layer_list <- c('OSVZ','ISVZ','VZ')
DE_list <- base::lapply(layer_list,function(x){
  formula_char <- paste0('~','donor','+',as.character(x))
  formula_char <- as.formula(formula_char)
  temp <- DESeqDataSetFromMatrix(countData = macaque_express_matrix,colData = macaque_coldata,design = formula_char)
  temp <- DESeq(temp)
  temp <- results(temp,contrast = c(as.character(x),'YES','NO'),lfcThreshold = 0.5,alpha = 0.05)
  temp <- temp[!is.na(temp$padj),]
  temp <- temp[!is.na(temp$log2FoldChange),]
  temp <- temp[(temp$log2FoldChange > 0.5) & (temp$padj < 0.05),]
  return(temp)
})
names(DE_list) <- layer_list

temp <- macaque_RNA_seurat[,macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_1'] > 0 & macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_2'] > -2.5]

#OSVZ
OSVZ_marker <- DE_list$OSVZ
OSVZ_marker <- OSVZ_marker[order(OSVZ_marker$padj,decreasing = FALSE),]
OSVZ_marker <- rownames(OSVZ_marker)[1:12]
OSVZ_marker <- unique(OSVZ_marker)
table(OSVZ_marker %in% macaque_annotation$gene_id)
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% macaque_annotation$gene_id]
OSVZ_marker <- macaque_annotation$gene_name[base::match(OSVZ_marker,table = macaque_annotation$gene_id)]
table(OSVZ_marker %in% rownames(macaque_RNA_seurat))
OSVZ_marker <- OSVZ_marker[OSVZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(OSVZ_marker))
OSVZ_marker <- unique(OSVZ_marker)

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_bulk_OSVZ_marker_featureplot.pdf',width = 9,height = 6)
geneSetAveragePlot(genes = OSVZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE,num.panel.rows = 3)
dev.off()

#ISVZ
ISVZ_marker <- DE_list$ISVZ
ISVZ_marker <- ISVZ_marker[order(ISVZ_marker$padj,decreasing = FALSE),]
ISVZ_marker <- rownames(ISVZ_marker)[1:12]
ISVZ_marker <- unique(ISVZ_marker)
table(ISVZ_marker %in% macaque_annotation$gene_id)
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% macaque_annotation$gene_id]
ISVZ_marker <- macaque_annotation$gene_name[base::match(ISVZ_marker,table = macaque_annotation$gene_id)]
table(ISVZ_marker %in% rownames(macaque_RNA_seurat))
ISVZ_marker <- ISVZ_marker[ISVZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(ISVZ_marker))
ISVZ_marker <- unique(ISVZ_marker)

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_bulk_ISVZ_marker_featureplot.pdf',width = 9,height = 6)
geneSetAveragePlot(genes = ISVZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE,num.panel.rows = 3)
dev.off()

#VZ
VZ_marker <- DE_list$VZ
VZ_marker <- VZ_marker[order(VZ_marker$padj,decreasing = FALSE),]
VZ_marker <- rownames(VZ_marker)[1:12]
VZ_marker <- unique(VZ_marker)
table(VZ_marker %in% macaque_annotation$gene_id)
VZ_marker <- VZ_marker[VZ_marker %in% macaque_annotation$gene_id]
VZ_marker <- macaque_annotation$gene_name[base::match(VZ_marker,table = macaque_annotation$gene_id)]
table(VZ_marker %in% rownames(macaque_RNA_seurat))
VZ_marker <- VZ_marker[VZ_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(VZ_marker))
VZ_marker <- unique(VZ_marker)

pdf(file = './res/step_5_fig_210727/macaque_RNA_seurat_bulk_VZ_marker_featureplot.pdf',width = 9,height = 6)
geneSetAveragePlot(genes = VZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE,num.panel.rows = 3)
dev.off()

# make heatmap to validate the marker
macaque_annotation <- macaque_annotation[,c('gene_name','gene_id')]
macaque_annotation <- cbind(macaque_annotation,macaque_annotation)

macaque_express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = macaque_express_matrix,anno = macaque_annotation,filter_anno = FALSE)
macaque_express_matrix <- macaque_express_matrix[c(OSVZ_marker,ISVZ_marker,VZ_marker),c('10068A_ISVZ','10068B_ISVZ','11002A_ISVZ','11002B_ISVZ',
                                                                                        '10068A_OSVZ','10068B_OSVZ','11002A_OSVZ','11002B_OSVZ',
                                                                                        '10068A_VZ','10068B_VZ','11002A_VZ','11002B_VZ')]
macaque_express_matrix <- log1p(macaque_express_matrix)
macaque_express_matrix <- t(base::scale(t(macaque_express_matrix),center = TRUE,scale = TRUE))

#cluster sample
col_dend <- as.dendrogram(hclust(dist(t(macaque_express_matrix)),method = 'complete'))
col_dend <- reorder(col_dend,c(1,1,1,1,3,3,3,3,2,2,2,2))
#cluster gene
row_dend <- as.dendrogram(hclust(dist(macaque_express_matrix),method = 'complete'))
#create annotation
temp_col <- ggsci::pal_aaas()(3)
names(temp_col) <- c('VZ','ISVZ','OSVZ')
sample_anno <- HeatmapAnnotation(sample = sub(pattern = '^.*_',replacement = '',colnames(macaque_express_matrix)),
                                 col = list(sample = temp_col))
gene_anno <- rowAnnotation(marker = c(rep('OSVZ',12),rep('ISVZ',12),rep('VZ',12)),
                           col = list(marker = temp_col))
#plot
col_fun <- colorRamp2(c(-2,0,2),c('#3B4992FF','white','#EE0000FF'))

pdf(file = './res/step_5_fig_210727/macaque_bulk_RNA_marker_heatmap.pdf',width = 8,height = 9)
Heatmap(macaque_express_matrix,cluster_rows = row_dend,cluster_columns = col_dend,show_row_names = TRUE,show_column_names = TRUE,
        top_annotation = sample_anno,left_annotation = gene_anno,height = unit(6,units = 'inches'),width = unit(3,units = 'inches'),
        col = col_fun,rect_gp = gpar(col = 'black'),name = 'z score')
dev.off()
#conclusion: no vRG in macaque_RNA_seurat

###################################
##investigate G1 proportion in RG##
###################################
# #preprocess mouse data
# mouse_RNA_seurat <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/mouse_RNA_seurat.rds')
# mouse_RNA_seurat <- NormalizeData(object = mouse_RNA_seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)
# mouse_RNA_seurat <- FindVariableFeatures(object = mouse_RNA_seurat,selection.method = 'vst',assay = 'RNA',nfeatures = 3000)
# mouse_RNA_seurat <- ScaleData(object = mouse_RNA_seurat,features = VariableFeatures(mouse_RNA_seurat),
#                               vars.to.regress = c('nCount_RNA','nFeature_RNA','percent_mito','CC_Difference'),do.scale = TRUE,do.center = TRUE)
# mouse_RNA_seurat <- RunPCA(object = mouse_RNA_seurat,assay = 'RNA',features = VariableFeatures(mouse_RNA_seurat),npcs = 50)
# 
# ElbowPlot(mouse_RNA_seurat,ndims = 50)
# 
# mouse_RNA_seurat <- FindNeighbors(object = mouse_RNA_seurat,dims = 1:50)
# mouse_RNA_seurat <- RunUMAP(object = mouse_RNA_seurat,dims = 1:39)
# DimPlot(mouse_RNA_seurat,group.by = 'New_cellType',label = FALSE)
# 
# saveRDS(mouse_RNA_seurat,file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/mouse_RNA_seurat.rds')


#load mouse data
mouse_RNA_seurat <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/mouse_RNA_seurat.rds')
mouse_RNA_seurat$Age <- as.character(sub(pattern = '_S.',replacement = '',mouse_RNA_seurat$biosample_id))
DefaultAssay(mouse_RNA_seurat) <- 'RNA'

pdf(file = './res/step_5_fig_210727/mouse_RNA_seurat_cell_type_dimplot.pdf',width = 12,height = 10)
DimPlot(mouse_RNA_seurat,group.by = 'New_cellType',label = TRUE,repel = TRUE)
dev.off()

DimPlot(mouse_RNA_seurat,group.by = 'Phase',label = FALSE)

#load human data
human_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
DefaultAssay(human_RNA_seurat) <- 'converted'

pdf(file = './res/step_5_fig_210727/greenleaf_RNA_seurat_cell_type_dimplot.pdf',width = 12,height = 10)
DimPlot(human_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP')
dev.off()

g2m.signature <- cc.genes$g2m.genes
g2m.signature <- g2m.signature[g2m.signature %in% rownames(human_RNA_seurat)]
s.signature <- cc.genes$s.genes
s.signature <- s.signature[s.signature %in% rownames(human_RNA_seurat)]

human_RNA_seurat <- CellCycleScoring(object = human_RNA_seurat,s.features = s.signature,g2m.features = g2m.signature)
DimPlot(human_RNA_seurat,group.by = 'Phase',label = FALSE,reduction = 'UMAP')

#load macaque_RNA_seurat
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'

# RG cell number in different time
temp <- mouse_RNA_seurat[,mouse_RNA_seurat$New_cellType == 'Apical progenitors']
table(temp$donor_id)

temp <- human_RNA_seurat[,human_RNA_seurat$cell_type %in% c('Early RG','Late RG','tRG')]
table(temp$Age)

temp <- My_Cell_Proportion(seu.obj = mouse_RNA_seurat,group.by = 'New_cellType',split.by = 'Age')
ggplot(data = temp, mapping = aes(x = Age,y = Proportion,fill = New_cellType))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.95) + 
  theme_classic() + 
  theme(aspect.ratio = 2.5,
        axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) + 
  labs(fill = 'cell type') + xlab('') + ylab('proportion') + 
  coord_flip()

temp <- My_Cell_Proportion(seu.obj = human_RNA_seurat,group.by = 'cell_type',split.by = 'Age')
ggplot(data = temp, mapping = aes(x = Age,y = Proportion,fill = cell_type))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.95) + 
  theme_classic() + 
  theme(aspect.ratio = 2.5,
        axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) + 
  labs(fill = 'cell type') + xlab('') + ylab('proportion') + 
  coord_flip()

temp <- temp[temp$cell_type %in% c('Early RG','Late RG','tRG'),]

for (i in c("pcw16","pcw20","pcw21","pcw24")) {
  j <- temp[temp$Age == i,]
  j <- sum(j$Proportion)
  print(j)
}

#cell cycle Proportion in different time
#human
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

temp <- human_RNA_seurat[,human_RNA_seurat$cell_type %in% c('Early RG','Late RG','tRG')]
temp <- My_Cell_Proportion(seu.obj = temp,group.by = 'Phase',split.by = 'Age')

pdf(file = './res/step_5_fig_210727/greenleaf_RG_seurat_cell_cycle_Proportion.pdf',width = 6,height = 2.5)
ggplot(data = temp, mapping = aes(x = Age,y = Proportion,fill = Phase))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.75) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6,
        axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) + 
  labs(fill = 'cell type') + xlab('') + ylab('proportion') + 
  coord_flip()
dev.off()

#mouse
temp <- mouse_RNA_seurat[,mouse_RNA_seurat$New_cellType == 'Apical progenitors']
temp <- My_Cell_Proportion(seu.obj = temp,group.by = 'Phase',split.by = 'Age')

pdf(file = './res/step_5_fig_210727/mouse_RG_seurat_cell_cycle_proportion.pdf',width = 6,height = 3)
ggplot(data = temp, mapping = aes(x = Age,y = Proportion,fill = Phase))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.75) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 0.6,
        axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12)) + 
  labs(fill = 'cell type') + xlab('') + ylab('proportion') + 
  coord_flip()
dev.off()

# combine
temp <- human_RNA_seurat[,human_RNA_seurat$cell_type %in% c('Early RG','Late RG','tRG')]
temp <- My_Cell_Proportion(seu.obj = temp,group.by = 'Phase',split.by = 'Age')

p1 <- ggplot(data = temp, mapping = aes(x = Age,y = Proportion,fill = Phase))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.75) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        aspect.ratio = 0.6) + 
  labs(fill = 'cell type') + xlab('') + ylab('Human RG cell cycle') + 
  coord_flip()

temp <- mouse_RNA_seurat[,mouse_RNA_seurat$New_cellType == 'Apical progenitors']
temp <- My_Cell_Proportion(seu.obj = temp,group.by = 'Phase',split.by = 'Age')

p2 <- ggplot(data = temp, mapping = aes(x = Age,y = Proportion,fill = Phase))+
  geom_bar(stat = 'identity',position = 'stack',width = 0.75) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(family="Times",face = 'bold',size = 14),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        aspect.ratio = 1.05) + 
  labs(fill = 'cell type') + xlab('') + ylab('Mouse RG cell cycle') + 
  coord_flip()

pdf(file = './res/step_5_fig_210727/combined_greenleaf_mouse_RG_cell_cycle_proportion.pdf',width = 4.5,height = 5)
p1+p2+plot_layout(nrow = 2)
dev.off()

# G1 Phase RG in VZ or OSVZ?
SVZ_marker <- c('ACSBG1','BTBD17','PPP1R17','LYN','PLAC9','PLAGL1','PREX2','TNC','EOMES','PAX6','CDK6','GAS1','JAG1','PTH2R','SETD7')
VZ_marker <- c('ST13','MGARP','CA3','CASP6','CD44','GLT8D2','PGAM2','POU4F1','SFRP2','GFAP')

temp <- human_RNA_seurat[,human_RNA_seurat$cell_type %in% c('Early RG','Late RG','tRG')]

temp <- My_add_module_score(seu.obj = temp,assay = 'converted',features = SVZ_marker,meta_var = 'SVZ_marker',scale = FALSE,center = FALSE)
temp <- My_add_module_score(seu.obj = temp,assay = 'converted',features = VZ_marker,meta_var = 'VZ_marker',scale = FALSE,center = FALSE)
temp <- temp@meta.data

p1 <- ggplot(temp,aes(x=Age,y=SVZ_marker,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'SVZ layer marker') + xlab('Age') + ylab('Module score')

p2 <- ggplot(temp,aes(x=Age,y=VZ_marker,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'VZ layer marker') + xlab('Age') + ylab('Module score')

pdf(file = './res/step_5_fig_210727/greenleaf_RG_seurat_layer_marker_score.pdf',width = 8,height = 4)
p1+p2+plot_layout(ncol = 2)
dev.off()

#G1 phase RG proliferation or differentiation?
#do the GO
#G1 enrich
temp <- human_RNA_seurat[,human_RNA_seurat$cell_type %in% c('Early RG','Late RG','tRG')]
G1_enrich <- FindMarkers(object = temp,ident.1 = 'G1',ident.2 = c('S','G2M'),group.by = 'Phase',assay = 'converted',
                         logfc.threshold = 0.25,test.use = 'bimod',min.pct = 0.4,only.pos = TRUE)

gene_list <- rownames(human_RNA_seurat@assays$converted@counts)
gene_list <- c(gene_list %in% rownames(G1_enrich))
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(human_RNA_seurat@assays$converted@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

GO_G1_enrich <- new("topGOdata",
                    description = "UP regulated",
                    ontology = "BP",
                    allGenes = gene_list,
                    nodeSize = 10,annotationFun = annFUN.org,
                    mapping = 'org.Hs.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_G1_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_G1_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = allRes$Term)

p1 <- ggplot(allRes[c(1,5,6,8,14,18,29,41,42),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#D62728FF',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.8,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('G1 enrich')
#S/G2M enrich
division_enrich <- FindMarkers(object = temp,ident.2 = 'G1',ident.1 = c('S','G2M'),group.by = 'Phase',assay = 'converted',
                               logfc.threshold = 0.25,test.use = 'bimod',min.pct = 0.4,only.pos = TRUE)

gene_list <- rownames(human_RNA_seurat@assays$converted@counts)
gene_list <- c(gene_list %in% rownames(division_enrich))
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(human_RNA_seurat@assays$converted@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

GO_division_enrich <- new("topGOdata",
                          description = "UP regulated",
                          ontology = "BP",
                          allGenes = gene_list,
                          nodeSize = 10,annotationFun = annFUN.org,
                          mapping = 'org.Hs.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_division_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_division_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = allRes$Term)

p2 <- ggplot(allRes[c(3,6:9,14:17),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#2CA02CFF',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.8,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('S/G2M enrich')


pdf(file = './res/step_5_fig_210727/greenleaf_RG_seurat_cell_cycle_enrich_GO.pdf',width = 7,height = 7)
p1+p2+plot_layout(nrow = 2)
dev.off()

#get GO list from mouse
developmental_process_gene_list <- getGOgeneSet(x = 'developmental process',OrgDb = 'org.Mm.eg.db',ont = 'BP',keytype = 'SYMBOL')
anatomical_structure_development_gene_list <- getGOgeneSet(x = 'anatomical structure development',OrgDb = 'org.Mm.eg.db',ont = 'BP',keytype = 'SYMBOL')
localization_gene_list <- getGOgeneSet(x = 'localization',OrgDb = 'org.Mm.eg.db',ont = 'BP',keytype = 'SYMBOL')
neurogenesis_gene_list <- getGOgeneSet(x = 'neurogenesis',OrgDb = 'org.Mm.eg.db',ont = 'BP',keytype = 'SYMBOL')
cell_adhesion_gene_list <- getGOgeneSet(x = 'cell adhesion',OrgDb = 'org.Mm.eg.db',ont = 'BP',keytype = 'SYMBOL')
cell_differentiation_gene_list <- getGOgeneSet(x = 'cell differentiation',OrgDb = 'org.Mm.eg.db',ont = 'BP',keytype = 'SYMBOL')
export_from_cell_gene_list <- getGOgeneSet(x = 'export from cell',OrgDb = 'org.Mm.eg.db',ont = 'BP',keytype = 'SYMBOL')

temp <- mouse_RNA_seurat[,mouse_RNA_seurat$New_cellType == 'Apical progenitors']

temp <- My_add_module_score(seu.obj = temp,assay = 'RNA',features = cell_adhesion_gene_list,meta_var = 'cell_adhesion',scale = FALSE,center = FALSE)
temp <- My_add_module_score(seu.obj = temp,assay = 'RNA',features = anatomical_structure_development_gene_list,meta_var = 'anatomical_structure_development',scale = FALSE,center = FALSE)
temp <- My_add_module_score(seu.obj = temp,assay = 'RNA',features = localization_gene_list,meta_var = 'localization',scale = FALSE,center = FALSE)
temp <- My_add_module_score(seu.obj = temp,assay = 'RNA',features = neurogenesis_gene_list,meta_var = 'neurogenesis',scale = FALSE,center = FALSE)
temp <- My_add_module_score(seu.obj = temp,assay = 'RNA',features = cell_differentiation_gene_list,meta_var = 'cell_differentiation',scale = FALSE,center = FALSE)
temp <- My_add_module_score(seu.obj = temp,assay = 'RNA',features = export_from_cell_gene_list,meta_var = 'export_from_cell',scale = FALSE,center = FALSE)

temp <- temp@meta.data

p1 <- ggplot(temp,aes(x=Age,y=cell_adhesion,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'cell adhesion GO list') + xlab('Age') + ylab('Module score')

p2 <- ggplot(temp,aes(x=Age,y=anatomical_structure_development,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'anatomical structure development GO list') + xlab('Age') + ylab('Module score')

p3 <- ggplot(temp,aes(x=Age,y=localization,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'localization GO list') + xlab('Age') + ylab('Module score')

p4 <- ggplot(temp,aes(x=Age,y=neurogenesis,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'neurogenesis GO list') + xlab('Age') + ylab('Module score')

p5 <- ggplot(temp,aes(x=Age,y=cell_differentiation,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'cell differentiation GO list') + xlab('Age') + ylab('Module score')

p6 <- ggplot(temp,aes(x=Age,y=export_from_cell,fill=Phase)) + 
  geom_boxplot(position = 'dodge2',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('G1','S','G2M')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'export from cell GO list') + xlab('Age') + ylab('Module score')

pdf(file = './res/step_5_fig_210727/mouse_RG_GO_list_module_score_boxplot.pdf',width = 18,height = 10)
p1+p2+p3+p4+p5+p6+plot_layout(nrow = 2)
dev.off()