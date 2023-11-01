#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: using monocle2 to investigate Excite Neuron                     ##
## Data: 2021.08.05                                                                ##
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
library(ggsci)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')


# load data ---------------------------------------------------------------
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
#subset the Ex
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('IP','Ex','Ex-U','Ex-SP')]
p1 <- DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)

p1+p2

# sub cell type -----------------------------------------------------------
my_temp_function <- function(seu.obj,cluster,cell_type){
  seu.obj@meta.data[seu.obj$seurat_clusters %in% cluster,"sub_cell_type"] <- cell_type
  return(seu.obj)
}

macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 30,cell_type = 'Ex-1')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 9,cell_type = 'Ex-2')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 1,cell_type = 'Ex-3')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 4,cell_type = 'Ex-4')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 12,cell_type = 'Ex-5')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 22,cell_type = 'Ex-6')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 2,cell_type = 'Ex-7')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 14,cell_type = 'Ex-8')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 8,cell_type = 'Ex-9')

macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 16,cell_type = 'Ex-U-1')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 15,cell_type = 'Ex-U-2')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 24,cell_type = 'Ex-U-3')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 6,cell_type = 'Ex-U-4')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 18,cell_type = 'Ex-U-5')

macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 19,cell_type = 'Ex-SP-1')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 20,cell_type = 'Ex-SP-2')

DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)


# monocle2 for all Ex -----------------------------------------------------
## create celldataset object -----------------------------------------------
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


## find variable gene ------------------------------------------------------
macaque_RNA_monocle <- detectGenes(macaque_RNA_monocle,min_expr = 3)
head(fData(macaque_RNA_monocle))
expressed_genes <- row.names(subset(fData(macaque_RNA_monocle),num_cells_expressed >= 50))

diff_test_res <- differentialGeneTest(macaque_RNA_monocle[expressed_genes,],fullModelFormulaStr = "~cell_type+donor+nCount_RNA",reducedModelFormulaStr = "~donor+nCount_RNA")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
macaque_RNA_monocle <- setOrderingFilter(macaque_RNA_monocle,ordering_genes = ordering_genes)
plot_ordering_genes(macaque_RNA_monocle)

## reduce dimension --------------------------------------------------------
macaque_RNA_monocle <- reduceDimension(macaque_RNA_monocle, max_components = 2, method = 'DDRTree',verbose = TRUE,norm_method = 'log')

## visulize ----------------------------------------------------------------
macaque_RNA_monocle <- orderCells(macaque_RNA_monocle)

plot_cell_trajectory(macaque_RNA_monocle, color_by = "cell_type",show_branch_points = FALSE,cell_size = 0.5) + 
  theme_classic() + 
  theme(legend.position = 'right',
        aspect.ratio = 1,
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(family = 'Times',face = 'bold',size = 12)) + 
  guides(colour = guide_legend(override.aes = list(shape = 15,size = 5))) + 
  xlab('Dimension 1') + ylab('Dimension 2') + labs(colour = 'Cell type')

#maybe a dead end

# layer marker ------------------------------------------------------------
temp <- macaque_RNA_seurat[,macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_1'] > -6 & macaque_RNA_seurat@reductions$umap@cell.embeddings[,'UMAP_2'] < 10]
DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE)

## human marker ------------------------------------------------------------
CPo_marker <- c('MPC1','DACT1','GALNT3','MAN1C1','MYCN','RHOU','THG1L','TMOD1')
CPi_marker <- c('LMO4','CD1D','CYTIP','FOXP1','HS3ST2','PCDH19','PTPRK','SAMD5','TRAF5')
SP_marker <- c('ADTRP','CDH18','HAS3','HS3ST4','NXPH3','SLC35F3','TMEM178A','ZFPM2')
IZ_marker <- c('ADM','COL20A1','NEU4','ENSMMUG00000056728','PDGFRA','PPP1R3G','SP5')
OSVZ_marker <- c('ACSBG1','BTBD17','PPP1R17','LYN','PLAC9','PLAGL1','PREX2','TNC','EOMES')
ISVZ_marker <- c('PAX6','CDK6','GAS1','JAG1','PTH2R','SETD7')
VZ_marker <- c('ST13','MGARP','CA3','CASP6','CD44','GLT8D2','PGAM2','POU4F1','SFRP2','GFAP')

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_human_CPo_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = CPo_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_human_CPi_marker_featureplot.pdf',width = 12,height = 9)
geneSetAveragePlot(genes = CPi_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE,num.panel.rows = 3)
dev.off()

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_human_SP_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = SP_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_human_IZ_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = IZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_human_OSVZ_marker_featureplot.pdf',width = 12,height = 9)
geneSetAveragePlot(genes = OSVZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_human_ISVZ_marker_featureplot.pdf',width = 12,height = 7)
geneSetAveragePlot(genes = ISVZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_human_VZ_marker_featureplot.pdf',width = 12,height = 9)
geneSetAveragePlot(genes = VZ_marker,object = temp,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE)
dev.off()


## macaque bulk layer score ------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('End','Per','Astrocyte','Mic'))]
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

meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")]

# meta_data <- t(base::scale(t(meta_data)))

#create annotation
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

cell_type_annotation <- HeatmapAnnotation(cell_type=macaque_RNA_seurat$cell_type,
                                          col = list(cell_type = temp_col[names(table(macaque_RNA_seurat$cell_type))]),
                                          show_annotation_name = TRUE,which = 'row',name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

#plot
col_fun <- colorRamp2(c(0,0.5,1),c("#440154FF","#21908CFF","#FDE725FF"))

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_macaque_bulk_module_heatmap.pdf',width = 8,height = 8)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(macaque_RNA_seurat$cell_type,
                           levels = c('InPSB','InMGE','InCGE','Ex-SP','Ex-U','Ex','IP','Cycling','oRG','vRG','OPC')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(6,units = 'inches'),width = unit(4,units = 'inches'),
        name = 'expression',border = TRUE,col = col_fun)
dev.off()


### only Ex -----------------------------------------------------------------
# color paramater
# color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
# temp_col <- data.frame(id=c('Ex-1','Ex-2','Ex-3','Ex-4','Ex-5','Ex-6','Ex-7','Ex-8','Ex-9'),
#                        col=c("#82C9BB","#74C1B4","#65B9AE","#57B1A7","#48A9A1","#39A19A","#2B9993","#1C918D","#0E8986"))
# color_paramater <- rbind(color_paramater,temp_col)
# temp_col <- data.frame(id=c('Ex-U-1','Ex-U-2','Ex-U-3','Ex-U-4','Ex-U-5'),
#                        col=c("#007881","#006E83","#006485","#005A87","#005089"))
# color_paramater <- rbind(color_paramater,temp_col)
# temp_col <- data.frame(id=c('Ex-SP-1','Ex-SP-2'),
#                        col=c("#002E5C","#00172E"))
# color_paramater <- rbind(color_paramater,temp_col)
# 
# write.csv(color_paramater,file = './data/parameter/color_paramater.csv')
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id


macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('IP','Ex','Ex-U','Ex-SP')]
my_temp_function <- function(seu.obj,cluster,cell_type){
  seu.obj@meta.data[seu.obj$seurat_clusters %in% cluster,"sub_cell_type"] <- cell_type
  return(seu.obj)
}

macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 30,cell_type = 'Ex-1')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 9,cell_type = 'Ex-2')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 1,cell_type = 'Ex-3')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 4,cell_type = 'Ex-4')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 12,cell_type = 'Ex-5')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 22,cell_type = 'Ex-6')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 2,cell_type = 'Ex-7')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 14,cell_type = 'Ex-8')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 8,cell_type = 'Ex-9')

macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 16,cell_type = 'Ex-U-1')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 15,cell_type = 'Ex-U-2')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 24,cell_type = 'Ex-U-3')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 6,cell_type = 'Ex-U-4')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 18,cell_type = 'Ex-U-5')

macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 19,cell_type = 'Ex-SP-1')
macaque_RNA_seurat <- my_temp_function(macaque_RNA_seurat,cluster = 20,cell_type = 'Ex-SP-2')
meta_data <- macaque_RNA_seurat@meta.data
meta_data <- meta_data[,c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")]

cell_type_annotation <- HeatmapAnnotation(cell_type=macaque_RNA_seurat$sub_cell_type,
                                          col=list(cell_type=temp_col[c('IP','Ex-1','Ex-2','Ex-3','Ex-4','Ex-5','Ex-6','Ex-7','Ex-8','Ex-9','Ex-U-1','Ex-U-2','Ex-U-3','Ex-U-4','Ex-U-5','Ex-SP-1','Ex-SP-2')]),
                                          show_annotation_name = TRUE,which = 'row',name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

pdf(file = './res/step_6_fig_210805/macaque_RNA_seurat_macaque_bulk_module_heatmap_only_Ex.pdf',width = 8,height = 8)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(macaque_RNA_seurat$sub_cell_type,
                           levels = c('IP','Ex-1','Ex-2','Ex-3','Ex-4','Ex-5','Ex-6','Ex-7','Ex-8','Ex-9','Ex-U-1','Ex-U-2','Ex-U-3','Ex-U-4','Ex-U-5','Ex-SP-1','Ex-SP-2')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(6,units = 'inches'),width = unit(4,units = 'inches'),
        name = 'expression',border = TRUE,col = col_fun)
dev.off()

# try to annotate Ex sub cell type ----------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
DimPlot(macaque_RNA_seurat,group.by = 'sub_cell_type',label = TRUE,repel = TRUE)
mouse_RNA_seurat <- readRDS(file = './data/public/Molecular_logic_of_cellular_diversification_in_the_mouse_cerebral_cortex/mouse_RNA_seurat.rds')
p1 <- DimPlot(mouse_RNA_seurat,group.by = 'New_cellType',label = FALSE)
p2 <- DimPlot(mouse_RNA_seurat,label = FALSE,cells.highlight = colnames(mouse_RNA_seurat[,mouse_RNA_seurat$New_cellType %in% c('Immature neurons')]))
p1+p2
human_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
greenleaf_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
DimPlot(greenleaf_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'UMAP')
## label transfer ----------------------------------------------------------
### human_RNA_seurat --------------------------------------------------------

#convert gene
temp <- macaque_RNA_seurat@assays$RNA@counts
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv',header = TRUE)
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp
DefaultAssay(macaque_RNA_seurat) <- 'converted'

#pre-process
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'converted',nfeatures = 5000,npcs = 50,preprocess = TRUE)
human_RNA_seurat <- my_process_seurat(object = human_RNA_seurat,assay = 'RNA',nfeatures = 5000,npcs = 50,preprocess = TRUE)

neuron_anchor <- FindTransferAnchors(reference = human_RNA_seurat,query = macaque_RNA_seurat,reference.assay = 'RNA',query.assay = 'converted',dims = 1:30)
neuron_prediction <- TransferData(anchorset = neuron_anchor,refdata = human_RNA_seurat$Cluster,dims = 1:30)
macaque_RNA_seurat <- AddMetaData(object = macaque_RNA_seurat,metadata = neuron_prediction)
table(macaque_RNA_seurat$predicted.id)
scibet::Confusion_heatmap(ori = macaque_RNA_seurat$cell_type,prd = macaque_RNA_seurat$predicted.id)

DimPlot(macaque_RNA_seurat,group.by = 'predicted.id',label = TRUE,reduction = 'umap',repel = TRUE)
macaque_RNA_seurat$human_prediction <- macaque_RNA_seurat$predicted.id

### mouse_RNA_seurat --------------------------------------------------------
#convert gene
temp <- mouse_RNA_seurat@assays$RNA@counts
mouse_to_macaque_anno <- read.csv(file = './data/reference/GRCm39_to_Mmul_10.csv',header = TRUE)
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = mouse_to_macaque_anno,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
mouse_RNA_seurat[['converted']] <- temp

#pre-process
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',nfeatures = 5000,npcs = 50,preprocess = TRUE)
mouse_RNA_seurat <- my_process_seurat(object = mouse_RNA_seurat,assay = 'converted')