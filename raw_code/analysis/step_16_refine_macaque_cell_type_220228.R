#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: refine macaque cell type                                        ##
## Data: 2022.02.28                                                                ##
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
library(ArchR)
library(parallel)
library(ggrepel)
library(org.Mmu.eg.db)
library(org.Hs.eg.db)
library(topGO)
library(ggsci)
library(scales)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#load data
macaque_integration_Seurat <- readRDS(file = './processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_220226.rds')
DimPlot(macaque_integration_Seurat,group.by = 'cell_type',split.by = 'tech',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP') + theme(aspect.ratio = 1)
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')
macaque_multiome_ArchR$RNA_cell_type <- macaque_integration_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]

# refine RG from gene expression---------------------------------------------------------------
## RG-3 --------------------------------------------------------------------
#dimplot
temp_RNA <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']@reductions$harmonyUMAP@cell.embeddings
p1 <- my_dimplot(embedding = temp_RNA,meta_data = macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']@meta.data,
                 group.by = 'cell_type') + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1.5),
        axis.title = element_blank())
p1
temp_ATAC <- as.data.frame(macaque_multiome_ArchR@embeddings$UMAP@listData$df)
p2 <- my_dimplot(embedding = temp_ATAC,meta_data = macaque_multiome_ArchR@cellColData,
                 group.by = 'RNA_cell_type') + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1.5),
        axis.title = element_blank())
p2

pdf(file = './res/step_16_fig_220228/macaque_RNA_ATAC_dimplot.pdf',width = 10,height = 5)
p1+p2+plot_layout(ncol = 2)
dev.off()

p1 <- my_dimplot(embedding = temp_RNA,meta_data = macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']@meta.data,
                 group.by = 'cell_type',
                 cells.highlight = colnames(macaque_integration_Seurat)[macaque_integration_Seurat$cell_type == 'RG-3']) + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1.5),
        axis.title = element_blank())
p1
p2 <- my_dimplot(embedding = temp_ATAC,meta_data = macaque_multiome_ArchR@cellColData,
                 group.by = 'RNA_cell_type',
                 cells.highlight = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$RNA_cell_type == 'RG-3']) + 
  theme(aspect.ratio = 1,
        legend.position = 'right',
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1.5),
        axis.title = element_blank())
p2

pdf(file = './res/step_16_fig_220228/macaque_RNA_ATAC_dimplot_highlight_RG_3.pdf',width = 10,height = 5)
p1+p2+plot_layout(ncol = 2)
dev.off()

#gene expression
pdf(file = './res/step_16_fig_220228/macaque_multiome_gene_express_vlnplot_show_RG_3.pdf',width = 12,height = 8)
VlnPlot(macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome'],group.by = 'cell_type',
        features = c('SOX10','EGFR','DLX5','DLX2','GAD1','GAD2'),
        pt.size = 0,assay = 'RNA',slot = 'data',ncol = 3)
dev.off()

pdf(file = './res/step_16_fig_220228/macaque_multiome_gene_express_feature_plot_show_RG_3.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = c('SOX10','EGFR','DLX5','DLX2','GAD1','GAD2'),object = macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome'],
                   object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',reduction_key = 'UMAP_',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,
                   trim = c(0,0.99),print = TRUE,num.panel.rows = 2)
dev.off()

#find marker
RG_3_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-3',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_3_marker_list <- RG_3_marker_list[order(RG_3_marker_list$pct.1,decreasing = TRUE),]
RG_3_marker_list <- RG_3_marker_list[RG_3_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_3_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

RG_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-3',ident.2 = c('RG-2','RG-1'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_subset_marker <- RG_subset_marker[order(RG_subset_marker$pct.1,decreasing = TRUE),]
RG_subset_marker <- RG_subset_marker[RG_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect RG-3 marker
gene_list <- dplyr::intersect(rownames(RG_3_marker_list),rownames(RG_subset_marker))
VlnPlot(object = macaque_integration_Seurat,features = gene_list[13:24],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#save RG-3 marker
saveRDS(gene_list,file = './res/step_16_fig_220228/macaque_integration_RG_3_marker_list.rds')

## RG-1 --------------------------------------------------------------------
RG_1_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-1',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_1_marker_list <- RG_1_marker_list[order(RG_1_marker_list$pct.1,decreasing = TRUE),]
RG_1_marker_list <- RG_1_marker_list[RG_1_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_1_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

RG_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-1',ident.2 = c('RG-2','RG-3'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_subset_marker <- RG_subset_marker[order(RG_subset_marker$pct.1,decreasing = TRUE),]
RG_subset_marker <- RG_subset_marker[RG_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect RG-1 marker
gene_list <- dplyr::intersect(rownames(RG_1_marker_list),rownames(RG_subset_marker))
VlnPlot(object = macaque_integration_Seurat,features = gene_list[13:24],assay = 'RNA',group.by = 'cell_type',pt.size = 0,slot = 'data')

#save RG-1 marker
saveRDS(gene_list,file = './res/step_16_fig_220228/macaque_integration_RG_1_marker_list.rds')

#Go enrichment
temp <- gene_list
gene_list <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- c(gene_list %in% temp)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

#BP
GO_enrich <- new("topGOdata",
                 description = "RG-1",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Mmu.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$Term[1] <- 'regulation of cellular response to growth factor stimulus'

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = allRes$Term)

p1 <- ggplot(allRes[1:5,],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#FF7F0EFF',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.9,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('RG-1 enriched BP')

pdf(file = './res/step_16_fig_220228/macaque_integration_RG_1_marker_list_GO_BP_barplot.pdf',width = 6,height = 3)
p1
dev.off()

#BMP signal gene
BMP_gene <- getGOgeneSet(x = 'GO:0030509',OrgDb = 'org.Mmu.eg.db',ont = 'BP',keytype = 'SYMBOL')
gene_list <- dplyr::intersect(rownames(RG_1_marker_list),rownames(RG_subset_marker))
gene_list <- dplyr::intersect(gene_list,BMP_gene)

pdf(file = './res/step_16_fig_220228/macaque_integration_RG_1_marker_list_BMP_gene_vlnplot.pdf',width = 12,height = 4)
VlnPlot(macaque_integration_Seurat,features = gene_list,pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data') + 
  theme(aspect.ratio = 0.5)
dev.off()

## RG-2 --------------------------------------------------------------------
RG_2_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-2',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_2_marker_list <- RG_2_marker_list[order(RG_2_marker_list$pct.1,decreasing = TRUE),]
RG_2_marker_list <- RG_2_marker_list[RG_2_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_2_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

RG_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-2',ident.2 = c('RG-1','RG-3'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_subset_marker <- RG_subset_marker[order(RG_subset_marker$pct.1,decreasing = TRUE),]
RG_subset_marker <- RG_subset_marker[RG_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect RG-2 marker
gene_list <- dplyr::intersect(rownames(RG_2_marker_list),rownames(RG_subset_marker))
VlnPlot(object = macaque_integration_Seurat,features = gene_list[1:12],assay = 'RNA',group.by = 'cell_type',pt.size = 0,slot = 'data')

#save RG-2 marker
saveRDS(gene_list,file = './res/step_16_fig_220228/macaque_integration_RG_2_marker_list.rds')

#Go enrichment
temp <- gene_list
gene_list <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- c(gene_list %in% temp)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

#BP
GO_enrich <- new("topGOdata",
                 description = "RG-2",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Mmu.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$Term[3] <- 'regulation of amyloid precursor protein catabolic process'

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = allRes$Term)

p1 <- ggplot(allRes[1:5,],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#FF7F0EFF',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.9,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('RG-2 enriched BP')

pdf(file = './res/step_16_fig_220228/macaque_integration_RG_2_marker_list_GO_BP_barplot.pdf',width = 6,height = 3)
p1
dev.off()

#amyloid formation gene
amyloid_formation_gene <- getGOgeneSet(x = 'GO:0034205',OrgDb = 'org.Mmu.eg.db',ont = 'BP',keytype = 'SYMBOL')
gene_list <- dplyr::intersect(rownames(RG_2_marker_list),rownames(RG_subset_marker))
amyloid_formation_gene <- dplyr::intersect(gene_list,amyloid_formation_gene)

pdf(file = './res/step_16_fig_220228/macaque_integration_RG_2_marker_list_amyloid_formation_gene_vlnplot.pdf',width = 8,height = 4)
VlnPlot(macaque_integration_Seurat,features = amyloid_formation_gene,pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data') + 
  theme(aspect.ratio = 0.5)
dev.off()

# refine Ex neuron --------------------------------------------------------
DimPlot(macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

dotplot_matrix <- my_dotplot(macaque_integration_Seurat,assay = 'RNA', 
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
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Ependymal','Mic','Per','End'))
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

## Ex-4 --------------------------------------------------------------------
#Ex-4 enrich SP marker
macaque_annotation <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[macaque_annotation$type == 'gene',]
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]
layer_marker <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')

#layer marker
temp <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome' & !(macaque_integration_Seurat$cell_type %in% c('End','Per','Mic','Ependymal','OPC'))]
temp@meta.data[temp$cell_type %in% c('Cyc-S','Cyc-G2M'),"cell_type"] <- 'Cycling'
temp@meta.data[temp$cell_type %in% c('RG-1','RG-2','RG-3'),"cell_type"] <- 'RG'

my_temp_function <- function(layer_marker,layer){
  marker <- layer_marker[[layer]]
  marker <- unique(marker)
  print('marker in annotation list')
  print(table(marker %in% macaque_annotation$gene_id))
  marker <- marker[marker %in% macaque_annotation$gene_id]
  print('marker trans id to name')
  marker <- macaque_annotation$gene_name[base::match(marker,table = macaque_annotation$gene_id)]
  print('marker in express matrix')
  print(table(marker %in% rownames(temp)))
  marker <- marker[marker %in% rownames(temp)]
  print('duplicated markers')
  print(table(duplicated(marker)))
  marker <- unique(marker)
  return(marker)
}

names(layer_marker)

for (i in c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")){
  temp_marker <- my_temp_function(layer_marker = layer_marker,layer = i)
  temp <- My_add_module_score(seu.obj = temp,assay = 'RNA',features = temp_marker,meta_var = i,scale = FALSE,center = FALSE)
}

meta_data <- temp@meta.data
meta_data <- meta_data[,c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")]

#create annotation
cell_type_annotation <- HeatmapAnnotation(cell_type=temp$cell_type,show_annotation_name = TRUE,
                                          which = 'row',name = 'cell type',annotation_label = 'cell type')

depth_annotation <- HeatmapAnnotation(depth = 1:8,col = list(depth = colorRamp2(c(1,8),c('lightgrey','black'))))

#heatmap plot
col_fun <- colorRamp2(c(0,0.5,1),c("#440154FF","#21908CFF","#FDE725FF"))

pdf(file = './res/step_16_fig_220228/macaque_multiome_gene_express_layer_marker_heatmap.pdf',width = 8,height = 8)
Heatmap(meta_data,cluster_columns = FALSE,cluster_rows = FALSE,show_column_names = TRUE,show_row_names = FALSE,
        row_split = factor(temp$cell_type,
                           levels = c('InMGE','InCGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cycling','RG')),
        column_split = factor(c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ"),
                              levels = c("MZ","CPo","CPi","SP","IZ","OSVZ","ISVZ","VZ")),
        top_annotation = depth_annotation,left_annotation = cell_type_annotation,
        height = unit(6,units = 'inches'),width = unit(4,units = 'inches'),
        name = 'expression',border = TRUE,col = col_fun)
dev.off()

#feature plot
temp_marker <- my_temp_function(layer_marker = layer_marker,layer = 'SP')
gene_list <- c('ADTRP','CDH18','HAS3','HS3ST4','NXPH3','SLC35F3','TMEM178A','ZFPM2')
gene_list <- gene_list[gene_list %in% rownames(macaque_integration_Seurat)]

pdf(file = './res/step_16_fig_220228/macaque_integration_Seurat_SP_marker_express_featureplot.pdf',width = 12,height = 6)
geneSetAveragePlot(genes = gene_list,object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.01,0.99),print = TRUE,
                   num.panel.rows = 2,plot.title = 'SP classic marker')
dev.off()

#find marker
Ex_4_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-4',group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_4_marker_list <- Ex_4_marker_list[order(Ex_4_marker_list$pct.1,decreasing = TRUE),]
Ex_4_marker_list <- Ex_4_marker_list[Ex_4_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_4_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

Ex_subset_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-4',ident.2 = c('Ex-1','Ex-2','Ex-3'),group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                     logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_subset_marker_list <- Ex_subset_marker_list[order(Ex_subset_marker_list$pct.1,decreasing = TRUE),]
Ex_subset_marker_list <- Ex_subset_marker_list[Ex_subset_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_subset_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#intersect Ex-4 marker
gene_list <- dplyr::intersect(rownames(Ex_4_marker_list),rownames(Ex_subset_marker_list))
VlnPlot(macaque_integration_Seurat,features = gene_list[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#save Ex_4 marker gene
saveRDS(gene_list,file = './res/step_16_fig_220228/macaque_integration_Ex_4_marker_list.rds')

#Go enrichment
temp <- gene_list
gene_list <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- c(gene_list %in% temp)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

#BP
GO_enrich <- new("topGOdata",
                 description = "Ex-4",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Mmu.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$Term[1] <- 'regulation of small GTPase mediated signal transduction'

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = allRes$Term)

p1 <- ggplot(allRes[1:5,],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#FF7F0EFF',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.9,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('Ex-4 enriched BP')

pdf(file = './res/step_16_fig_220228/macaque_integration_Ex_4_marker_list_GO_BP_barplot.pdf',width = 6,height = 3)
p1
dev.off()

## Ex-3 --------------------------------------------------------------------
#CP marker feature plot
gene_list <- c('MPC1','DACT1','GALNT3','MAN1C1','MYCN','RHOU','THG1L','TMOD1','LMO4')
pdf(file = './res/step_16_fig_220228/macaque_integration_Seurat_CPo_marker_express_featureplot.pdf')
geneSetAveragePlot(genes = gene_list,object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.01,0.99),print = TRUE,
                   num.panel.rows = 2,plot.title = 'CPo classic marker')
dev.off()

gene_list <- c('CD1D','CYTIP','FOXP1','HS3ST2','PCDH19','PTPRK','SAMD5','TRAF5')
pdf(file = './res/step_16_fig_220228/macaque_integration_Seurat_CPi_marker_express_featureplot.pdf')
geneSetAveragePlot(genes = gene_list,object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.01,0.99),print = TRUE,
                   num.panel.rows = 2,plot.title = 'CPi classic marker')
dev.off()

#find marker
Ex_3_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-3',group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_3_marker_list <- Ex_3_marker_list[order(Ex_3_marker_list$pct.1,decreasing = TRUE),]
Ex_3_marker_list <- Ex_3_marker_list[Ex_3_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_3_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

Ex_subset_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-3',ident.2 = c('Ex-1','Ex-2','Ex-4'),group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                     logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_subset_marker_list <- Ex_subset_marker_list[order(Ex_subset_marker_list$pct.1,decreasing = TRUE),]
Ex_subset_marker_list <- Ex_subset_marker_list[Ex_subset_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_subset_marker_list)[13:24],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#intersect Ex-3 marker
gene_list <- dplyr::intersect(rownames(Ex_3_marker_list),rownames(Ex_subset_marker_list))
VlnPlot(macaque_integration_Seurat,features = gene_list[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#save Ex-3 marker
saveRDS(gene_list,file = './res/step_16_fig_220228/macaque_integration_Ex_3_marker_list.rds')

#Go enrichment
temp <- gene_list
gene_list <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- c(gene_list %in% temp)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

#BP
GO_enrich <- new("topGOdata",
                 description = "Ex-3",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Mmu.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = allRes$Term)

p1 <- ggplot(allRes[1:5,],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#FF7F0EFF',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.9,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('Ex-3 enriched BP')

pdf(file = './res/step_16_fig_220228/macaque_integration_Ex_3_marker_list_GO_BP_barplot.pdf',width = 6,height = 3)
p1
dev.off()

## Ex-2 --------------------------------------------------------------------
#find marker
Ex_2_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-2',group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_2_marker_list <- Ex_2_marker_list[order(Ex_2_marker_list$pct.1,decreasing = TRUE),]
Ex_2_marker_list <- Ex_2_marker_list[Ex_2_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_2_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

Ex_subset_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-2',ident.2 = c('Ex-1','Ex-3','Ex-4'),group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                     logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_subset_marker_list <- Ex_subset_marker_list[order(Ex_subset_marker_list$pct.1,decreasing = TRUE),]
Ex_subset_marker_list <- Ex_subset_marker_list[Ex_subset_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_subset_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#intersect Ex-2 marker
gene_list <- dplyr::intersect(rownames(Ex_2_marker_list),rownames(Ex_subset_marker_list))
VlnPlot(macaque_integration_Seurat,features = gene_list[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#save gene list
saveRDS(gene_list,file = './res/step_16_fig_220228/macaque_integration_Ex_2_marker_list.rds')

#feature plot
geneSetAveragePlot(genes = gene_list[1:12],object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0,0.99),print = TRUE,num.panel.rows = 3)

#cluster?
DimPlot(object = macaque_integration_Seurat,group.by = c('cell_type','seurat_clusters'),label = TRUE,repel = TRUE)
#can not see clearly cluster

#merge Ex-1 and Ex-2
macaque_integration_Seurat <- readRDS(file = './processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_220226.rds')
macaque_integration_Seurat$old_cell_type <- macaque_integration_Seurat$cell_type
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$cell_type %in% c('Ex-1','Ex-2'),"cell_type"] <- 'Ex-1'
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$cell_type %in% c('Ex-3'),"cell_type"] <- 'Ex-2'
macaque_integration_Seurat@meta.data[macaque_integration_Seurat$cell_type %in% c('Ex-4'),"cell_type"] <- 'Ex-3'
DimPlot(macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
saveRDS(macaque_integration_Seurat,file = './processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_annotated_220302.rds')

## Ex-1 --------------------------------------------------------------------
#find marker
Ex_1_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-1',group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_1_marker_list <- Ex_1_marker_list[order(Ex_1_marker_list$pct.1,decreasing = TRUE),]
Ex_1_marker_list <- Ex_1_marker_list[Ex_1_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_1_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

Ex_subset_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-1',ident.2 = c('Ex-2','Ex-3'),group.by = 'cell_type',assay = 'RNA',slot = 'data',
                                     logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,verbose = TRUE,only.pos = TRUE)
Ex_subset_marker_list <- Ex_subset_marker_list[order(Ex_subset_marker_list$pct.1,decreasing = TRUE),]
Ex_subset_marker_list <- Ex_subset_marker_list[Ex_subset_marker_list$pct.2 < 0.4,]
VlnPlot(macaque_integration_Seurat,features = rownames(Ex_subset_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#intersect Ex-1 marker
gene_list <- dplyr::intersect(rownames(Ex_1_marker_list),rownames(Ex_subset_marker_list))
VlnPlot(macaque_integration_Seurat,features = gene_list[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')

#save gene list
save(gene_list,file = './res/step_16_fig_220228/macaque_integration_Ex_1_marker_list.rds')
