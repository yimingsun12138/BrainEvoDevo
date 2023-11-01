#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: investigate InPSB                                               ##
## Data: 2021.07.02                                                                ##
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
library(DESeq2)
library(rliger)
library(patchwork)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
#InPSB
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

#SP marker heatmap
SP_classic_marker <- c('ADTRP','CDH18','HAS3','HS3ST4','NXPH3','SLC35F3','TMEM178A','ZFPM2')
table(SP_classic_marker %in% rownames(macaque_RNA_seurat))

pdf(file = './res/step_2_fig_210702/InPSB_SP_classic_marker_vlnplot.pdf',width = 15,height = 10)
VlnPlot(macaque_RNA_seurat,features = SP_classic_marker,group.by = 'cell_type',assay = 'RNA',slot = 'data')
dev.off()

#SP layer marker by bulk RNAsea
layer_marker <- readRDS(file = './res/step_3_fig_210624/macaque_layer_DE_list_sunym_210628.rds')

macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

SP_marker <- layer_marker[['SP']]
SP_marker <- unique(SP_marker)
table(SP_marker %in% macaque_annotation$gene_id)
SP_marker <- SP_marker[SP_marker %in% macaque_annotation$gene_id]
SP_marker <- macaque_annotation$gene_name[base::match(SP_marker,table = macaque_annotation$gene_id)]
table(SP_marker %in% rownames(macaque_RNA_seurat))
SP_marker <- SP_marker[SP_marker %in% rownames(macaque_RNA_seurat)]
table(duplicated(SP_marker))
SP_marker <- unique(SP_marker)
macaque_RNA_seurat <- My_add_module_score(macaque_RNA_seurat,features = SP_marker,meta_var = 'SP_marker',scale = FALSE,center = TRUE)

#ggplot
temp <- macaque_RNA_seurat@meta.data
temp$cell_type <- factor(temp$cell_type,levels = c('End','Per','Mic','Astrocyte','vRG','oRG','OPC','Cycling','IP',
                                                   'Ex','Ex-U','Ex-SP','InCGE','InMGE','InPSB'))
color_paramater <- read.csv(file = './data/parameter/color_paramater.csv',row.names = 1)
temp_col <- color_paramater$col
names(temp_col) <- color_paramater$id

pdf(file = './res/step_2_fig_210702/InPSB_SP_layer_signature_boxplot.pdf',width = 5,height = 5)
ggplot(temp,aes(x=cell_type,y=SP_marker,fill=cell_type)) + 
  geom_boxplot(position = 'identity',outlier.alpha = 0) + 
  scale_fill_manual(values = temp_col[c('End','Per','Mic','Astrocyte','vRG','oRG','OPC','Cycling','IP',
                                        'Ex','Ex-U','Ex-SP','InCGE','InMGE','InPSB')]) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1),
        axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1)) + 
  ylim(c(-0.15,0.9)) + CenterTitle() + 
  xlab('cell type') + ylab('SP signature') + labs(title = 'centered log1p(mean counts)') + 
  geom_hline(yintercept = 0,linetype = 'dashed',size = 1,color = 'brown')
dev.off()

#make pseudobulk to see the correlation with layer bulk
InPSB_pseudobulk <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type == 'InPSB']@assays$RNA@counts
table(rownames(InPSB_pseudobulk) %in% macaque_annotation$gene_name)
InPSB_pseudobulk <- InPSB_pseudobulk[rownames(InPSB_pseudobulk) %in% macaque_annotation$gene_name,]
rownames(InPSB_pseudobulk) <- macaque_annotation[base::match(x=rownames(InPSB_pseudobulk),table = macaque_annotation$gene_name),"gene_id"]

layer_bulk <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')

gene_list <- dplyr::intersect(rownames(layer_bulk),rownames(InPSB_pseudobulk))
layer_bulk <- layer_bulk[gene_list,]
InPSB_pseudobulk <- InPSB_pseudobulk[gene_list,]

InPSB_pseudobulk_rep1 <- sample(colnames(InPSB_pseudobulk),size = 103,replace = FALSE)
InPSB_pseudobulk_rep2 <- sample(colnames(InPSB_pseudobulk)[!(colnames(InPSB_pseudobulk) %in% InPSB_pseudobulk_rep1)],size = 103,replace = FALSE)
InPSB_pseudobulk_rep3 <- colnames(InPSB_pseudobulk)[!(colnames(InPSB_pseudobulk) %in% c(InPSB_pseudobulk_rep1,InPSB_pseudobulk_rep2))]

InPSB_pseudobulk <- data.frame(InPSB_rep1 = rowSums(InPSB_pseudobulk[,InPSB_pseudobulk_rep1]),
                               InPSB_rep2 = rowSums(InPSB_pseudobulk[,InPSB_pseudobulk_rep2]),
                               InPSB_rep3 = rowSums(InPSB_pseudobulk[,InPSB_pseudobulk_rep3]))

layer_bulk <- cbind(layer_bulk,InPSB_pseudobulk)

#TPM transform
for(i in colnames(layer_bulk)){
  all_counts <- sum(layer_bulk[,i])
  layer_bulk[,i] <- layer_bulk[,i]/all_counts
}
layer_bulk <- log1p(layer_bulk * 10000)

cor_matrix <- cor(layer_bulk,method = 'pearson')
pdf(file = './res/step_2_fig_210702/InPSB_layer_cor.pdf',width = 8,height = 8)
Heatmap(cor_matrix,cluster_rows = TRUE,cluster_columns = TRUE,height = unit(5,units = 'inches'),width = unit(5,units = 'inches'),name = 'cor')
dev.off()

#InPSB marker
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
InPSB_marker <- FindMarkers(macaque_RNA_seurat,ident.1 = 'InPSB',group.by = 'cell_type',assay = 'RNA',
                            slot = 'data',test.use = 'bimod',only.pos = TRUE,min.diff.pct = 0.2,logfc.threshold = 0.1)
InPSB_marker <- InPSB_marker[order(InPSB_marker$pct.2,decreasing = FALSE),]
InPSB_marker <- InPSB_marker[order(InPSB_marker$p_val_adj,decreasing = FALSE),]

geneSetAveragePlot(genes = rownames(InPSB_marker)[1:9],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

geneSetAveragePlot(genes = rownames(InPSB_marker)[10:18],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

InPSB_marker <- InPSB_marker[InPSB_marker$p_val_adj < 0.01,]
InPSB_marker <- InPSB_marker[order(InPSB_marker$pct.2,decreasing = FALSE),]
InPSB_marker <- InPSB_marker[InPSB_marker$pct.1 >= 0.4,]

geneSetAveragePlot(genes = rownames(InPSB_marker)[1:9],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

geneSetAveragePlot(genes = rownames(InPSB_marker)[10:18],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

geneSetAveragePlot(genes = rownames(InPSB_marker)[19:27],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

geneSetAveragePlot(genes = rownames(InPSB_marker)[28:36],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

geneSetAveragePlot(genes = rownames(InPSB_marker)[37:45],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

#marker
#'FOXP2','FOXP4','EYA2','TSHZ1','PBX3','LHFPL6','MEIS2','PLXNC1','PTPRM','ALK'

#try more genes
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
InPSB_marker <- FindMarkers(macaque_RNA_seurat,ident.1 = 'InPSB',group.by = 'cell_type',assay = 'RNA',
                            slot = 'data',test.use = 'bimod',only.pos = TRUE,min.pct = 0,logfc.threshold = 0)

temp <- InPSB_marker[InPSB_marker$p_val_adj < 0.01,]
temp <- temp[temp$pct.1 > 0.4,]
temp <- temp[temp$pct.2 < 0.2,]

geneSetAveragePlot(genes = rownames(temp)[1:9],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

geneSetAveragePlot(genes = rownames(temp)[10:18],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

geneSetAveragePlot(genes = rownames(temp)[19:25],object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

pdf(file = './res/step_2_fig_210702/InPSB_marker_featureplot.pdf',width = 10,height = 7)
geneSetAveragePlot(genes = c('FOXP2','FOXP4','EYA2','TSHZ1','PBX3','LHFPL6','MEIS2','PLXNC1','PTPRM','ALK','ROR2'),object = macaque_RNA_seurat,object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))
dev.off()

pdf(file = './res/step_2_fig_210702/InPSB_marker_vlnplot_by_celltype.pdf',width = 15,height = 10)
VlnPlot(macaque_RNA_seurat,features = c('FOXP2','FOXP4','EYA2','TSHZ1','PBX3','LHFPL6','MEIS2','PLXNC1','PTPRM','ALK','ROR2'),
        group.by = 'cell_type',assay = 'RNA',slot = 'data')
dev.off()

#InPSB marker from In
InPSB_marker <- FindMarkers(macaque_RNA_seurat,ident.1 = 'InPSB',ident.2 = c('InCGE','InMGE'),group.by = 'cell_type',assay = 'RNA',
                            slot = 'data',test.use = 'bimod',only.pos = TRUE,min.pct = 0,logfc.threshold = 0)

temp <- InPSB_marker[InPSB_marker$p_val_adj < 0.01,]
temp <- temp[temp$pct.1 > 0.4,]
temp <- temp[temp$pct.2 < 0.2,]

geneSetAveragePlot(genes = rownames(temp)[1:9],object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE')],object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

#FOXP2,LHFPL6,EYA2,TSHZ1,FOXP4,ALK,PLXNC1,FRAS1

geneSetAveragePlot(genes = rownames(temp)[10:18],object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE')],object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

#PRKCB,TRPM3,CDIN1,HS6ST2

geneSetAveragePlot(genes = rownames(temp)[19:27],object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE')],object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

#ESRRG

geneSetAveragePlot(genes = rownames(temp)[28:36],object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE')],object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))

pdf(file = './res/step_2_fig_210702/InPSB_marker_featureplot_only_In.pdf',width = 15,height = 10)
geneSetAveragePlot(genes = c('FOXP2','LHFPL6','EYA2','TSHZ1','FOXP4','ALK','PLXNC1','FRAS1','PRKCB','TRPM3','CDIN1','HS6ST2','ESRRG'),
                   object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE')],object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 3,color.palette = c("lightgrey", "blue"))
dev.off()

temp <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE')]

temp$sub_cell_type <- paste(temp$cell_type,temp$seurat_clusters,sep = '_')


pdf(file = './res/step_2_fig_210702/InPSB_marker_vlnplot_by_celltype_only_In.pdf',width = 15,height = 10)
VlnPlot(temp,features = c('FOXP2','LHFPL6','EYA2','TSHZ1','FOXP4','ALK','PLXNC1','FRAS1','PRKCB','TRPM3','CDIN1','HS6ST2','ESRRG'),
        group.by = 'sub_cell_type',assay = 'RNA',slot = 'data')
dev.off()


#create the maker gene table of InPSB
temp <- FindMarkers(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE','InPSB')],ident.1 = 'InPSB',
                    group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',logfc.threshold = 0,min.pct = 0,only.pos = TRUE,
                    features = c("FOXP2","LHFPL6","EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","PRKCB","TRPM3",
                                 "CDIN1","HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"))

temp$InPSB_pct <- unlist(base::lapply(rownames(temp),function(g){
  temp <- sum(macaque_RNA_seurat@assays$RNA@counts[g,macaque_RNA_seurat$cell_type == 'InPSB'] > 0)/sum(macaque_RNA_seurat$cell_type == 'InPSB')
}))

temp$InMGE_pct <- unlist(base::lapply(rownames(temp),function(g){
  temp <- sum(macaque_RNA_seurat@assays$RNA@counts[g,macaque_RNA_seurat$cell_type == 'InMGE'] > 0)/sum(macaque_RNA_seurat$cell_type == 'InMGE')
}))

temp$InCGE_pct <- unlist(base::lapply(rownames(temp),function(g){
  temp <- sum(macaque_RNA_seurat@assays$RNA@counts[g,macaque_RNA_seurat$cell_type == 'InCGE'] > 0)/sum(macaque_RNA_seurat$cell_type == 'InCGE')
}))

write.csv(temp,file = './res/step_2_fig_210702/InPSB_marker_list.csv')


#InPSB marker express in macaque bulk layer RNAseq
macaque_layer_bulk <- readRDS(file = './processed_data/macaque_bulk_layer_RNAseq_sunym_210626.rds')

macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

table(rownames(macaque_layer_bulk) %in% macaque_annotation$gene_id)
macaque_to_macaque_anno <- data.frame(gene_id_1=macaque_annotation$gene_id,gene_id_2=macaque_annotation$gene_id,
                                      gene_name_1=macaque_annotation$gene_name,gene_name_2=macaque_annotation$gene_name)

macaque_layer_bulk <- My_Convert_Homology_Gene_ID(express_matrix = macaque_layer_bulk,anno = macaque_to_macaque_anno)
macaque_layer_bulk <- as.data.frame(macaque_layer_bulk)
for(i in colnames(macaque_layer_bulk)){
  macaque_layer_bulk[,i] <- (macaque_layer_bulk[,i])/sum(macaque_layer_bulk[,i])
}
macaque_layer_bulk <- log1p(macaque_layer_bulk*10000)
macaque_layer_bulk <- macaque_layer_bulk[c("FOXP2","LHFPL6","EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","PRKCB",
                                           "TRPM3","CDIN1","HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),]
# macaque_layer_bulk <- t(scale(t(macaque_layer_bulk),center = TRUE))

temp <- data.frame(sample=colnames(macaque_layer_bulk),donor=colnames(macaque_layer_bulk))
temp_order <- data.frame(layer=c('VZ','ISVZ','OSVZ','IZ','SP','CPi','CPo','MZ'),order=1:8)
rownames(temp_order) <- temp_order$layer
rownames(temp) <- temp$sample
temp$donor <- sub(pattern = '-.*$',replacement = '',fixed = FALSE,temp$sample)
temp$layer <- sub(pattern = '^.*-',replacement = '',fixed = FALSE,temp$sample)
temp$order <- temp_order[temp$layer,"order"]
temp <- temp[order(temp$donor,decreasing = FALSE),]
temp <- temp[order(temp$order,decreasing = FALSE),]

macaque_layer_bulk <- macaque_layer_bulk[,temp$sample]
row_dend = as.dendrogram(hclust(dist(macaque_layer_bulk),method = 'complete'))

pdf(file = './res/step_2_fig_210702/InPSB_marker_layer_express.pdf',width = 14,height = 7)
Heatmap(macaque_layer_bulk,cluster_rows = row_dend,cluster_columns = FALSE,show_row_names = TRUE,show_column_names = TRUE,
        height = unit(5.1,units = 'inches'),width = unit(9.6,units = 'inches'),name = 'normalized express')
dev.off()

macaque_layer_bulk <- t(scale(t(macaque_layer_bulk),center = TRUE))
pdf(file = './res/step_2_fig_210702/InPSB_marker_layer_express_scaled.pdf',width = 14,height = 7)
Heatmap(macaque_layer_bulk,cluster_rows = row_dend,cluster_columns = FALSE,show_row_names = TRUE,show_column_names = TRUE,
        height = unit(5.1,units = 'inches'),width = unit(9.6,units = 'inches'),name = 'normalized express')
dev.off()

# In marker express in cell types
gene_list <- list(InPSB=c("FOXP2","LHFPL6","EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","PRKCB","TRPM3","CDIN1","HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
                  In=c('DLX1',"DLX2","DLX5",'DLX6','GAD1','GAD2'),
                  InMGE=c('LHX6','SST'),
                  InCGE=c('SP8','NR2F2'))

pdf(file = './res/step_2_fig_210702/InPSB_marker_express_cell_type_210719.pdf',width = 10,height = 5)
my_dotplot(object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE','InPSB')],
           assay = 'RNA',features = gene_list,scale = FALSE,group.by = 'cell_type',return_data_plot = FALSE) + 
  theme_classic() + 
  theme(axis.line = element_blank(),
        panel.border = element_rect(fill = NA,colour = 'black',size = 1),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))
dev.off()

express_matrix_raw <- my_dotplot(object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE','InPSB')],
                                 assay = 'RNA',features = gene_list,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
write.csv(express_matrix_raw,file = './res/step_2_fig_210702/InPSB_marker_express_cell_type_210719.csv')

# new marker
# MEIS2,ETV1,SP8
c('MEIS2','ETV1','SP8') %in% rownames(macaque_RNA_seurat)

pdf(file = './res/step_2_fig_210702/InPSB_featureplot_210723.pdf',width = 15,height = 7)
geneSetAveragePlot(genes = c('MEIS2','ETV1','SP8'),
                   object = macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InMGE','InCGE')],object.class = 'seurat',
                   embedding = 'umap',plot.type = 'panels',scaled = FALSE,trim = NULL,lims = NULL,
                   aspectratio = 1,print = TRUE,num.panel.rows = 1,color.palette = c("lightgrey", "blue"))
dev.off()

pdf(file = './res/step_2_fig_210702/InPSB_vlnplot_210723.pdf',width = 15,height = 7)
VlnPlot(macaque_RNA_seurat,group.by = 'cell_type',features = c('MEIS2','ETV1','SP8')) + 
  theme(aspect.ratio = 1)
dev.off()
#public data to find InPSB

##########################
##compare with greenleaf##
##########################

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
human_RNA_seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat.rds')
DefaultAssay(macaque_RNA_seurat) <- 'RNA'
DefaultAssay(human_RNA_seurat) <- 'originalexp'

#settle the converted gene symbol assay
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')
temp <- human_RNA_seurat@assays$originalexp@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
human_RNA_seurat[['converted']] <- temp
rm(temp)
DefaultAssay(human_RNA_seurat) <- 'converted'
human_RNA_seurat <- NormalizeData(human_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
DefaultAssay(human_RNA_seurat) <- 'originalexp'
human_RNA_seurat <- NormalizeData(human_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
human_RNA_seurat$greenleaf_cluster <- human_RNA_seurat$seurat_clusters
#all settled

#In group in greenleaf data
DimPlot(human_RNA_seurat,group.by = 'greenleaf_cluster',label = TRUE,repel = TRUE)
VlnPlot(human_RNA_seurat,group.by = 'greenleaf_cluster',assay = 'converted',slot = 'data',
        features = c('DLX2','DLX5','GAD2','LHX6','SST','SP8','NR2F2','MEIS2','PAX6','ETV1'))
human_In_seurat <- human_RNA_seurat[,human_RNA_seurat$greenleaf_cluster %in% c('c1','c3')]
DefaultAssay(human_In_seurat) <- 'originalexp'
human_In_seurat <- NormalizeData(human_In_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
human_In_seurat <- FindVariableFeatures(human_In_seurat,selection.method = 'vst',nfeatures = 2000)
human_In_seurat <- ScaleData(human_In_seurat)
human_In_seurat <- RunPCA(human_In_seurat,assay = 'originalexp',features = VariableFeatures(human_In_seurat),npcs = 50)
ElbowPlot(human_In_seurat,ndims = 50)
human_In_seurat <- FindNeighbors(human_In_seurat,dims = 1:30)
human_In_seurat <- FindClusters(human_In_seurat,resolution = 0.4)
human_In_seurat <- RunUMAP(human_In_seurat,dims = 1:30,reduction.name = 'Dim',reduction.key = 'Dim_')

p1 <- DimPlot(human_In_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'Dim') + theme(aspect.ratio = 1)
p2 <- DimPlot(human_In_seurat,group.by = 'Age',label = FALSE,reduction = 'Dim') + theme(aspect.ratio = 1)

DimPlot(human_In_seurat,group.by = 'DF_classification',label = FALSE,reduction = 'Dim')
DimPlot(human_In_seurat,group.by = 'greenleaf_cluster',label = FALSE,reduction = 'Dim')

pdf(file = './res/step_2_fig_210702/human_In_seurat_all_marker_featureplot.pdf',width = 15,height = 10)
geneSetAveragePlot(human_In_seurat,genes = c('DLX2','DLX5','GAD2','LHX6','SST','SP8','NR2F2',"FOXP2",
                                             "LHFPL6","EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","PRKCB","TRPM3",
                                             "CDIN1","HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
                   object.class = 'seurat',assay = 'converted',embedding = 'Dim',reduction_key = 'Dim',plot.type = 'panels',
                   scaled = FALSE,color.palette = c("lightgrey", "blue"),aspectratio = 1,print = TRUE,num.panel.rows = 4)
dev.off()

#label transfer
macaque_In_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InCGE','InMGE')]
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
macaque_to_human_anno[which(macaque_to_human_anno[,1] == ''),1] <- macaque_to_human_anno[which(macaque_to_human_anno[,1] == ''),2]

temp <- macaque_In_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_In_seurat[['converted']] <- temp
rm(temp)
DefaultAssay(macaque_In_seurat) <- 'converted'
DefaultAssay(human_In_seurat) <- 'converted'
macaque_In_seurat <- NormalizeData(macaque_In_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
human_In_seurat <- NormalizeData(human_In_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
macaque_In_seurat <- FindVariableFeatures(macaque_In_seurat,selection.method = 'vst',nfeatures = 2000)
human_In_seurat <- FindVariableFeatures(human_In_seurat,selection.method = 'vst',nfeatures = 2000)
macaque_In_seurat <- ScaleData(macaque_In_seurat)
human_In_seurat <- ScaleData(human_In_seurat)

In.anchors <- FindTransferAnchors(reference = macaque_In_seurat,query = human_In_seurat,
                                  reference.assay = 'converted',query.assay = 'converted',dims = 1:30,verbose = TRUE)
predictions <- TransferData(anchorset = In.anchors,refdata = macaque_In_seurat$cell_type,dims = 1:30)
human_In_seurat <- AddMetaData(human_In_seurat,metadata = predictions)

p1 <- DimPlot(human_In_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE,reduction = 'Dim') + theme(aspect.ratio = 1)
p2 <- DimPlot(human_In_seurat,group.by = 'Age',label = FALSE,reduction = 'Dim') + theme(aspect.ratio = 1)
p3 <- DimPlot(human_In_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE,reduction = 'Dim') + theme(aspect.ratio = 1)
p4 <- DimPlot(human_In_seurat,group.by = 'DF_classification',label = FALSE,reduction = 'Dim') + theme(aspect.ratio = 1)
p5 <- DimPlot(human_In_seurat,group.by = 'greenleaf_cluster',label = FALSE,reduction = 'Dim') + theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_210702/human_In_seurat_dimplot.pdf',width = 15,height = 10)
p1+p2+p3+p4+p5+plot_layout(nrow = 2)
dev.off()


########################
##compare with PDhuman##
########################

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')

#preprocess of PDhuman_RNA_seurat
temp <- PDhuman_RNA_seurat@assays$RNA@counts
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
PDhuman_RNA_seurat[['converted']] <- temp

DefaultAssay(PDhuman_RNA_seurat) <- 'RNA'
PDhuman_RNA_seurat <- NormalizeData(PDhuman_RNA_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
PDhuman_RNA_seurat <- SCTransform(object = PDhuman_RNA_seurat,assay = 'RNA',new.assay.name = 'SCT',
                                  variable.features.n = 3000,vars.to.regress = NULL,verbose = TRUE)
DefaultAssay(PDhuman_RNA_seurat) <- 'SCT'

PDhuman_RNA_seurat <- RunPCA(PDhuman_RNA_seurat,assay = 'SCT',npcs = 50)
ElbowPlot(PDhuman_RNA_seurat,ndims = 50)

PDhuman_RNA_seurat <- FindNeighbors(PDhuman_RNA_seurat,dims = 1:40)
PDhuman_RNA_seurat <- FindClusters(PDhuman_RNA_seurat,resolution = 0.8)
PDhuman_RNA_seurat <- RunUMAP(PDhuman_RNA_seurat,dims = 1:40)

p1 <- DimPlot(PDhuman_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(PDhuman_RNA_seurat,group.by = 'Cluster',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_210702/PDhuman_RNA_seurat_dimplot.pdf',width = 12,height = 6)
p1+p2
dev.off()

pdf(file = './res/step_2_fig_210702/PDhuman_RNA_seurat_all_marker_featureplot.pdf',width = 15,height = 10)
geneSetAveragePlot(PDhuman_RNA_seurat,genes = c('DLX2','DLX5','GAD2','LHX6','SST','SP8','NR2F2',"FOXP2",
                                                "EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","PRKCB","TRPM3",
                                                "HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
                   object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c("lightgrey", "blue"),aspectratio = 1,print = TRUE,num.panel.rows = 4)
dev.off()


VlnPlot(object = PDhuman_RNA_seurat,features = c('DLX2','DLX5','GAD2','LHX6','SST','SP8','NR2F2',"FOXP2",
                                                 "EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","PRKCB","TRPM3",
                                                 "HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
        assay = 'RNA',slot = 'data',group.by = 'Cluster')

#subset PDhuman_In_seurat
PDhuman_RNA_seurat <- readRDS(file = './data/public/A_Single_Cell_Transcriptomic_Atlas_of_Human_Neocortical_Development_during_Mid_gestation/PD_human_RNA_seurat_210312.rds')
PDhuman_In_seurat <- PDhuman_RNA_seurat[,PDhuman_RNA_seurat$Cluster %in% c('InCGE','InMGE')]
DefaultAssay(PDhuman_In_seurat) <- 'RNA'
PDhuman_In_seurat <- NormalizeData(PDhuman_In_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
PDhuman_In_seurat <- FindVariableFeatures(PDhuman_In_seurat,selection.method = 'vst',nfeatures = 2000)
PDhuman_In_seurat <- ScaleData(PDhuman_In_seurat)

PDhuman_In_seurat <- RunPCA(PDhuman_In_seurat,assay = 'RNA',npcs = 50)
ElbowPlot(PDhuman_In_seurat,ndims = 50)

PDhuman_In_seurat <- FindNeighbors(PDhuman_In_seurat,dims = 1:20)
PDhuman_In_seurat <- FindClusters(PDhuman_In_seurat,resolution = 0.4)
PDhuman_In_seurat <- RunUMAP(PDhuman_In_seurat,dims = 1:20)

p1 <- DimPlot(PDhuman_In_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(PDhuman_In_seurat,group.by = 'Cluster',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(PDhuman_In_seurat,group.by = 'Donor',label = FALSE) + theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_210702/PDhuman_In_seurat_all_marker_featureplot.pdf',width = 15,height = 10)
geneSetAveragePlot(genes = c('DLX2','DLX5','GAD2','LHX6','SST','SP8','NR2F2',"FOXP2",
                             "EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","PRKCB","TRPM3",
                             "HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
                   object = PDhuman_In_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c("lightgrey", "blue"),aspectratio = 1,num.panel.rows = 4)
dev.off()

#transfer label
temp <- PDhuman_In_seurat@assays$RNA@counts
human_to_human_anno <- read.csv(file = './data/reference/GRCh38_to_GRCh38.csv')
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = human_to_human_anno)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
PDhuman_In_seurat[['converted']] <- temp

macaque_In_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InPSB','InCGE','InMGE')]
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
macaque_to_human_anno[which(macaque_to_human_anno[,1] == ''),1] <- macaque_to_human_anno[which(macaque_to_human_anno[,1] == ''),2]
temp <- macaque_In_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_to_human_anno)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_In_seurat[['converted']] <- temp
rm(temp)

DefaultAssay(macaque_In_seurat) <- 'converted'
DefaultAssay(PDhuman_In_seurat) <- 'converted'


macaque_In_seurat <- NormalizeData(macaque_In_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
PDhuman_In_seurat <- NormalizeData(PDhuman_In_seurat,normalization.method = 'LogNormalize',scale.factor = 10000)
macaque_In_seurat <- FindVariableFeatures(macaque_In_seurat,selection.method = 'vst',nfeatures = 2000)
PDhuman_In_seurat <- FindVariableFeatures(PDhuman_In_seurat,selection.method = 'vst',nfeatures = 2000)
macaque_In_seurat <- ScaleData(macaque_In_seurat)
PDhuman_In_seurat <- ScaleData(PDhuman_In_seurat)

In.anchors <- FindTransferAnchors(reference = macaque_In_seurat,query = PDhuman_In_seurat,
                                  reference.assay = 'converted',query.assay = 'converted',dims = 1:7,verbose = TRUE)
predictions <- TransferData(anchorset = In.anchors,refdata = macaque_In_seurat$cell_type,dims = 1:7)
PDhuman_In_seurat <- AddMetaData(PDhuman_In_seurat,metadata = predictions)

p1 <- DimPlot(PDhuman_In_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p2 <- DimPlot(PDhuman_In_seurat,group.by = 'Cluster',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)
p3 <- DimPlot(PDhuman_In_seurat,group.by = 'Donor',label = FALSE) + theme(aspect.ratio = 1)
p4 <- DimPlot(PDhuman_In_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1)

pdf(file = './res/step_2_fig_210702/PDhuman_In_dimplot.pdf',width = 15,height = 10)
p1+p2+p3+p4+plot_layout(nrow = 2)
dev.off()


############################
##compare with 110 PCD DFC##
############################

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
macaque_DFC_seurat <- readRDS(file = './data/public/Spatiotemporal_transcriptomic_divergence_across_human_and_macaque_brain_development/DFC_scRNAseq/macaque_DFC_seurat.rds')

#prepare data
macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

macaque_annotation <- data.frame(gene_name=macaque_annotation$gene_name,gene_id=macaque_annotation$gene_id,
                                 gene_name_1=macaque_annotation$gene_name,gene_id_1=macaque_annotation$gene_id)

temp <- macaque_DFC_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_annotation,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

macaque_DFC_seurat[['converted']] <- temp
DefaultAssay(macaque_DFC_seurat) <- 'RNA'

p1 <- DimPlot(macaque_DFC_seurat,group.by = 'subtype',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_DFC_seurat[,macaque_DFC_seurat$ctype == 'IntN'],group.by = 'subtype',label = TRUE,repel = TRUE)

p1+p2

geneSetAveragePlot(genes = c('GAD2','LHX6','SST','SP8','NR2F2',"FOXP2",
                             "EYA2","TSHZ1","FOXP4","PLXNC1","FRAS1","TRPM3",
                             "HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
                   object = macaque_DFC_seurat[,macaque_DFC_seurat$ctype == 'IntN'],object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c("lightgrey", "blue"),aspectratio = 1,num.panel.rows = 4)

#subset interneuron
macaque_DFC_seurat <- macaque_DFC_seurat[,macaque_DFC_seurat$ctype == 'IntN']
macaque_DFC_seurat <- my_process_seurat(object = macaque_DFC_seurat,preprocess = TRUE,dim_to_use = 10)
macaque_DFC_seurat <- my_process_seurat(object = macaque_DFC_seurat,preprocess = FALSE,dim_to_use = 10)

p1 <- DimPlot(macaque_DFC_seurat,group.by = 'subtype',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_DFC_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p3 <- DimPlot(macaque_DFC_seurat,group.by = 'donor',label = FALSE)

pdf(file = './res/step_2_fig_210702/macaque_DFC_seurat_dimplot.pdf',width = 15,height = 5)
p1+p2+p3
dev.off()

pdf(file = './res/step_2_fig_210702/macaque_DFC_seurat_IntN_featureplot.pdf',width = 10,height = 7.5)
geneSetAveragePlot(genes = c('GAD2','LHX6','SST','SP8','NR2F2',"FOXP2",
                             "EYA2","TSHZ1","FOXP4","PLXNC1","FRAS1","TRPM3",
                             "HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
                   object = macaque_DFC_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c("lightgrey", "blue"),aspectratio = 1,num.panel.rows = 4)
dev.off()

#label transfer
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE','InPSB')]
temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_annotation,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)

macaque_RNA_seurat[['converted']] <- temp

macaque_DFC_seurat <- my_process_seurat(object = macaque_DFC_seurat,preprocess = TRUE,assay = 'converted')
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,preprocess = TRUE,assay = 'converted')

In.anchors <- FindTransferAnchors(reference = macaque_RNA_seurat,query = macaque_DFC_seurat,
                                  reference.assay = 'converted',query.assay = 'converted',dims = 1:10,verbose = TRUE)
predictions <- TransferData(anchorset = In.anchors,refdata = macaque_RNA_seurat$cell_type,dims = 1:10)
macaque_DFC_seurat <- AddMetaData(macaque_DFC_seurat,metadata = predictions)

p1 <- DimPlot(macaque_DFC_seurat,group.by = 'subtype',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_DFC_seurat,group.by = 'predicted.id',label = TRUE,repel = TRUE)

pdf(file = './res/step_2_fig_210702/macaque_DFC_seurat_label_transfer.pdf',width = 10,height = 5)
p1+p2
dev.off()




#####################################
##compare with 110 PCD ,macaque STR##
#####################################

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/macaque_RNA_seurat_annotated_210629.rds')
macaque_STR_seurat <- readRDS(file = './data/public/Spatiotemporal_transcriptomic_divergence_across_human_and_macaque_brain_development/STR_scRNAseq/macaque_STR_seurat.rds')

#prepare converted assay
macaque_annotation <- rtracklayer::import(con = './data/reference/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_annotation <- rtracklayer::as.data.frame(macaque_annotation)
macaque_annotation <- macaque_annotation[,c('gene_id','gene_name')]
macaque_annotation <- macaque_annotation[!duplicated(macaque_annotation$gene_id),]
macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_name"] <- macaque_annotation[which(is.na(macaque_annotation$gene_name)),"gene_id"]

macaque_annotation <- data.frame(gene_name=macaque_annotation$gene_name,gene_id=macaque_annotation$gene_id,
                                 gene_name_1=macaque_annotation$gene_name,gene_id_1=macaque_annotation$gene_id)

temp <- macaque_RNA_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_annotation,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_RNA_seurat[['converted']] <- temp

temp <- macaque_STR_seurat@assays$RNA@counts
temp <- My_Convert_Homology_Gene_ID(express_matrix = temp,anno = macaque_annotation,filter_anno = TRUE)
temp <- CreateAssayObject(counts = temp,min.cells = 0,min.features = 0)
macaque_STR_seurat[['converted']] <- temp

#marker gene express
p1 <- DimPlot(macaque_STR_seurat,group.by = 'ctype',label = TRUE,repel = TRUE)
p2 <- DimPlot(macaque_STR_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)

pdf(file = './res/step_2_fig_210702/macaque_STR_seurat_dimplot.pdf',width = 12,height = 5)
p1+p2
dev.off()

pdf(file = './res/step_2_fig_210702/macaque_STR_seurat_featureplot.pdf',width = 15,height = 10)
geneSetAveragePlot(genes = c('GAD2','LHX6','SST','SP8','NR2F2',"FOXP2",
                             "EYA2","TSHZ1","FOXP4","ALK","PLXNC1","FRAS1","TRPM3",
                             "HS6ST2","ESRRG","PBX3","MEIS2","PTPRM","ROR2"),
                   object = macaque_STR_seurat,object.class = 'seurat',assay = 'RNA',embedding = 'umap',reduction_key = 'UMAP',
                   plot.type = 'panels',scaled = FALSE,color.palette = c("lightgrey", "blue"),aspectratio = 1,num.panel.rows = 4)
dev.off()

pdf(file = './res/step_2_fig_210702/macaque_STR_seurat_vlnplot.pdf',width = 18,height = 8)
VlnPlot(macaque_STR_seurat,features = c('FOXP2','EYA2','TSHZ1','FOXP4','PBX3','MEIS2'),
        assay = 'RNA',group.by = 'subtype')
dev.off()

#integration
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE','InPSB')]
macaque_STR_seurat <- macaque_STR_seurat[,macaque_STR_seurat$ctype %in% c('NascIntN','MSN','NascMSN')]

macaque_GABA_liger <- createLiger(list(cortex=macaque_RNA_seurat@assays$converted@counts,
                                       STR=macaque_STR_seurat@assays$converted@counts))

macaque_GABA_liger <- normalize(macaque_GABA_liger)
macaque_GABA_liger <- selectGenes(macaque_GABA_liger)
macaque_GABA_liger <- scaleNotCenter(macaque_GABA_liger)

macaque_GABA_liger <- optimizeALS(macaque_GABA_liger,k = 20,lambda = 0.5)
macaque_GABA_liger <- quantile_norm(macaque_GABA_liger)
macaque_GABA_liger <- louvainCluster(macaque_GABA_liger,resolution = 0.4)
macaque_GABA_liger <- runUMAP(macaque_GABA_liger)
all.plots <- plotByDatasetAndCluster(macaque_GABA_liger, axis.labels = c('UMAP 1', 'UMAP 2'), return.plots = T)
all.plots[[1]]+all.plots[[2]]

macaque_GABA_seurat <- ligerToSeurat(macaque_GABA_liger)
macaque_GABA_seurat@meta.data[1:3,]
macaque_GABA_seurat$cell_type <- 'NA'

macaque_GABA_seurat@meta.data[paste('cortex',colnames(macaque_RNA_seurat),sep = '_'),"cell_type"] <- as.character(macaque_RNA_seurat$cell_type)
macaque_GABA_seurat@meta.data[paste('STR',colnames(macaque_STR_seurat),sep = '_'),"cell_type"] <- as.character(macaque_STR_seurat$subtype)

p1 <- DimPlot(macaque_GABA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)

macaque_GABA_seurat$InPSB <- c(macaque_GABA_seurat$cell_type == 'InPSB')
macaque_GABA_seurat$NascMSN2 <- c(macaque_GABA_seurat$cell_type == 'NascMSN2')

p2 <- DimPlot(macaque_GABA_seurat,group.by = 'InPSB',label = FALSE,repel = TRUE,cols = c('lightgrey','black'))
p3 <- DimPlot(macaque_GABA_seurat,group.by = 'NascMSN2',label = FALSE,repel = TRUE,cols = c('lightgrey','black'))

pdf(file = './res/step_2_fig_210702/macaque_STR_seurat_integration.pdf',width = 5,height = 15)
p1+p2+p3+plot_layout(nrow = 1)
dev.off()