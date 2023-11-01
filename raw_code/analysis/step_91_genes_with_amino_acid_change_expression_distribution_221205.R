#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: genes with amino acid change expression distribution            ##
## Data: 2022.12.05                                                                ##
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

# load data ---------------------------------------------------------------
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

gene_list_backup <- read.csv(file = './res/step_91_fig_221205/genes_with_amino_acid_change_in_human.csv')
gene_list <- gene_list_backup$Gene.name
gene_list <- unique(gene_list)

#refine gene_list
temp_gene_list <- gene_list[gene_list %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)]
temp <- gene_list[!(gene_list %in% rownames(Greenleaf_RNA_Seurat@assays$RNA@counts))]
gene_list <- c(temp_gene_list,c('COQ8A','ELOA','ITPRID2','CFAP65','CFAP44','LRRC2','ACY1','ABHD18','ADGRF1','ELAPOR2',
                                'CFAP157','PRXL2C','SLF2','CFAP46','COMMD3','ZNF22-AS1','MFRP','MIF','PLAAT5','CRACR2B',
                                'LINC01559','CRACR2A','KNL1','CEMIP','MEAK7','MYOCD','TSPOAP1','ACE','MRPL58','CEP131',
                                'TEPSIN','FAAP100','PPAN','PTBP1','EPPIN','MMP11','MRTFA','RTL9','CFAP47'))
gene_list <- unique(gene_list)

saveRDS(object = gene_list,file = './res/step_91_fig_221205/human_gene_list.rds')

# find homology gene in mouse and macaque ---------------------------------
#load annotation
human_to_macaque_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCh38_to_Mmul10.csv')
human_to_mouse_anno <- read.csv(file = './data/reference/BioMart_release_105/GRCh38_to_GRCm39.csv')

#find homology gene
human_to_macaque_anno <- human_to_macaque_anno[human_to_macaque_anno[,2] %in% gene_list,]
human_to_macaque_anno <- human_to_macaque_anno[human_to_macaque_anno[,4] != '',]
human_to_macaque_anno <- human_to_macaque_anno[,c(2,4)]
colnames(human_to_macaque_anno) <- c('human','macaque')
rownames(human_to_macaque_anno) <- human_to_macaque_anno$human

human_to_mouse_anno <- human_to_mouse_anno[human_to_mouse_anno[,2] %in% gene_list,]
human_to_mouse_anno <- human_to_mouse_anno[human_to_mouse_anno[,4] != '',]
duplicated_gene_list <- unique(human_to_mouse_anno[duplicated(human_to_mouse_anno[,2]),2])
temp <- unlist(base::lapply(X = duplicated_gene_list,FUN = function(x){
  temp <- human_to_mouse_anno[human_to_mouse_anno[,2] == x,4]
  for (j in temp) {
    if(j %in% rownames(mouse_multiome_Seurat@assays$RNA@counts)){
      return(j)
    }
  }
  return(NA)
}))
i <- human_to_mouse_anno[!(human_to_mouse_anno[,2] %in% duplicated_gene_list),]
i <- i[,c(2,4)]
colnames(i) <- c('human','mouse')
human_to_mouse_anno <- data.frame(human = duplicated_gene_list,mouse = temp)
human_to_mouse_anno <- rbind(i,human_to_mouse_anno)
rownames(human_to_mouse_anno) <- human_to_mouse_anno$human

#macaque homology gene
temp <- gene_list[gene_list %in% rownames(human_to_macaque_anno)]
homology_gene_list <- human_to_macaque_anno[temp,"macaque"]
names(homology_gene_list) <- temp

temp <- gene_list[!(gene_list %in% rownames(human_to_macaque_anno))]
temp <- temp[temp %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)]
names(temp) <- temp

macaque_homology_gene_list <- c(homology_gene_list,temp)
saveRDS(macaque_homology_gene_list,file = './res/step_91_fig_221205/macaque_homology_gene_list.rds')

#mouse homology gene
temp <- gene_list[gene_list %in% rownames(human_to_mouse_anno)]
homology_gene_list <- human_to_mouse_anno[temp,"mouse"]
names(homology_gene_list) <- temp

temp <- gene_list[!(gene_list %in% rownames(human_to_mouse_anno))]
i <- temp
temp <- str_to_title(temp)
names(temp) <- i
temp <- temp[temp %in% rownames(mouse_multiome_Seurat@assays$RNA@counts)]

mouse_homology_gene_list <- c(homology_gene_list,temp)
saveRDS(object = mouse_homology_gene_list,file = './res/step_91_fig_221205/mouse_homology_gene_list.rds')

# human gene plot ---------------------------------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,!(Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('End','Per'))]
gene_list <- readRDS(file = './res/step_91_fig_221205/human_gene_list.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#filter gene with no express
temp <- Greenleaf_RNA_Seurat@assays$RNA@counts[gene_list,]
temp <- temp[rowSums(temp) > 10,]
gene_list <- rownames(temp)

plot_matrix <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = gene_list,cols = c("lightgrey","blue"),
                          scale = FALSE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)

gene_list <- unlist(base::lapply(X = gene_list,FUN = function(x){
  temp <- plot_matrix[plot_matrix$features.plot == x,]
  if(max(temp$pct.exp) > 20){
    return(x)
  }else{
    return(NA)
  }
}))
gene_list <- gene_list[!is.na(gene_list)]

#dotplot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_RNA_Seurat$ReAnno_celltype]
Greenleaf_RNA_Seurat$ReAnno_celltype <- factor(x = Greenleaf_RNA_Seurat$ReAnno_celltype,levels = cell_type_list)
i <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = gene_list,
             group.by = 'ReAnno_celltype',scale = TRUE) + 
  coord_flip() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))


pdf(file = '/home/sunym/temp/human_dotplot.pdf',width = 10,height = 0.5*length(gene_list))
print(i)
dev.off()

#featureplot
FeaturePlot(object = Greenleaf_RNA_Seurat,features = 'VCAM1',pt.size = 0.1,order = FALSE,slot = 'data') + theme(aspect.ratio = 1)

#class genes
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

#plot by class
i <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = Progenitor_gene,
             group.by = 'ReAnno_celltype',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/human_Progenitor_gene_dotplot.pdf',width = 6,height = 0.3*length(Progenitor_gene))
print(i)
dev.off()

i <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = Ex_gene,
             group.by = 'ReAnno_celltype',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/human_Ex_gene_dotplot.pdf',width = 6,height = 0.3*length(Ex_gene))
print(i)
dev.off()

i <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = Glial_gene,
             group.by = 'ReAnno_celltype',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/human_Glial_gene_dotplot.pdf',width = 6,height = 0.3*length(Glial_gene))
print(i)
dev.off()

i <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = Global_gene,
             group.by = 'ReAnno_celltype',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/human_Global_gene_dotplot.pdf',width = 6,height = 0.3*length(Global_gene))
print(i)
dev.off()

i <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = In_gene,
             group.by = 'ReAnno_celltype',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/human_In_gene_dotplot.pdf',width = 6,height = 6)
print(i)
dev.off()

# macaque gene express ----------------------------------------------------

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
gene_list <- readRDS(file = './res/step_91_fig_221205/macaque_homology_gene_list.rds')
gene_list <- gene_list[gene_list %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)]

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_RNA_Seurat$ReAnno_celltype & cell_type_list %in% macaque_multiome_Seurat$cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

## correlation with human --------------------------------------------------
#use all genes
human_express_matrix <- AverageExpression(object = Greenleaf_RNA_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'ReAnno_celltype',slot = 'data',verbose = TRUE)
macaque_express_matrix <- AverageExpression(object = macaque_multiome_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'data',verbose = TRUE)

cor_list <- unlist(base::lapply(X = names(gene_list),FUN = function(x){
  temp_human <- human_express_matrix$RNA[x,cell_type_list]
  temp_macaque <- macaque_express_matrix$RNA[gene_list[x],cell_type_list]
  return(cor(x = temp_human,y = temp_macaque,method = 'pearson'))
}))
names(cor_list) <- names(gene_list)
cor_list <- cor_list[!(is.na(cor_list))]

hist(cor_list)
all_cor_list <- cor_list

#only use genes expressed in human
gene_list_backup <- gene_list
gene_list <- names(gene_list)
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,!(Greenleaf_RNA_Seurat$ReAnno_celltype %in% c('End','Per'))]
temp <- Greenleaf_RNA_Seurat@assays$RNA@counts[gene_list,]
temp <- temp[rowSums(temp) > 10,]
gene_list <- rownames(temp)

plot_matrix <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = gene_list,cols = c("lightgrey","blue"),
                          scale = FALSE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)
gene_list <- unlist(base::lapply(X = gene_list,FUN = function(x){
  temp <- plot_matrix[plot_matrix$features.plot == x,]
  if(max(temp$pct.exp) > 20){
    return(x)
  }else{
    return(NA)
  }
}))
gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list_backup[gene_list]

cor_list <- unlist(base::lapply(X = names(gene_list),FUN = function(x){
  temp_human <- human_express_matrix$RNA[x,cell_type_list]
  temp_macaque <- macaque_express_matrix$RNA[gene_list[x],cell_type_list]
  return(cor(x = temp_human,y = temp_macaque,method = 'pearson'))
}))
names(cor_list) <- names(gene_list)
cor_list <- cor_list[!(is.na(cor_list))]

hist(cor_list)
expressed_cor_list <- cor_list

#correlation distribution
temp <- data.frame(cor = c(all_cor_list,expressed_cor_list),
                   group = c(rep('all',times = length(all_cor_list)),
                             rep('expressed',times = length(expressed_cor_list))))

pdf(file = './res/step_91_fig_221205/human_macaque_gene_expression_correlation_density_plot.pdf',width = 4,height = 4)
ggplot(data = temp,aes(x = cor,color = group)) + 
  geom_density() + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 12,face = 'bold',hjust = 0.5)) + 
  labs(title = 'expression corretion')
dev.off()

gene_list <- names(all_cor_list)
temp <- rowSums(Greenleaf_RNA_Seurat@assays$RNA@counts[gene_list,])
temp <- log1p(temp)
temp <- data.frame(expression = temp,cor = all_cor_list)

ggplot(data = temp,aes(x = expression,y = cor)) + 
  geom_point() + 
  theme_bw() + 
  theme(aspect.ratio = 1)

#well, not high correlation between expression and correlation

## human identified gene expression ----------------------------------------
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

#modify annotation
gene_list <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene,Global_gene)
human_to_macaque_anno <- readRDS(file = './res/step_91_fig_221205/macaque_homology_gene_list.rds')
temp <- gene_list[!(gene_list %in% names(human_to_macaque_anno))]
names(temp) <- temp
human_to_macaque_anno <- c(human_to_macaque_anno,temp)

#add gene to macaque expression matrix
gene_list <- human_to_macaque_anno[gene_list]
gene_list <- gene_list[!(gene_list %in% rownames(macaque_multiome_Seurat))]
temp <- matrix(data = 0,nrow = length(gene_list),ncol = ncol(macaque_multiome_Seurat@assays$RNA@counts))
rownames(temp) <- gene_list
colnames(temp) <- colnames(macaque_multiome_Seurat@assays$RNA@counts)
temp <- as(temp,'dgCMatrix')

macaque_express_matrix <- rbind(macaque_multiome_Seurat@assays$RNA@counts,temp)
macaque_express_matrix <- CreateSeuratObject(counts = macaque_express_matrix,project = 'temp',assay = 'RNA',meta.data = macaque_multiome_Seurat@meta.data,min.cells = 0,min.features = 0)
macaque_express_matrix <- NormalizeData(object = macaque_express_matrix,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)

macaque_express_matrix <- macaque_express_matrix[,macaque_express_matrix$cell_type %in% cell_type_list]
macaque_express_matrix$cell_type <- factor(macaque_express_matrix$cell_type,levels = cell_type_list)

#progenitor gene dot plot
Progenitor_gene <- human_to_macaque_anno[Progenitor_gene]
names(Progenitor_gene) <- NULL
i <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = Progenitor_gene,
             group.by = 'cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/macaque_Progenitor_gene_dotplot.pdf',width = 6,height = 0.3*length(Progenitor_gene))
print(i)
dev.off()

#Ex gene dot plot
Ex_gene <- human_to_macaque_anno[Ex_gene]
names(Ex_gene) <- NULL
i <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = Ex_gene,
             group.by = 'cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/macaque_Ex_gene_dotplot.pdf',width = 6,height = 0.3*length(Ex_gene))
print(i)
dev.off()

#In gene dot plot
In_gene <- human_to_macaque_anno[In_gene]
names(In_gene) <- NULL
i <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = In_gene,
             group.by = 'cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/macaque_In_gene_dotplot.pdf',width = 6,height = 6)
print(i)
dev.off()

#Glial gene dot plot
Glial_gene <- human_to_macaque_anno[Glial_gene]
names(Glial_gene) <- NULL
i <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = Glial_gene,
             group.by = 'cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/macaque_Glial_gene_dotplot.pdf',width = 6,height = 0.3*length(Glial_gene))
print(i)
dev.off()

#Global gene dot plot
Global_gene <- human_to_macaque_anno[Global_gene]
names(Global_gene) <- NULL
i <- DotPlot(object = macaque_express_matrix,assay = 'RNA',features = Global_gene,
             group.by = 'cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/macaque_Global_gene_dotplot.pdf',width = 6,height = 0.3*length(Global_gene))
print(i)
dev.off()

# mouse gene express ------------------------------------------------------

## correlation with human gene expression ----------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

human_to_macaque_anno <- readRDS(file = './res/step_91_fig_221205/macaque_homology_gene_list.rds')
human_to_mouse_anno <- readRDS(file = './res/step_91_fig_221205/mouse_homology_gene_list.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_RNA_Seurat$ReAnno_celltype & cell_type_list %in% macaque_multiome_Seurat$cell_type & cell_type_list %in% mouse_multiome_Seurat$macaque_cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

#filter cell type
Greenleaf_RNA_Seurat <- Greenleaf_RNA_Seurat[,Greenleaf_RNA_Seurat$ReAnno_celltype %in% cell_type_list]
macaque_multiome_Seurat <- macaque_multiome_Seurat[,macaque_multiome_Seurat$cell_type %in% cell_type_list]
mouse_multiome_Seurat <- mouse_multiome_Seurat[,mouse_multiome_Seurat$macaque_cell_type %in% cell_type_list]

#get human expressed gene list
gene_list <- readRDS(file = './res/step_91_fig_221205/human_gene_list.rds')
temp <- Greenleaf_RNA_Seurat@assays$RNA@counts[gene_list,]
temp <- temp[rowSums(temp) > 10,]
gene_list <- rownames(temp)

plot_matrix <- my_dotplot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = gene_list,cols = c("lightgrey","blue"),
                          scale = FALSE,group.by = 'ReAnno_celltype',return_data_plot = TRUE)

gene_list <- unlist(base::lapply(X = gene_list,FUN = function(x){
  temp <- plot_matrix[plot_matrix$features.plot == x,]
  if(max(temp$pct.exp) > 20){
    return(x)
  }else{
    return(NA)
  }
}))
gene_list <- gene_list[!is.na(gene_list)]

#get averaged gene expression
human_express_matrix <- AverageExpression(object = Greenleaf_RNA_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'ReAnno_celltype',slot = 'data',verbose = TRUE)
macaque_express_matrix <- AverageExpression(object = macaque_multiome_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'data',verbose = TRUE)
mouse_express_matrix <- AverageExpression(object = mouse_multiome_Seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'macaque_cell_type',slot = 'data',verbose = TRUE)

#macaque human expression cor
macaque_gene_list <- gene_list[gene_list %in% names(human_to_macaque_anno)]
temp <- human_to_macaque_anno[macaque_gene_list]
temp <- temp[temp %in% rownames(macaque_multiome_Seurat@assays$RNA@counts)]
macaque_gene_list <- names(temp)

temp_human <- human_express_matrix$RNA[macaque_gene_list,cell_type_list]
macaque_gene_list <- human_to_macaque_anno[macaque_gene_list]
temp_macaque <- macaque_express_matrix$RNA[macaque_gene_list,cell_type_list]

macaque_cor <- unlist(base::lapply(X = 1:length(macaque_gene_list),FUN = function(x){
  temp <- cor(x = temp_human[x,],y = temp_macaque[x,],method = 'pearson')
  return(temp)
}))

#mouse human expression cor
mouse_gene_list <- gene_list[gene_list %in% names(human_to_mouse_anno)]
temp <- human_to_mouse_anno[mouse_gene_list]
temp <- temp[temp %in% rownames(mouse_multiome_Seurat)]
mouse_gene_list <- names(temp)

temp_human <- human_express_matrix$RNA[mouse_gene_list,cell_type_list]
mouse_gene_list <- human_to_mouse_anno[mouse_gene_list]
temp_mouse <- mouse_express_matrix$RNA[mouse_gene_list,cell_type_list]

mouse_cor <- unlist(base::lapply(X = 1:length(mouse_gene_list),FUN = function(x){
  temp <- cor(x = temp_human[x,],y = temp_mouse[x,],method = 'pearson')
  return(temp)
}))

#plot correlation distribution
temp <- data.frame(cor = c(macaque_cor,mouse_cor),
                   group = c(rep('macaque',times = length(macaque_cor)),
                             rep('mouse',times = length(mouse_cor))))

pdf(file = './res/step_91_fig_221205/human_macaque_cor_human_mouse_cor_density_plot.pdf',width = 4,height = 4)
ggplot(data = temp,aes(x = cor,color = group)) + 
  geom_density() + 
  scale_color_manual(values = color_param$species[c('macaque','mouse')]) + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 12,face = 'bold',hjust = 0.5)) + 
  labs(title = 'expression corretion')
dev.off()

## human identified gene expression ----------------------------------------

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

human_to_mouse_anno <- readRDS(file = './res/step_91_fig_221205/mouse_homology_gene_list.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_RNA_Seurat$ReAnno_celltype & cell_type_list %in% mouse_multiome_Seurat$macaque_cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

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

#modify annotation
gene_list <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene,Global_gene)
temp <- gene_list[!(gene_list %in% names(human_to_mouse_anno))]
i <- temp
temp <- str_to_title(temp)
names(temp) <- i
human_to_mouse_anno <- c(human_to_mouse_anno,temp)

#add gene to mouse expression matrix
gene_list <- human_to_mouse_anno[gene_list]
gene_list <- gene_list[!(gene_list %in% rownames(mouse_multiome_Seurat))]
temp <- matrix(data = 0,nrow = length(gene_list),ncol = ncol(mouse_multiome_Seurat@assays$RNA@counts))
rownames(temp) <- gene_list
colnames(temp) <- colnames(mouse_multiome_Seurat@assays$RNA@counts)
temp <- as(temp,'dgCMatrix')

mouse_express_matrix <- rbind(mouse_multiome_Seurat@assays$RNA@counts,temp)
mouse_express_matrix <- CreateSeuratObject(counts = mouse_express_matrix,project = 'temp',assay = 'RNA',meta.data = mouse_multiome_Seurat@meta.data,min.cells = 0,min.features = 0)
mouse_express_matrix <- NormalizeData(object = mouse_express_matrix,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)

mouse_express_matrix <- mouse_express_matrix[,mouse_express_matrix$macaque_cell_type %in% cell_type_list]
mouse_express_matrix$macaque_cell_type <- factor(mouse_express_matrix$macaque_cell_type,levels = cell_type_list)

#progenitor gene dot plot
Progenitor_gene <- human_to_mouse_anno[Progenitor_gene]
names(Progenitor_gene) <- NULL
i <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = Progenitor_gene,
             group.by = 'macaque_cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/mouse_Progenitor_gene_dotplot.pdf',width = 6,height = 0.3*length(Progenitor_gene))
print(i)
dev.off()

#Ex gene dot plot
Ex_gene <- human_to_mouse_anno[Ex_gene]
names(Ex_gene) <- NULL
i <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = Ex_gene,
             group.by = 'macaque_cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/mouse_Ex_gene_dotplot.pdf',width = 6,height = 0.3*length(Ex_gene))
print(i)
dev.off()

#In gene dot plot
In_gene <- human_to_mouse_anno[In_gene]
names(In_gene) <- NULL
i <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = In_gene,
             group.by = 'macaque_cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/mouse_In_gene_dotplot.pdf',width = 6,height = 6)
print(i)
dev.off()

#Glial gene dot plot
Glial_gene <- human_to_mouse_anno[Glial_gene]
names(Glial_gene) <- NULL
i <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = Glial_gene,
             group.by = 'macaque_cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/mouse_Glial_gene_dotplot.pdf',width = 6,height = 0.3*length(Glial_gene))
print(i)
dev.off()

#Global gene dot plot
Global_gene <- human_to_mouse_anno[Global_gene]
names(Global_gene) <- NULL
i <- DotPlot(object = mouse_express_matrix,assay = 'RNA',features = Global_gene,
             group.by = 'macaque_cell_type',scale = TRUE) + 
  coord_flip() + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1))

pdf(file = './res/step_91_fig_221205/mouse_Global_gene_dotplot.pdf',width = 6,height = 0.3*length(Global_gene))
print(i)
dev.off()

# generate species combined dotplot ---------------------------------------

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

human_to_macaque_anno <- readRDS(file = './res/step_91_fig_221205/macaque_homology_gene_list.rds')
human_to_mouse_anno <- readRDS(file = './res/step_91_fig_221205/mouse_homology_gene_list.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

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
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)
col.min <- min(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#human plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 3) + 
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

pdf(file = './res/step_91_fig_221205/species_combined_dotplot_Progenitor_gene.pdf',width = 10,height = 8)
temp_human+temp_macaque+temp_mouse+plot_layout(ncol = 3)
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
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)
col.min <- min(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#human plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 2) + 
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

pdf(file = './res/step_91_fig_221205/species_combined_dotplot_Ex_gene.pdf',width = 10,height = 8)
temp_human+temp_macaque+temp_mouse+plot_layout(ncol = 3)
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
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)
col.min <- min(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#human plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 0.25) + 
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

pdf(file = './res/step_91_fig_221205/species_combined_dotplot_In_gene.pdf',width = 10,height = 8)
temp_human+temp_macaque+temp_mouse+plot_layout(ncol = 3)
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
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)
col.min <- min(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#human plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 2) + 
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

pdf(file = './res/step_91_fig_221205/species_combined_dotplot_Glial_gene.pdf',width = 10,height = 8)
temp_human+temp_macaque+temp_mouse+plot_layout(ncol = 3)
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
temp_macaque <- my_dotplot(object = macaque_express_matrix,assay = 'RNA',features = mk_gene_macaque,scale = TRUE,group.by = 'cell_type',return_data_plot = TRUE)
temp_mouse <- my_dotplot(object = mouse_express_matrix,assay = 'RNA',features = mk_gene_mouse,scale = TRUE,group.by = 'macaque_cell_type',return_data_plot = TRUE)

col.max <- max(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)
col.min <- min(temp_human$avg.exp.scaled,temp_macaque$avg.exp.scaled,temp_mouse$avg.exp.scaled)

scale.max <- max(temp_human$pct.exp,temp_macaque$pct.exp,temp_mouse$pct.exp)
scale.min <- 0

#human plot
temp_human <- DotPlot(object = Greenleaf_RNA_Seurat,assay = 'RNA',features = mk_gene,cols = c('lightgrey',color_param$species['human']),
                      col.min = col.min,col.max = col.max,group.by = 'ReAnno_celltype',scale = TRUE,scale.min = scale.min,scale.max = scale.max) + 
  theme_bw() + coord_flip() + theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + NoLegend() + 
  theme(aspect.ratio = 3) + 
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

pdf(file = './res/step_91_fig_221205/species_combined_dotplot_Global_gene.pdf',width = 10,height = 8)
temp_human+temp_macaque+temp_mouse+plot_layout(ncol = 3)
dev.off()


# generate species combined heatmap plot --------------------------------------

#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
mouse_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/mouse_multiome_Seurat_221009.rds')

human_to_macaque_anno <- readRDS(file = './res/step_91_fig_221205/macaque_homology_gene_list.rds')
human_to_mouse_anno <- readRDS(file = './res/step_91_fig_221205/mouse_homology_gene_list.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

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

## human plot --------------------------------------------------------------

#set param
species_list <- 'human'
express_matrix <- Greenleaf_RNA_Seurat
group.by <- 'ReAnno_celltype'

#get cell type list and mk gene
cell_type_list <- levels(express_matrix@meta.data[,group.by])
mk_gene <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene)
names(mk_gene) <- c(rep('Progenitor',times = length(Progenitor_gene)),rep('Ex',times = length(Ex_gene)),
                    rep('In',times = length(In_gene)),rep('Glial',times = length(Glial_gene)))

#aggregate express matrix
express_matrix <- AverageExpression(object = express_matrix,assays = 'RNA',return.seurat = FALSE,group.by = group.by,slot = 'data',verbose = TRUE)
express_matrix <- express_matrix$RNA[mk_gene,]
express_matrix <- t(scale(t(express_matrix)))
express_matrix <- express_matrix[mk_gene,cell_type_list]

#add annotation
top_anno <- HeatmapAnnotation(cell_type = colnames(express_matrix),
                              col = list(cell_type = color_param$celltype[cell_type_list]),
                              show_legend = TRUE,which = 'column',border = FALSE,show_annotation_name = FALSE)

col_fun <- colorRamp2(breaks = c(-2.5,0,2.5),colors = c('lightgrey','white',color_param$species[species_list]))

#plot
h_human <- Heatmap(matrix = express_matrix,show_column_names = FALSE,show_row_names = FALSE,cluster_columns = FALSE,
                   cluster_rows = FALSE,row_split = factor(names(mk_gene),levels = c('Progenitor','Ex','In','Glial')),
                   border = TRUE,top_annotation = top_anno,col = col_fun,name = 'human',
                   width = unit(0.125*length(cell_type_list),'inches'),height = unit(0.125*length(mk_gene),'inches'))

## macaque plot --------------------------------------------------------------

#set param
species_list <- 'macaque'
express_matrix <- macaque_express_matrix
group.by <- 'cell_type'

#get cell type list and mk gene
cell_type_list <- levels(express_matrix@meta.data[,group.by])
mk_gene <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene)
mk_gene <- human_to_macaque_anno[mk_gene]
names(mk_gene) <- c(rep('Progenitor',times = length(Progenitor_gene)),rep('Ex',times = length(Ex_gene)),
                    rep('In',times = length(In_gene)),rep('Glial',times = length(Glial_gene)))

#aggregate express matrix
express_matrix <- AverageExpression(object = express_matrix,assays = 'RNA',return.seurat = FALSE,group.by = group.by,slot = 'data',verbose = TRUE)
express_matrix <- express_matrix$RNA[mk_gene,]
express_matrix <- t(scale(t(express_matrix)))
express_matrix <- express_matrix[mk_gene,cell_type_list]

#add annotation
top_anno <- HeatmapAnnotation(cell_type = colnames(express_matrix),
                              col = list(cell_type = color_param$celltype[cell_type_list]),
                              show_legend = FALSE,which = 'column',border = FALSE,show_annotation_name = FALSE)

col_fun <- colorRamp2(breaks = c(-2.5,0,2.5),colors = c('lightgrey','white',color_param$species[species_list]))

#plot
h_macaque <- Heatmap(matrix = express_matrix,show_column_names = FALSE,show_row_names = FALSE,cluster_columns = FALSE,
                     cluster_rows = FALSE,row_split = factor(names(mk_gene),levels = c('Progenitor','Ex','In','Glial')),
                     border = TRUE,top_annotation = top_anno,col = col_fun,name = 'macaque',row_order = h_human@row_order,
                     width = unit(0.125*length(cell_type_list),'inches'),height = unit(0.125*length(mk_gene),'inches'))

## mouse plot --------------------------------------------------------------

#set param
species_list <- 'mouse'
express_matrix <- mouse_express_matrix
group.by <- 'macaque_cell_type'

#get cell type list and mk gene
cell_type_list <- levels(express_matrix@meta.data[,group.by])
cell_type_list <- cell_type_list[cell_type_list %in% express_matrix@meta.data[,group.by]]
mk_gene <- c(Progenitor_gene,Ex_gene,In_gene,Glial_gene)
mk_gene <- human_to_mouse_anno[mk_gene]
names(mk_gene) <- c(rep('Progenitor',times = length(Progenitor_gene)),rep('Ex',times = length(Ex_gene)),
                    rep('In',times = length(In_gene)),rep('Glial',times = length(Glial_gene)))

#aggregate express matrix
express_matrix <- AverageExpression(object = express_matrix,assays = 'RNA',return.seurat = FALSE,group.by = group.by,slot = 'data',verbose = TRUE)
express_matrix <- express_matrix$RNA[mk_gene,]
express_matrix <- t(scale(t(express_matrix)))
express_matrix <- express_matrix[mk_gene,cell_type_list]

#add annotation
top_anno <- HeatmapAnnotation(cell_type = colnames(express_matrix),
                              col = list(cell_type = color_param$celltype[cell_type_list]),
                              show_legend = FALSE,which = 'column',border = FALSE,show_annotation_name = FALSE)

col_fun <- colorRamp2(breaks = c(-2.5,0,2.5),colors = c('lightgrey','white',color_param$species[species_list]))

#plot
h_mouse <- Heatmap(matrix = express_matrix,show_column_names = FALSE,show_row_names = FALSE,cluster_columns = FALSE,
                   cluster_rows = FALSE,row_split = factor(names(mk_gene),levels = c('Progenitor','Ex','In','Glial')),
                   border = TRUE,top_annotation = top_anno,col = col_fun,name = 'mouse',row_order = h_human@row_order,
                   width = unit(0.125*length(cell_type_list),'inches'),height = unit(0.125*length(mk_gene),'inches'))

## combine plot ------------------------------------------------------------
pdf(file = './res/step_91_fig_221205/species_combined_Heatmap_all_gene.pdf',width = 6.5,height = 8.5)
print(h_human + h_macaque + h_mouse)
dev.off()