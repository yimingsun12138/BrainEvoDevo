#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: genes with amino acid change expression distribution            ##
## Data: 2022.12.08                                                                ##
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

# gene expressed in human -------------------------------------------------
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

#save data
saveRDS(object = gene_list,file = './res/step_91_fig_221208/human_expressed_gene_list.rds')

# GO analysis of human expressed gene -------------------------------------
gene_list <- readRDS(file = './res/step_91_fig_221208/human_expressed_gene_list.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

## GO enrichment -----------------------------------------------------------
#BP
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human expressed gene BP analysis",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[1:10,]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
p_BP <- ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = pal_d3(alpha = 0.5)(3)[1],colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.5,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('BP')

#CC
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human expressed gene CC analysis",
                 ontology = "CC",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[1:10,]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
p_CC <- ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = pal_d3(alpha = 0.5)(3)[2],colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.5,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('CC')

#MF
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human expressed gene MF analysis",
                 ontology = "MF",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!(duplicated(allRes$Term)),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[1:10,]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
p_MF <- ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = pal_d3(alpha = 0.5)(3)[3],colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.5,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('MF')

#combined plot
pdf(file = './res/step_91_fig_221208/human_expressed_gene_GO_barplot.pdf',width = 8,height = 8)
p_BP + p_CC + p_MF + plot_layout(ncol = 1)
dev.off()

## GO structure ------------------------------------------------------------
#BP
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human expressed gene BP analysis",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot
pdf(file = './res/step_91_fig_221208/human_expressed_gene_BP_GO_structure.pdf',width = 4,height = 3)
showSigOfNodes(GOdata = GO_enrich,termsP.value = score(resultFisher),firstSigNodes = 5,useInfo = 'all')
dev.off()

# gene number enriched in each cell type ----------------------------------
#load data
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
gene_list <- readRDS(file = './res/step_91_fig_221208/human_expressed_gene_list.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_RNA_Seurat$ReAnno_celltype]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

temp <- c('Progenitor','Progenitor','Progenitor','Progenitor','Ex','Ex','Ex','Ex','In','In','Glial','Glial')
names(temp) <- cell_type_list
cell_type_list <- temp

#generate cell type marker
group_list <- c('Progenitor','Ex','In','Glial')
marker_table <- do.call(what = rbind,args = base::lapply(X = group_list,FUN = function(x){
  #get cell type
  cell_type_1 <- names(cell_type_list)[cell_type_list == x]
  cell_type_2 <- names(cell_type_list)[cell_type_list != x]
  
  #get express matrix
  cell_type_1 <- Greenleaf_RNA_Seurat@assays$RNA@data[gene_list,Greenleaf_RNA_Seurat$ReAnno_celltype %in% cell_type_1]
  cell_type_2 <- Greenleaf_RNA_Seurat@assays$RNA@data[gene_list,Greenleaf_RNA_Seurat$ReAnno_celltype %in% cell_type_2]
  cell_type_1 <- expm1(cell_type_1)
  cell_type_2 <- expm1(cell_type_2)
  
  #do wilcox
  DF_matrix <- my_DF_wilcox_test(mat1 = cell_type_1,mat2 = cell_type_2,alternative = 'greater',paired = FALSE,workers = 1,future.globals.maxSize = 200*(1024^3))
  
  #modify DF matrix
  DF_matrix$gene <- rownames(DF_matrix)
  DF_matrix$group <- x
  
  #return
  return(DF_matrix)
}))

#get marker gene list
mk_gene <- base::lapply(X = group_list,FUN = function(x){
  DF_matrix <- marker_table[marker_table$group == x,]
  temp <- which(DF_matrix$fdr < 0.01 & DF_matrix$log2FC > 0.5)
  if(length(temp) == 0){
    return(NA)
  }else{
    return(DF_matrix$gene[temp])
  }
})
names(mk_gene) <- group_list

#barplot
temp <- unlist(base::lapply(X = group_list,FUN = function(x){
  temp <- mk_gene[[x]]
  return(length(temp))
}))
temp <- data.frame(group = group_list,num = temp)
temp$group <- factor(temp$group,levels = group_list)

p <- ggplot(data = temp,aes(x = group,y = num)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'lightgrey',color = 'black') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1)) + 
  xlab('cell type category') + ylab('marker gene number')

pdf(file = './res/step_91_fig_221208/human_expressed_gene_cell_type_marker_num_barplot.pdf',width = 3,height = 3)
print(p)
dev.off()

# GO analysis of different group of marker gene ---------------------------
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
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

## Progenitor --------------------------------------------------------------
#set param
gene_list <- Progenitor_gene

#GO
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human Progenitor BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

allRes <- allRes[1:5,]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
p_Progenitor <- ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'lightgrey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('Progenitor BP')

## Ex --------------------------------------------------------------
#set param
gene_list <- Ex_gene

#GO
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human Ex BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

allRes <- allRes[1:5,]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
p_Ex <- ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'lightgrey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('Ex BP')

## Glial --------------------------------------------------------------
#set param
gene_list <- Glial_gene

#GO
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human Glial BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

allRes <- allRes[c(1,2,3,5,6),]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
p_Glial <- ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'lightgrey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('Glial BP')

## Global --------------------------------------------------------------
#set param
gene_list <- Global_gene

#GO
human_all_gene <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human Global BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", 
                   numChar = 100, topNodes = length(resultFisher@score))
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

allRes <- allRes[1:5,]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
p_Global <- ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'lightgrey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('Global BP')

## combine plot ------------------------------------------------------------
pdf(file = './res/step_91_fig_221208/cell_type_marker_BP_barplot.pdf',width = 10,height = 8)
p_Progenitor + p_Ex + p_Glial + p_Global + plot_layout(ncol = 1)
dev.off()

# genes worth further investigate -----------------------------------------
Progenitor_gene <- c('COL11A1','AKR7A2','PGM1','CHEK1','LAMA1')
Ex_gene <- c('CSMD2','PTPN4','SLC8A1','TENM4')

#change dir
setwd(dir = '/home/sunym/temp')

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = '/home/sunym/temp/merged_Greenleaf_ATAC_ArchR_221018/')
macaque_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_multiome_ArchR_221011/')
mouse_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/mouse_multiome_ArchR_221009/')

human_to_macaque_anno <- readRDS(file = '/content/data/sunym/project/Brain/res/step_91_fig_221205/macaque_homology_gene_list.rds')
human_to_macaque_anno <- c(human_to_macaque_anno,c(LAMA1 = 'LAMA1'))
human_to_mouse_anno <- readRDS(file = '/content/data/sunym/project/Brain/res/step_91_fig_221205/mouse_homology_gene_list.rds')
human_to_mouse_anno <- c(human_to_mouse_anno,c(AKR7A2 = 'Akr7a5'))

#impute
Greenleaf_ATAC_ArchR <- addImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR)
macaque_multiome_ArchR <- addImputeWeights(ArchRProj = macaque_multiome_ArchR)
mouse_multiome_ArchR <- addImputeWeights(ArchRProj = mouse_multiome_ArchR)

#marker gene
mk_gene <- c(Progenitor_gene,Ex_gene,Glial_gene)
mk_gene_macaque <- human_to_macaque_anno[mk_gene]
mk_gene_mouse <- human_to_mouse_anno[mk_gene]

mk_gene_macaque <- mk_gene_macaque[!(mk_gene_macaque %in% c('PGM1','LAMA1'))]

#featureplot
p_human <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,embedding = 'UMAP',colorBy = 'GeneScoreMatrix',name = mk_gene,
                         imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),size = 0.1,plotAs = 'points')
p_macaque <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,embedding = 'UMAP',colorBy = 'GeneScoreMatrix',name = mk_gene_macaque,
                           imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),size = 0.1,plotAs = 'points')
p_mouse <- plotEmbedding(ArchRProj = mouse_multiome_ArchR,embedding = 'UMAP',colorBy = 'GeneScoreMatrix',name = mk_gene_mouse,
                         imputeWeights = getImputeWeights(ArchRProj = mouse_multiome_ArchR),size = 0.1,plotAs = 'points')

#check accessibility

#COL11A1 do trackplot
temp_huamn <- p_human$COL11A1 + theme(legend.position = 'none')
temp_macaque <- p_macaque$COL11A1 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Col11a1 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

#AKR7A2 do trackplot
temp_huamn <- p_human$AKR7A2 + theme(legend.position = 'none')
temp_macaque <- p_macaque$AKR7A2 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Akr7a5 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

#CHEK1 do trackplot
temp_huamn <- p_human$CHEK1 + theme(legend.position = 'none')
temp_macaque <- p_macaque$CHEK1 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Chek1 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

#CSMD2
temp_huamn <- p_human$CSMD2 + theme(legend.position = 'none')
temp_macaque <- p_macaque$CSMD2 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Csmd2 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

#PTPN4 do trackplot
temp_huamn <- p_human$PTPN4 + theme(legend.position = 'none')
temp_macaque <- p_macaque$PTPN4 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Ptpn4 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

#SLC8A1 do trackplot
temp_huamn <- p_human$SLC8A1 + theme(legend.position = 'none')
temp_macaque <- p_macaque$SLC8A1 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Slc8a1 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

#TENM4
temp_huamn <- p_human$TENM4 + theme(legend.position = 'none')
temp_macaque <- p_macaque$TENM4 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Tenm4 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

#ADAM28 do trackplot
temp_huamn <- p_human$ADAM28 + theme(legend.position = 'none')
temp_macaque <- p_macaque$ADAM28 + theme(legend.position = 'none')
temp_mouse <- p_mouse$Adam28 + theme(legend.position = 'none')

temp_huamn + temp_macaque + temp_mouse + plot_layout(ncol = 3)

# trackplot ---------------------------------------------------------------
#load data
human_to_macaque_anno <- readRDS(file = '/content/data/sunym/project/Brain/res/step_91_fig_221205/macaque_homology_gene_list.rds')
human_to_macaque_anno <- c(human_to_macaque_anno,c(LAMA1 = 'LAMA1'))
human_to_mouse_anno <- readRDS(file = '/content/data/sunym/project/Brain/res/step_91_fig_221205/mouse_homology_gene_list.rds')
human_to_mouse_anno <- c(human_to_mouse_anno,c(AKR7A2 = 'Akr7a5'))

#gene list
mk_gene <- c('COL11A1','AKR7A2','PTPN4','SLC8A1')
mk_gene_macaque <- human_to_macaque_anno[mk_gene]
mk_gene_mouse <- human_to_mouse_anno[mk_gene]

#plot function
plot_track <- function(gene_list = gene_list,species = 'human'){
  #load data
  macaque_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_multiome_ArchR_221011/')
  Greenleaf_ATAC_ArchR <- loadArchRProject(path = '/home/sunym/temp/merged_Greenleaf_ATAC_ArchR_221018/')
  mouse_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/mouse_multiome_ArchR_221009/')
  
  color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
  human_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
  human_anno <- rtracklayer::as.data.frame(x = human_anno)
  macaque_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
  macaque_anno <- rtracklayer::as.data.frame(x = macaque_anno)
  mouse_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Mus_musculus.GRCm38.102.gtf',format = 'gtf')
  mouse_anno <- rtracklayer::as.data.frame(x = mouse_anno)
  
  #set parameter
  cell_type_list <- c('RG-1','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE')
  
  cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)
  extend_size <- 100000
  ylim_track <- 0.03
  
  #load bw
  if(species == 'human'){
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type/')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }else if(species == 'macaque'){
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }else{
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/mouse/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }
  
  #for loop
  for (i in gene_list) {
    #set gene name
    gene_name <- i
    
    #get target site
    if(species == 'human'){
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = Greenleaf_ATAC_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }else if(species == 'macaque'){
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = macaque_multiome_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }else{
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = mouse_multiome_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }
    
    TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]
    if(length(TSS_site) > 1){
      stop('TSS number error!')
    }
    chrom <- as.character(TSS_site@seqnames)
    start_site <- TSS_site@ranges@start - extend_size
    end_site <- TSS_site@ranges@start + extend_size
    
    #plot gene annotation
    if(species == 'human'){
      p_gene_anno <- trancriptVis(gtfFile = human_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }else if(species == 'macaque'){
      p_gene_anno <- trancriptVis(gtfFile = macaque_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }else{
      p_gene_anno <- trancriptVis(gtfFile = mouse_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }
    
    #plot peak
    if(species == 'human'){
      temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
    }else if(species == 'macaque'){
      temp <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
    }else{
      temp <- getPeakSet(ArchRProj = mouse_multiome_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
    }
    
    p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed'),
                          chr = chrom,region.min = start_site - 10000,region.max = end_site + 10000,
                          show.legend = FALSE,fill = ggsci::pal_d3()(1),track.width = 0.3)
    
    #plot track
    if(species == 'human'){
      temp_anno <- human_anno
    }else if(species == 'macaque'){
      temp_anno <- macaque_anno
    }else{
      temp_anno <- mouse_anno
    }
    
    p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = temp_anno,
                             chr = chrom,region.min = start_site,region.max = end_site,
                             sample.order = cell_type_list,space.y = 0,
                             y.max = ylim_track,color = color_param$celltype[cell_type_list],
                             theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)
    
    #save plot
    char <- paste0('./res/step_91_fig_221208/',species,'_',gene_name,'_trackplot.pdf')
    pdf(file = char,width = 8,height = 8)
    print(p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.05) %>% insert_bottom(plot = p_gene_anno,height = 0.05))
    dev.off()
  }
  return(NULL)
}

#plot
plot_track(gene_list = mk_gene,species = 'human')
gc()
plot_track(gene_list = mk_gene_macaque,species = 'macaque')
gc()
plot_track(gene_list = mk_gene_mouse,species = 'mouse')
gc()

my_send_sms('plot done!')

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}