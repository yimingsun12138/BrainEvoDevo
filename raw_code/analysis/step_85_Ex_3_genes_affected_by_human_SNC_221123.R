#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Ex-3 genes affected by human SNC                                ##
## Data: 2022.11.23                                                                ##
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

# filter Ex-3 gene --------------------------------------------------------
cell_type <- 'Ex-3'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
sequence_list <- readRDS(file = './res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')

affinity_difference <- readRDS(file = './res/step_77_fig_221108/affinity_difference_data_frame.rds')
human_fimo <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
motif_info <- read.csv(file = './data/reference/motif_DB/CISBP/human_221108/TF_information.csv')

human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')

#load gene list
cell_type_dot <- gsub(pattern = '-',replacement = '.',fixed = TRUE,x = cell_type)
gene_list <- list.files(path = './res/step_79_fig_221112/up_gene_regulated_by_human_specific_SNC')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/up_gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(file = gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#modify affinity difference
affinity_difference$human_SNC <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))

#filter gene list
gene_list <- base::lapply(X = gene_list,FUN = function(x){
  #filter human SNC overlap with p2g
  temp <- which(p2g@metadata$geneSet$name == x)
  temp <- which(p2g$idxRNA %in% temp)
  temp <- p2g$idxATAC[temp]
  temp <- p2g@metadata$peakSet[temp]
  
  names(human_SNC) <- human_SNC$human_SNC
  SNC_list <- countOverlaps(query = human_SNC,subject = temp)
  SNC_list <- names(SNC_list)[SNC_list > 0]
  
  #see the affinity change
  temp <- affinity_difference[affinity_difference$human_SNC %in% SNC_list & affinity_difference$cell_type == cell_type,]
  
  #return
  if(nrow(temp) > 0){
    return(x)
  }else{
    return(NA)
  }
})

gene_list <- unlist(gene_list)
gene_list <- gene_list[!is.na(gene_list)]

#save gene list
saveRDS(object = gene_list,file = './res/step_85_fig_221123/gene_list.rds')

# RTN4RL2 -----------------------------------------------------------------
#load param
cell_type <- 'Ex-3'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
gene_name <- 'RTN4RL2'
extend_size <- 100000
resolution <- 1000
ylim_track <- 0.03

## RNA expression ----------------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
Idents(macaque_multiome_Seurat) <- 'cell_type'
Idents(Greenleaf_RNA_Seurat) <- 'ReAnno_celltype'

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

char <- paste0('./res/step_85_fig_221123/',gene_name,'_human_macaque_RNA_express_featureplot.pdf')
pdf(file = char,width = 5,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

## see motif changes -------------------------------------------------------
#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
sequence_list <- readRDS(file = './res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')

#load affinity change
affinity_difference <- readRDS(file = './res/step_77_fig_221108/affinity_difference_data_frame.rds')
human_fimo <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
motif_info <- read.csv(file = './data/reference/motif_DB/CISBP/human_221108/TF_information.csv')

#find human SNC related to CSTB
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

temp <- which(p2g@metadata$geneSet$name == gene_name)
print(paste('gene length:',length(temp),sep = ' '))
temp <- which(p2g$idxRNA %in% temp)
temp <- p2g$idxATAC[temp]
temp <- p2g@metadata$peakSet[temp]

names(human_SNC) <- human_SNC$human_SNC
SNC_list <- countOverlaps(query = human_SNC,subject = temp)
SNC_list <- names(SNC_list)[SNC_list > 0]

#see the affinity change
affinity_difference$human_SNC <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))
affinity_difference <- affinity_difference[affinity_difference$human_SNC %in% SNC_list & affinity_difference$cell_type == cell_type,]
SNC_list <- SNC_list[SNC_list %in% affinity_difference$human_SNC]
print(paste('SNC list length:',length(SNC_list),sep = ' '))

#motif
motif_name <- 'M07764_2.00'
motif_1 <- list.files(path = './data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_1 <- motif_1[grep(pattern = motif_name,x = motif_1,fixed = TRUE)]
motif_1 <- paste('./data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_1,sep = '/')
motif_1 <- read_cisbp(file = motif_1)
motif_1@name <- motif_name

char <- paste0('./res/step_85_fig_221123/',motif_name,'_motif_view.pdf')
pdf(file = char,width = 12,height = 2.5)
view_motifs(motifs = motif_1,show.names = TRUE,names.pos = 'right')
dev.off()

#seems like a dead end

# ATP1A3 -----------------------------------------------------------------
#load param
cell_type <- 'Ex-3'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
gene_name <- 'ATP1A3'
extend_size <- 100000
resolution <- 1000
ylim_track <- 0.03

## RNA expression ----------------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
Idents(macaque_multiome_Seurat) <- 'cell_type'
Idents(Greenleaf_RNA_Seurat) <- 'ReAnno_celltype'

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

char <- paste0('./res/step_85_fig_221123/',gene_name,'_human_macaque_RNA_express_featureplot.pdf')
pdf(file = char,width = 5,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

## see motif changes -------------------------------------------------------
#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
sequence_list <- readRDS(file = './res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')

#load affinity change
affinity_difference <- readRDS(file = './res/step_77_fig_221108/affinity_difference_data_frame.rds')
human_fimo <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
motif_info <- read.csv(file = './data/reference/motif_DB/CISBP/human_221108/TF_information.csv')

#find human SNC related to CSTB
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

temp <- which(p2g@metadata$geneSet$name == gene_name)
print(paste('gene length:',length(temp),sep = ' '))
temp <- which(p2g$idxRNA %in% temp)
temp <- p2g$idxATAC[temp]
temp <- p2g@metadata$peakSet[temp]

names(human_SNC) <- human_SNC$human_SNC
SNC_list <- countOverlaps(query = human_SNC,subject = temp)
SNC_list <- names(SNC_list)[SNC_list > 0]

#see the affinity change
affinity_difference$human_SNC <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))
affinity_difference <- affinity_difference[affinity_difference$human_SNC %in% SNC_list & affinity_difference$cell_type == cell_type,]
SNC_list <- SNC_list[SNC_list %in% affinity_difference$human_SNC]
print(paste('SNC list length:',length(SNC_list),sep = ' '))

#motif
motif_name <- 'M07764_2.00'
motif_1 <- list.files(path = './data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_1 <- motif_1[grep(pattern = motif_name,x = motif_1,fixed = TRUE)]
motif_1 <- paste('./data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_1,sep = '/')
motif_1 <- read_cisbp(file = motif_1)
motif_1@name <- motif_name

char <- paste0('./res/step_85_fig_221123/',motif_name,'_motif_view.pdf')
pdf(file = char,width = 12,height = 2.5)
view_motifs(motifs = motif_1,show.names = TRUE,names.pos = 'right')
dev.off()

#seems like a dead end