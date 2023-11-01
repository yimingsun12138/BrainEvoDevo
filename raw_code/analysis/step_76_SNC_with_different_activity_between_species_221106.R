#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: SNC with different activity between species                     ##
## Data: 2022.10.23                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# active in human ---------------------------------------------------------

## modify human SNC GRannges -----------------------------------------------

#load data
human_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_human_SNC.rds')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')


#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]

#for loop
for (i in cell_type_list) {
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #get human DF
  DF_list <- list.files(path = './res/step_72_fig_221101/DF_list/')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_72_fig_221101/DF_list',DF_list,sep = '/')
  print(file.exists(DF_list))
  DF_list <- readRDS(file = DF_list)
  
  DF_list <- DF_list[DF_list$group == 'human_specific',]
  DF_list <- DF_list$human_SNC
  
  #add meta data
  human_SNC@elementMetadata[,cell_type_dot] <- 'NO'
  human_SNC@elementMetadata[human_SNC$human_SNC %in% DF_list,cell_type_dot] <- 'YES'
}

## human active SNC cell type distribution ---------------------------------
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)
temp <- data.frame(cell_type = cell_type_list,
                   SNC = unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
                     temp <- as.character(human_SNC@elementMetadata[,x])
                     return(sum(temp == 'YES'))
                   })))

temp$cell_type <- factor(temp$cell_type,levels = rev(temp$cell_type))

#plot
pdf(file = './res/step_76_fig_221106/human_active_SNC_cell_type_enrichment.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x = cell_type,y = SNC,fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8) + 
  scale_fill_manual(values = color_param$celltype[as.character(temp$cell_type)]) + 
  theme_bw() + 
  ylim(c(0,300)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(linewidth = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(linewidth = 0.8),
        aspect.ratio = 2.5) + 
  xlab('Cell type') + ylab('human DA SNC number') + 
  coord_flip()
dev.off()

## human SNC cell type enrichment intersect size ---------------------------
temp <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp <- human_SNC[human_SNC@elementMetadata[,x] == 'YES']
  temp <- temp$human_SNC
  return(temp)
})

names(temp) <- cell_type_list

#upset plot
pdf(file = './res/step_76_fig_221106/human_active_SNC_cell_type_enrichment_intersect_upsetplot.pdf',width = 12,height = 7)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#intersect size
temp <- rtracklayer::as.data.frame(x = human_SNC)
temp <- temp[,cell_type_list_dot]
rownames(temp) <- human_SNC$human_SNC

temp <- base::lapply(X = human_SNC$human_SNC,FUN = function(x){
  enrich_num <- temp[x,]
  enrich_num <- sum(enrich_num == 'YES')
  return(enrich_num)
})
temp <- unlist(temp)
names(temp) <- human_SNC$human_SNC

temp <- as.data.frame(table(temp))
colnames(temp) <- c('Intersect','percent')
temp$percent <- temp$percent/sum(temp$percent)*100
temp$Intersect <- factor(temp$Intersect,levels = temp$Intersect)

pdf(file = './res/step_76_fig_221106/human_active_SNC_cell_type_enrichment_intersect_size.pdf',width = 5,height = 3)
ggplot(data = temp,aes(x = Intersect,y = percent)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8,color = 'grey') + 
  theme_bw() + 
  scale_y_break(breaks = c(8,82)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),
        legend.position = 'none') + 
  xlab('cell type size') + ylab('Percent %')
dev.off()

## human active SNC genome distribution ------------------------------------
#make Txdb object
gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
gene_anno <- rtracklayer::as.data.frame(gene_anno)
gene_anno$seqnames <- as.character(gene_anno$seqnames)
gene_anno$seqnames <- unlist(base::lapply(X = gene_anno$seqnames,FUN = function(x){
  if(x %in% c(as.character(1:22),'X','Y','MT')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
}))
gene_anno <- makeGRangesFromDataFrame(df = gene_anno,keep.extra.columns = TRUE)
gene_anno <- GenomicFeatures::makeTxDbFromGRanges(gene_anno)

#each cell type
SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp <- human_SNC[human_SNC@elementMetadata[,x] == 'YES']
  return(temp)
})
names(SNC_list) <- cell_type_list

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_76_fig_221106/human_active_SNC_cell_type_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

#different intersect size
temp <- as.data.frame(human_SNC@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
human_SNC$intersect <- temp

SNC_list <- base::lapply(X = c(0:10,12:13),FUN = function(x){
  return(human_SNC[human_SNC$intersect == x])
})
names(SNC_list) <- c(0:10,12:13)

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_76_fig_221106/human_active_SNC_cell_type_intersect_size_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

## overlap with MPRA results -----------------------------------------------

#load data
human_peak_matrix <- readRDS(file = './res/step_72_fig_221101/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_72_fig_221101/macaque_peak_matrix.rds')

human_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_macaque_SNC.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

#do wilcox
human_cell_list <- rownames(Greenleaf_ATAC_ArchR@cellColData)[Greenleaf_ATAC_ArchR$cell_type %in% c('RG-1','RG-2','IP')]
macaque_cell_list <- rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$cell_type %in% c('RG-1','RG-2','IP')]

norm_factor <- c(Greenleaf_ATAC_ArchR@cellColData[human_cell_list,"ReadsInPeaks"],macaque_multiome_ArchR@cellColData[macaque_cell_list,"ReadsInPeaks"])
norm_factor <- median(norm_factor) * 10^4

subset_human_peak_matrix <- human_peak_matrix[,human_cell_list]
subset_human_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = human_cell_list,FUN = function(x){
  temp <- subset_human_peak_matrix[,x]/Greenleaf_ATAC_ArchR@cellColData[x,"ReadsInPeaks"]*norm_factor
  return(temp)
}))
colnames(subset_human_peak_matrix) <- human_cell_list
rownames(subset_human_peak_matrix) <- rownames(human_peak_matrix)

subset_macaque_peak_matrix <- macaque_peak_matrix[,macaque_cell_list]
subset_macaque_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = macaque_cell_list,FUN = function(x){
  temp <- subset_macaque_peak_matrix[,x]/macaque_multiome_ArchR@cellColData[x,"ReadsInPeaks"]*norm_factor
  return(temp)
}))
colnames(subset_macaque_peak_matrix) <- macaque_cell_list
rownames(subset_macaque_peak_matrix) <- rownames(macaque_peak_matrix)

names(human_SNC) <- paste(human_SNC@seqnames,as.character(human_SNC@ranges),sep = '-')
names(macaque_SNC) <- paste(macaque_SNC@seqnames,as.character(macaque_SNC@ranges),sep = '-')

rownames(subset_human_peak_matrix) <- human_SNC[rownames(subset_human_peak_matrix)]$human_SNC
rownames(subset_macaque_peak_matrix) <- macaque_SNC[rownames(subset_macaque_peak_matrix)]$human_SNC

subset_macaque_peak_matrix <- subset_macaque_peak_matrix[rownames(subset_human_peak_matrix),]

DF_list <- my_sparseMatWilcoxon(mat1 = subset_human_peak_matrix,mat2 = subset_macaque_peak_matrix)
DF_list$human_SNC <- rownames(DF_list)

#human_peak_file
RG.1_peak_file <- readRDS(file = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls/RG.1-reproduciblePeaks.gr.rds')
RG.2_peak_file <- readRDS(file = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls/RG.2-reproduciblePeaks.gr.rds')
IP_peak_file <- readRDS(file = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls/IP-reproduciblePeaks.gr.rds')
human_peak_file <- append(RG.1_peak_file,RG.2_peak_file)
human_peak_file <- append(human_peak_file,IP_peak_file)

#macaque_peak_file
RG.1_peak_file <- readRDS(file = './processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls/RG.1-reproduciblePeaks.gr.rds')
RG.2_peak_file <- readRDS(file = './processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls/RG.2-reproduciblePeaks.gr.rds')
IP_peak_file <- readRDS(file = './processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls/IP-reproduciblePeaks.gr.rds')
macaque_peak_file <- append(RG.1_peak_file,RG.2_peak_file)
macaque_peak_file <- append(macaque_peak_file,IP_peak_file)

#human specific
temp <- DF_list[DF_list$log2FC > 1 & DF_list$fdr < 0.01,]
names(human_SNC) <- human_SNC$human_SNC
temp <- human_SNC[temp$human_SNC]
temp <- countOverlaps(query = temp,subject = human_peak_file)
temp <- names(temp)[temp > 0]

DF_list[DF_list$human_SNC %in% temp,'group'] <- 'human_specific'

#macaque specific
temp <- DF_list[DF_list$log2FC < -1 & DF_list$fdr < 0.01,]
names(macaque_SNC) <- macaque_SNC$human_SNC
temp <- macaque_SNC[temp$human_SNC]
temp <- countOverlaps(query = temp,subject = macaque_peak_file)
temp <- names(temp)[temp > 0]

DF_list[DF_list$human_SNC %in% temp,'group'] <- 'macaque_specific'

#species_conserved
temp <- DF_list[abs(DF_list$log2FC) <= 1,]
names(human_SNC) <- human_SNC$human_SNC
temp <- human_SNC[temp$human_SNC]
temp <- countOverlaps(query = temp,subject = human_peak_file)
temp <- names(temp)[temp > 0]
temp_human <- temp

temp <- DF_list[abs(DF_list$log2FC) <= 1,]
names(macaque_SNC) <- macaque_SNC$human_SNC
temp <- macaque_SNC[temp$human_SNC]
temp <- countOverlaps(query = temp,subject = macaque_peak_file)
temp <- names(temp)[temp > 0]
temp_macaque <- temp

temp <- dplyr::intersect(x = temp_human,y = temp_macaque)

DF_list[DF_list$human_SNC %in% temp,"group"] <- 'species_conserved'

#add meta data
DF_list <- DF_list[!is.na(DF_list$group),]
table(DF_list$group)

names(human_SNC) <- human_SNC$human_SNC
human_SNC$group <- NA
human_SNC[DF_list$human_SNC]@elementMetadata[,"group"] <- DF_list$group

#venn plot of MPRA human active and ATAC human active
MPRA_human_specific <- human_SNC[human_SNC$NPC_differentially_active == 'Yes' & human_SNC$NPC_exp_lfc_modern_vs_archaic > 0]$human_SNC
ATAC_human_specific <- human_SNC[which(human_SNC$group == 'human_specific')]$human_SNC

pdf(file = './res/step_76_fig_221106/MPRA_ATAC_human_specific_SNC_overlap.pdf',width = 5,height = 5)
ggvenn(data = list(MPRA = MPRA_human_specific,ATAC = ATAC_human_specific),columns = c('MPRA','ATAC'))
dev.off()

#dot plot
human_SNC[is.na(human_SNC$group)]$group <- 'others'
table(human_SNC$group)

names(human_SNC) <- human_SNC$human_SNC
DF_list_temp <- my_sparseMatWilcoxon(mat1 = subset_human_peak_matrix,mat2 = subset_macaque_peak_matrix)
table(names(human_SNC) %in% rownames(DF_list_temp))

human_SNC <- human_SNC[rownames(DF_list_temp)]
human_SNC$ATAC_lfc <- DF_list_temp$log2FC
human_SNC <- human_SNC[!is.na(human_SNC$NPC_exp_lfc_modern_vs_archaic)]

temp <- as.data.frame(human_SNC@elementMetadata)

pdf(file = './res/step_76_fig_221106/MPRA_ATAC_SNC_NPC_lfc_dot_plot.pdf',width = 5,height = 5)
ggplot(data = temp[temp$group != 'others',],aes(x = NPC_exp_lfc_modern_vs_archaic, y = ATAC_lfc, color = group)) + 
  geom_point(size = 0.1) + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 12,face = 'bold',hjust = 0.5)) + 
  guides(color = guide_legend(override.aes = list(size=2))) + 
  labs(title = 'NPC')
dev.off()

# active in all species ---------------------------------------------------

## modify human SNC GRanges ------------------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_human_SNC.rds')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')


#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]

#for loop
for (i in cell_type_list) {
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #get human DF
  DF_list <- list.files(path = './res/step_72_fig_221101/DF_list/')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_72_fig_221101/DF_list',DF_list,sep = '/')
  print(file.exists(DF_list))
  DF_list <- readRDS(file = DF_list)
  
  DF_list <- DF_list[DF_list$group == 'species_conserved',]
  DF_list <- DF_list$human_SNC
  
  #add meta data
  human_SNC@elementMetadata[,cell_type_dot] <- 'NO'
  human_SNC@elementMetadata[human_SNC$human_SNC %in% DF_list,cell_type_dot] <- 'YES'
}

## all active SNC cell type enrichment distribution ------------------------
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)
temp <- data.frame(cell_type = cell_type_list,
                   SNC = unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
                     temp <- as.character(human_SNC@elementMetadata[,x])
                     return(sum(temp == 'YES'))
                   })))

temp$cell_type <- factor(temp$cell_type,levels = rev(temp$cell_type))

#plot
pdf(file = './res/step_76_fig_221106/all_active_SNC_cell_type_enrichment.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x = cell_type,y = SNC,fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8) + 
  scale_fill_manual(values = color_param$celltype[as.character(temp$cell_type)]) + 
  theme_bw() + 
  ylim(c(0,350)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(linewidth = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(linewidth = 0.8),
        aspect.ratio = 2.5) + 
  xlab('Cell type') + ylab('all active SNC number') + 
  coord_flip()
dev.off()

## all active SNC cell type enrichment intersect size ----------------------
temp <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp <- human_SNC[human_SNC@elementMetadata[,x] == 'YES']
  temp <- temp$human_SNC
  return(temp)
})

names(temp) <- cell_type_list

#upset plot
pdf(file = './res/step_76_fig_221106/all_active_SNC_cell_type_enrichment_intersect_upsetplot.pdf',width = 12,height = 7)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#intersect size
temp <- rtracklayer::as.data.frame(x = human_SNC)
temp <- temp[,cell_type_list_dot]
rownames(temp) <- human_SNC$human_SNC

temp <- base::lapply(X = human_SNC$human_SNC,FUN = function(x){
  enrich_num <- temp[x,]
  enrich_num <- sum(enrich_num == 'YES')
  return(enrich_num)
})
temp <- unlist(temp)
names(temp) <- human_SNC$human_SNC

temp <- as.data.frame(table(temp))
colnames(temp) <- c('Intersect','percent')
temp$percent <- temp$percent/sum(temp$percent)*100
temp$Intersect <- factor(temp$Intersect,levels = temp$Intersect)

pdf(file = './res/step_76_fig_221106/all_active_SNC_cell_type_enrichment_intersect_size.pdf',width = 5,height = 3)
ggplot(data = temp,aes(x = Intersect,y = percent)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8,color = 'grey') + 
  theme_bw() + 
  scale_y_break(breaks = c(5,87)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(linewidth = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(linewidth = 0.8),
        legend.position = 'none') + 
  xlab('cell type size') + ylab('Percent %')
dev.off()

## all active SNC genome distribution --------------------------------------
#make Txdb object
gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
gene_anno <- rtracklayer::as.data.frame(gene_anno)
gene_anno$seqnames <- as.character(gene_anno$seqnames)
gene_anno$seqnames <- unlist(base::lapply(X = gene_anno$seqnames,FUN = function(x){
  if(x %in% c(as.character(1:22),'X','Y','MT')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
}))
gene_anno <- makeGRangesFromDataFrame(df = gene_anno,keep.extra.columns = TRUE)
gene_anno <- GenomicFeatures::makeTxDbFromGRanges(gene_anno)

#each cell type
SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp <- human_SNC[human_SNC@elementMetadata[,x] == 'YES']
  return(temp)
})
names(SNC_list) <- cell_type_list

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_76_fig_221106/all_active_SNC_cell_type_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

#different intersect size
temp <- as.data.frame(human_SNC@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
human_SNC$intersect <- temp

SNC_list <- base::lapply(X = c(0:13),FUN = function(x){
  return(human_SNC[human_SNC$intersect == x])
})
names(SNC_list) <- c(0:13)

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_76_fig_221106/all_active_SNC_cell_type_intersect_size_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

# human active SNC posit at regulation region -----------------------------
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

#human active
human_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_human_SNC.rds')
for (i in cell_type_list) {
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #get human DF
  DF_list <- list.files(path = './res/step_72_fig_221101/DF_list/')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_72_fig_221101/DF_list',DF_list,sep = '/')
  print(file.exists(DF_list))
  DF_list <- readRDS(file = DF_list)
  
  DF_list <- DF_list[DF_list$group == 'human_specific',]
  DF_list <- DF_list$human_SNC
  
  #add meta data
  human_SNC@elementMetadata[,cell_type_dot] <- 'NO'
  human_SNC@elementMetadata[human_SNC$human_SNC %in% DF_list,cell_type_dot] <- 'YES'
}
human_active <- human_SNC

#all active
human_SNC <- readRDS(file = './res/step_72_fig_221101/macaque_lifted_SNC_after_seq_filter/extented_human_SNC.rds')
for (i in cell_type_list) {
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #get human DF
  DF_list <- list.files(path = './res/step_72_fig_221101/DF_list/')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_72_fig_221101/DF_list',DF_list,sep = '/')
  print(file.exists(DF_list))
  DF_list <- readRDS(file = DF_list)
  
  DF_list <- DF_list[DF_list$group == 'species_conserved',]
  DF_list <- DF_list$human_SNC
  
  #add meta data
  human_SNC@elementMetadata[,cell_type_dot] <- 'NO'
  human_SNC@elementMetadata[human_SNC$human_SNC %in% DF_list,cell_type_dot] <- 'YES'
}
all_active <- human_SNC

#make Txdb object
gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
gene_anno <- rtracklayer::as.data.frame(gene_anno)
gene_anno$seqnames <- as.character(gene_anno$seqnames)
gene_anno$seqnames <- unlist(base::lapply(X = gene_anno$seqnames,FUN = function(x){
  if(x %in% c(as.character(1:22),'X','Y','MT')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
}))
gene_anno <- makeGRangesFromDataFrame(df = gene_anno,keep.extra.columns = TRUE)
gene_anno <- GenomicFeatures::makeTxDbFromGRanges(gene_anno)

#boxplot
#human
SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  return(human_active[human_active@elementMetadata[,x] == 'YES'])
})
names(SNC_list) <- cell_type_list

peakAnno_human <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
peakAnno_human <- do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- as.data.frame(peakAnno_human[[x]]@annoStat)
  promoter_percent <- sum(temp[grepl(pattern = 'Promoter',x = temp$Feature,fixed = TRUE),'Frequency'])
  others_percent <- sum(temp[!grepl(pattern = 'Promoter',x = temp$Feature,fixed = TRUE),'Frequency'])
  temp <- data.frame(Feature = c('Promoter','Others'),
                     Percent = c(promoter_percent,others_percent),
                     cell_type = x,group = 'human active')
  return(temp)
}))

#all
SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  return(all_active[all_active@elementMetadata[,x] == 'YES'])
})
names(SNC_list) <- cell_type_list

peakAnno_all <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
peakAnno_all <- do.call(what = rbind,args = base::lapply(X = cell_type_list,FUN = function(x){
  temp <- as.data.frame(peakAnno_all[[x]]@annoStat)
  promoter_percent <- sum(temp[grepl(pattern = 'Promoter',x = temp$Feature,fixed = TRUE),'Frequency'])
  others_percent <- sum(temp[!grepl(pattern = 'Promoter',x = temp$Feature,fixed = TRUE),'Frequency'])
  temp <- data.frame(Feature = c('Promoter','Others'),
                     Percent = c(promoter_percent,others_percent),
                     cell_type = x,group = 'all active')
  return(temp)
}))

#plot
peakAnno <- rbind(peakAnno_human,peakAnno_all)
peakAnno <- peakAnno[peakAnno$Feature == 'Promoter',]
peakAnno$cell_type <- factor(peakAnno$cell_type,levels = cell_type_list)
peakAnno$group <- factor(peakAnno$group,levels = c('human active','all active'))

pdf(file = './res/step_76_fig_221106/human_active_posit_at_regulation_region.pdf',width = 3,height = 5)
ggplot(data = peakAnno,aes(x = group,y = Percent,fill = group)) + 
  geom_boxplot(outlier.alpha = 0,size = 1) + 
  geom_point(aes(color = cell_type),size = 0.5,alpha = 0.5) + 
  geom_path(aes(group = cell_type,color = 'lightgrey'),alpha = 0.5) + 
  scale_color_manual(values = color_param$celltype[cell_type_list]) + 
  scale_fill_manual(values = c('human active' = '#3FA4D9','all active' = '#B2C224')) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1.8,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        legend.position = 'none') + 
  ylab('Promoter percent %')
dev.off()

#merge
temp <- as.data.frame(human_active@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
human_active$intersect <- temp
human_active <- human_active[human_active$intersect != 0,]

temp <- as.data.frame(all_active@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
all_active$intersect <- temp
all_active <- all_active[all_active$intersect != 0,]

SNC_list <- list(`human active` = human_active,
                 `all active` = all_active)

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_76_fig_221106/merged_human_active_all_active_SNC_genome_distribution.pdf',width = 6,height = 3)
p + theme(aspect.ratio = 0.3)
dev.off()