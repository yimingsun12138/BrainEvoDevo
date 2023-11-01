#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: SNC with different activity between species                     ##
## Data: 2022.10.23                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#sleep
ii <- 1
while(1){
  i <- matrix(data = rnorm(n = 2000*2000),nrow = 2000,ncol = 2000)
  i <- i %*% i
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

# load data ---------------------------------------------------------------
SNC_human <- readRDS(file = './data/public/The_cis_regulatory_effects_of_modern_human_specific_variants/processed_data/SNC_GRanges.rds')
SNC_macaque <- readRDS(file = './res/step_69_fig_221022/SNC_macaque.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

# modify SNC --------------------------------------------------------------
names(SNC_human) <- paste(SNC_human@seqnames,as.character(SNC_human@ranges),sep = '-')
SNC_human <- SNC_human[SNC_macaque$ori_peak]

cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

#DA in human
for (i in cell_type_list_dot) {
  
  #load differentially active SNC
  DF_SNC <- list.files('./res/step_69_fig_221022/DF_SNC/')
  DF_SNC <- DF_SNC[grep(pattern = i,x = DF_SNC,fixed = TRUE)]
  DF_SNC <- paste('./res/step_69_fig_221022/DF_SNC',DF_SNC,sep = '/')
  DF_SNC <- readRDS(file = DF_SNC)
  DF_SNC <- DF_SNC$human
  
  #add meta data
  SNC_human@elementMetadata[,i] <- 'NO'
  SNC_human[DF_SNC]@elementMetadata[,i] <- 'YES'
}

saveRDS(object = SNC_human,file = './res/step_70_fig_221023/SNC_human_DA_in_human.rds')

#DA in all
for (i in cell_type_list_dot) {
  
  #load differentially active SNC
  DF_SNC <- list.files('./res/step_69_fig_221022/DF_SNC/')
  DF_SNC <- DF_SNC[grep(pattern = i,x = DF_SNC,fixed = TRUE)]
  DF_SNC <- paste('./res/step_69_fig_221022/DF_SNC',DF_SNC,sep = '/')
  DF_SNC <- readRDS(file = DF_SNC)
  DF_SNC <- DF_SNC$all
  
  #add meta data
  SNC_human@elementMetadata[,i] <- 'NO'
  SNC_human[DF_SNC]@elementMetadata[,i] <- 'YES'
}

saveRDS(object = SNC_human,file = './res/step_70_fig_221023/SNC_human_active_in_all.rds')

# DA in human -------------------------------------------------------------
#load data
SNC_human <- readRDS(file = './res/step_70_fig_221023/SNC_human_DA_in_human.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

## cell type enrichment ----------------------------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- data.frame(cell_type = cell_type_list,
                   SNC = unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
                     sum(SNC_human@elementMetadata[,x] == 'YES')
                   })))
temp$cell_type <- factor(temp$cell_type,levels = rev(temp$cell_type))

pdf(file = './res/step_70_fig_221023/human_active_SNC_cell_type_enrichment.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x = cell_type,y = SNC,fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8) + 
  scale_fill_manual(values = color_param$celltype[as.character(temp$cell_type)]) + 
  theme_bw() + 
  ylim(c(0,450)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),
        aspect.ratio = 2.5) + 
  xlab('Cell type') + ylab('human DA SNC number') + 
  coord_flip()
dev.off()

## cell type enrichment intersect ------------------------------------------
#UPset plot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp <- names(SNC_human)[SNC_human@elementMetadata[,x] == 'YES']
  return(temp)
})

names(temp) <- cell_type_list

#plot
pdf(file = './res/step_70_fig_221023/human_active_SNC_cell_type_enrichment_intersect_upsetplot.pdf',width = 12,height = 7)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#overlap number
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- rtracklayer::as.data.frame(x = SNC_human)
temp <- temp[,cell_type_list_dot]
temp <- base::lapply(X = rownames(temp),FUN = function(x){
  return(sum(temp[x,] == 'YES'))
})
temp <- unlist(temp)

temp <- as.data.frame(table(temp))
colnames(temp) <- c('Intersect','percent')
temp$percent <- temp$percent/sum(temp$percent)*100
temp$Intersect <- factor(temp$Intersect,levels = temp$Intersect)

pdf(file = './res/step_70_fig_221023/human_active_SNC_number_with_intersect_size.pdf',width = 5,height = 3)
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

## overlap with NPC human DA SNC -------------------------------------------
temp <- list(NPC_MPRA = names(SNC_human)[SNC_human$NPC_differentially_active == 'Yes'],
             RG_ATAC = names(SNC_human)[SNC_human$RG.1 == 'YES'|SNC_human$RG.2 == 'YES'])

pdf(file = './res/step_70_fig_221023/human_active_SNC_in_RG_NPC_DA_MRPA_SNC_overlap.pdf',width = 4,height = 4)
ggvenn(data = temp)
dev.off()

#seems not high


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
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  return(SNC_human[SNC_human@elementMetadata[,x] == 'YES'])
})
names(SNC_list) <- cell_type_list

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_70_fig_221023/human_active_SNC_cell_type_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

#different intersect size
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- as.data.frame(SNC_human@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
SNC_human$intersect <- temp

SNC_list <- base::lapply(X = 0:13,FUN = function(x){
  return(SNC_human[SNC_human$intersect == x])
})
names(SNC_list) <- 0:13

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_70_fig_221023/human_active_SNC_with_intersect_size_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

# active in all species ---------------------------------------------------
#load data
SNC_human <- readRDS(file = './res/step_70_fig_221023/SNC_human_active_in_all.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

## cell type enrichment ----------------------------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- data.frame(cell_type = cell_type_list,
                   SNC = unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
                     sum(SNC_human@elementMetadata[,x] == 'YES')
                   })))
temp$cell_type <- factor(temp$cell_type,levels = rev(temp$cell_type))

pdf(file = './res/step_70_fig_221023/all_active_SNC_cell_type_enrichment.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x = cell_type,y = SNC,fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8) + 
  scale_fill_manual(values = color_param$celltype[as.character(temp$cell_type)]) + 
  theme_bw() + 
  ylim(c(0,350)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),
        aspect.ratio = 2.5) + 
  xlab('Cell type') + ylab('all active SNC number') + 
  coord_flip()
dev.off()

## cell type enrichment intersect ------------------------------------------
#UPset plot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp <- names(SNC_human)[SNC_human@elementMetadata[,x] == 'YES']
  return(temp)
})

names(temp) <- cell_type_list

#plot
pdf(file = './res/step_70_fig_221023/all_active_SNC_cell_type_enrichment_intersect_upsetplot.pdf',width = 12,height = 7)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#overlap number
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- rtracklayer::as.data.frame(x = SNC_human)
temp <- temp[,cell_type_list_dot]
temp <- base::lapply(X = rownames(temp),FUN = function(x){
  return(sum(temp[x,] == 'YES'))
})
temp <- unlist(temp)

temp <- as.data.frame(table(temp))
colnames(temp) <- c('Intersect','percent')
temp$percent <- temp$percent/sum(temp$percent)*100
temp$Intersect <- factor(temp$Intersect,levels = temp$Intersect)

pdf(file = './res/step_70_fig_221023/all_active_SNC_number_with_intersect_size.pdf',width = 5,height = 3)
ggplot(data = temp,aes(x = Intersect,y = percent)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8,color = 'grey') + 
  theme_bw() + 
  scale_y_break(breaks = c(6,86)) + 
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

## all active SNC genome distribution ------------------------------------
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
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  return(SNC_human[SNC_human@elementMetadata[,x] == 'YES'])
})
names(SNC_list) <- cell_type_list

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_70_fig_221023/all_active_SNC_cell_type_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

#different intersect size
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- as.data.frame(SNC_human@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
SNC_human$intersect <- temp

SNC_list <- base::lapply(X = 0:13,FUN = function(x){
  return(SNC_human[SNC_human$intersect == x])
})
names(SNC_list) <- 0:13

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_70_fig_221023/all_active_SNC_with_intersect_size_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

# human active SNC posit at regulation region -----------------------------
#load data
SNC_human <- readRDS(file = './res/step_70_fig_221023/SNC_human_DA_in_human.rds')
SNC_all <- readRDS(file = './res/step_70_fig_221023/SNC_human_active_in_all.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

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
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  return(SNC_human[SNC_human@elementMetadata[,x] == 'YES'])
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

SNC_list <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  return(SNC_all[SNC_all@elementMetadata[,x] == 'YES'])
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

pdf(file = './res/step_70_fig_221023/human_active_posit_at_regulation_region.pdf',width = 3,height = 5)
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
temp <- as.data.frame(SNC_human@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
SNC_human$intersect <- temp
SNC_human <- SNC_human[SNC_human$intersect != 0,]

temp <- as.data.frame(SNC_all@elementMetadata[,cell_type_list_dot])
temp <- rowSums(temp == 'YES')
SNC_all$intersect <- temp
SNC_all <- SNC_all[SNC_all$intersect != 0,]

SNC_list <- list(`human active` = SNC_human,
                 `all active` = SNC_all)

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_70_fig_221023/merged_human_active_all_active_SNC_genome_distribution.pdf',width = 6,height = 3)
p + theme(aspect.ratio = 0.3)
dev.off()