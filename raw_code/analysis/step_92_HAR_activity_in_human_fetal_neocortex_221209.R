#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: HAR activity in human fetal neocortex                           ##
## Data: 2022.12.09                                                                ##
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

# pre-process HAR data ----------------------------------------------------
#load data
HAR_region <- read.csv(file = './res/step_92_fig_221209/hg38_HAR.csv')

#turn to GRanges
temp <- HAR_region[,c("Chr","start","stop")]
colnames(temp) <- c('chrom','start','end')
temp <- as(temp,'GRanges')

temp$HAR_id <- HAR_region$HAR_ID
names(temp) <- temp$HAR_id

#save data
saveRDS(object = temp,file = './res/step_92_fig_221209/hg38_HAR.rds')

# Percentage of HAR active in human fetal neocortex -----------------------
#load data
HAR_region <- readRDS(file = './res/step_92_fig_221209/hg38_HAR.rds')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]

cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

#for loop count activity
for (i in cell_type_list_dot) {
  #load peakset
  peakset <- list.files(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls')
  peakset <- peakset[grep(pattern = i,x = peakset,fixed = TRUE)]
  peakset <- paste('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls',peakset,sep = '/')
  peakset <- readRDS(file = peakset)
  
  #add meta data
  HAR_region@elementMetadata[,i] <- 'NO'
  
  temp <- countOverlaps(query = HAR_region,subject = peakset)
  temp <- names(temp)[temp > 0]
  HAR_region@elementMetadata[HAR_region$HAR_id %in% temp,i] <- 'YES'
}

#calculate Brain sample activity
temp <- HAR_region@elementMetadata[,cell_type_list_dot]
temp <- as.data.frame(temp)
temp <- unlist(base::lapply(X = 1:length(HAR_region),FUN = function(x){
  temp_activity <- temp[x,]
  temp_activity <- sum(temp_activity == 'YES')
  if(temp_activity > 0){
    return('YES')
  }else{
    return('NO')
  }
}))
HAR_region$Brain_activity <- temp

temp <- HAR_region@elementMetadata[,cell_type_list_dot]
temp <- as.data.frame(temp)
temp <- unlist(base::lapply(X = 1:length(HAR_region),FUN = function(x){
  temp_activity <- temp[x,]
  temp_activity <- sum(temp_activity == 'YES')
  return(temp_activity)
}))
HAR_region$Brain_activity_size <- temp

#save data
saveRDS(object = HAR_region,file = './res/step_92_fig_221209/HAR_neocortex_activity.rds')

# investigate HAR MPRA activity distribution -----------------------------------
#load data
HAR_region <- readRDS(file = './res/step_92_fig_221209/HAR_neocortex_activity.rds')
HAR_activity <- read.csv(file = './res/step_92_fig_221209/CaMPRA_activity.csv',row.names = 1)

HAR_region <- HAR_region[HAR_region$HAR_id %in% rownames(HAR_activity)]
HAR_activity <- HAR_activity[HAR_region$HAR_id,]
HAR_activity <- log2(HAR_activity)

#determining active HAR using MPRA
MPRA_activity <- unlist(base::lapply(X = rownames(HAR_activity),FUN = function(x){
  #SHSY5Y_Human
  SHSY5Y_Human <- HAR_activity[x,grep(pattern = '^SHSY5Y_Human',x = colnames(HAR_activity),fixed = FALSE)]
  SHSY5Y_Human <- t.test(x = SHSY5Y_Human,mu = 0,alternative = 'greater')
  SHSY5Y_Human <- SHSY5Y_Human$p.value < 0.05
  
  #SHSY5Y_Chimpanzee
  SHSY5Y_Chimpanzee <- HAR_activity[x,grep(pattern = '^SHSY5Y_Chimpanzee',x = colnames(HAR_activity),fixed = FALSE)]
  SHSY5Y_Chimpanzee <- t.test(x = SHSY5Y_Chimpanzee,mu = 0,alternative = 'greater')
  SHSY5Y_Chimpanzee <- SHSY5Y_Chimpanzee$p.value < 0.05
  
  #N2A_Human
  N2A_Human <- HAR_activity[x,grep(pattern = '^N2A_Human',x = colnames(HAR_activity),fixed = FALSE)]
  N2A_Human <- t.test(x = N2A_Human,mu = 0,alternative = 'greater')
  N2A_Human <- N2A_Human$p.value < 0.05
  
  #N2A_Chimpanzee
  N2A_Chimpanzee <- HAR_activity[x,grep(pattern = '^N2A_Chimpanzee',x = colnames(HAR_activity),fixed = FALSE)]
  N2A_Chimpanzee <- t.test(x = N2A_Chimpanzee,mu = 0,alternative = 'greater')
  N2A_Chimpanzee <- N2A_Chimpanzee$p.value < 0.05
  
  #return
  if(SHSY5Y_Human | SHSY5Y_Chimpanzee | N2A_Human | N2A_Chimpanzee){
    return('YES')
  }else{
    return('NO')
  }
}))

HAR_region$MPRA_all_activity <- MPRA_activity

#determining human specific active HAR using MPRA
MPRA_activity <- unlist(base::lapply(X = rownames(HAR_activity),FUN = function(x){
  #SHSY5Y_Human
  SHSY5Y_Human <- HAR_activity[x,grep(pattern = '^SHSY5Y_Human',x = colnames(HAR_activity),fixed = FALSE)]
  SHSY5Y_Human <- t.test(x = SHSY5Y_Human,mu = 0,alternative = 'greater')
  SHSY5Y_Human <- SHSY5Y_Human$p.value < 0.05
  
  #N2A_Human
  N2A_Human <- HAR_activity[x,grep(pattern = '^N2A_Human',x = colnames(HAR_activity),fixed = FALSE)]
  N2A_Human <- t.test(x = N2A_Human,mu = 0,alternative = 'greater')
  N2A_Human <- N2A_Human$p.value < 0.05
  
  #return
  if(SHSY5Y_Human | N2A_Human){
    return('YES')
  }else{
    return('NO')
  }
}))

HAR_region$MPRA_human_activity <- MPRA_activity

#save data
saveRDS(object = HAR_region,file = './res/step_92_fig_221209/processed_HAR_region.rds')

# investigate HAR activity distribution -----------------------------------
#load data
HAR_region <- readRDS(file = './res/step_92_fig_221209/processed_HAR_region.rds')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

#active percentage
temp <- data.frame(table(HAR_region$Brain_activity))
colnames(temp) <- c('activity','num')
temp$proportion <- temp$num/sum(temp$num)
temp$group <- 'ATAC'
HAR_activity <- temp

temp <- data.frame(table(HAR_region$MPRA_human_activity))
colnames(temp) <- c('activity','num')
temp$proportion <- temp$num/sum(temp$num)
temp$group <- 'MPRA'
HAR_activity <- rbind(HAR_activity,temp)

HAR_activity$group <- factor(HAR_activity$group,levels = c('ATAC','MPRA'))

pdf(file = './res/step_92_fig_221209/HAR_region_activity_proportion_barplot.pdf',width = 4.5,height = 1.5)
ggplot(data = HAR_activity,aes(x = group,y = proportion,fill = activity)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7) + 
  scale_fill_manual(values = pal_d3()(2)) + coord_flip() + 
  theme_bw() + theme(aspect.ratio = 0.3)
dev.off()

#overlap with MPRA
temp <- list(ATAC = HAR_region$HAR_id[HAR_region$Brain_activity == 'YES'],
             MPRA = HAR_region$HAR_id[HAR_region$MPRA_human_activity == 'YES'])

pdf(file = './res/step_92_fig_221209/HAR_region_activity_ATAC_MPRA_overlap_vennplot.pdf',width = 4.5,height = 4.5)
ggvenn(data = temp,columns = c('ATAC','MPRA'))
dev.off()

#cell type distribution
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- cell_type_list[!(cell_type_list %in% c('End','Per'))]
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp_list <- HAR_region@elementMetadata[,x]
  return(sum(temp_list == 'YES'))
}))
temp <- data.frame(cell_type = cell_type_list,num = temp)
temp$cell_type <- factor(temp$cell_type,levels = cell_type_list)

pdf(file = './res/step_92_fig_221209/active_HAR_number_in_each_cell_type_barplot.pdf',width = 5,height = 3)
ggplot(data = temp,aes(x = cell_type,y = num,fill = cell_type)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7) + 
  scale_fill_manual(values = color_param$celltype[cell_type_list]) + 
  coord_flip() + theme_bw() + theme(aspect.ratio = 0.5,legend.position = 'none') + 
  xlab('') + ylab('active HAR number')
dev.off()

#cell type active HAR overalp
temp <- base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp_list <- HAR_region$HAR_id[HAR_region@elementMetadata[,x] == 'YES']
  return(temp_list)
})
names(temp) <- cell_type_list

pdf(file = './res/step_92_fig_221209/cell_type_active_HAR_overlap_upset_plot.pdf',width = 10,height = 5)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#HAR genome distribution
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

temp <- list(HAR_region,HAR_region[HAR_region$Brain_activity == 'YES'],HAR_region[HAR_region$MPRA_human_activity == 'YES'])
names(temp) <- c('Global','ATAC','MPRA')
peakAnno <- base::lapply(X = temp,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_92_fig_221209/HAR_region_genome_distribution_barplot.pdf',width = 5,height = 5)
p + theme_bw() + theme(aspect.ratio = 0.3,legend.position = 'bottom')
dev.off()

pdf(file = './res/step_92_fig_221209/HAR_region_genome_distribution_barplot_legend.pdf',width = 10,height = 5)
p + theme_bw() + theme(aspect.ratio = 0.3,legend.position = 'bottom')
dev.off()

temp <- base::lapply(X = cell_type_list,FUN = function(x){
  cell_type <- x
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  temp_region <- HAR_region[HAR_region@elementMetadata[,cell_type_dot] == 'YES']
  return(temp_region)
})
names(temp) <- cell_type_list
peakAnno <- base::lapply(X = temp,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Hs.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_92_fig_221209/cell_type_HAR_region_genome_distribution_barplot_legend.pdf',width = 10,height = 3)
p + theme_bw() + theme(aspect.ratio = 0.5)
dev.off()