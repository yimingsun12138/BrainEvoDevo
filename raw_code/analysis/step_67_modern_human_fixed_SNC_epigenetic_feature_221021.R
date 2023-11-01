#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: modern human fixed SNC epigenetic feature                       ##
## Data: 2022.10.21                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
#load ArchR project
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
SNC_GRanges <- readRDS(file = './data/public/The_cis_regulatory_effects_of_modern_human_specific_variants/processed_data/SNC_GRanges.rds')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

# get human cell type peaks -----------------------------------------------

#which cell type enrich these SNC
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

for (i in cell_type_list) {
  file_list <- paste0('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls/',i,'-reproduciblePeaks.gr.rds')
  print(file.exists(file_list))
  temp <- paste(i,'peak',sep = '.')
  assign(x = temp,value = readRDS(file_list))
  temp <- countOverlaps(query = SNC_GRanges,subject = get(temp))
  SNC_GRanges@elementMetadata[,i] <- 'NO'
  SNC_GRanges@elementMetadata[temp > 0,i] <- 'YES'
}

temp <- data.frame(cell_type = cell_type_list,
                   SNC = unlist(base::lapply(X = cell_type_list,FUN = function(x){
                     sum(SNC_GRanges@elementMetadata[,x] == 'YES')
                   })))
temp$cell_type <- gsub(pattern = '.',replacement = '-',x = temp$cell_type,fixed = TRUE)
temp$cell_type <- factor(temp$cell_type,levels = rev(temp$cell_type))

pdf(file = './res/step_67_fig_221021/human_fixed_SNC_cell_type_enrichment.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x = cell_type,y = SNC,fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8) + 
  scale_fill_manual(values = color_param$celltype[as.character(temp$cell_type)]) + 
  theme_bw() + 
  ylim(c(0,600)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),
        aspect.ratio = 2.5) + 
  xlab('Cell type') + ylab('SNC number') + 
  coord_flip()
dev.off()

# human cell type level coverage ------------------------------------------
#human (run on denglab server)
temp <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
ArrowFiles <- getArrowFiles(ArchRProj = temp)
file.exists(ArrowFiles)
Greenleaf_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = './Greenleaf_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)
Greenleaf_ATAC_ArchR <- Greenleaf_ATAC_ArchR[rownames(temp@cellColData)]
Greenleaf_ATAC_ArchR@cellColData <- temp@cellColData
table(Greenleaf_ATAC_ArchR$cell_type)

Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = getPeakSet(ArchRProj = temp))
Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR)
coverage_file <- getGroupBW(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInPeaks',tileSize = 25,maxCells = NULL,verbose = TRUE)

# modern human SNC cell type enrichment intersect -------------------------
#UPset plot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- base::lapply(X = cell_type_list,FUN = function(x){
  temp <- names(SNC_GRanges)[SNC_GRanges@elementMetadata[,x] == 'YES']
  return(temp)
})

cell_type_list <- gsub(pattern = '.',replacement = '-',x = cell_type_list,fixed = TRUE)
names(temp) <- cell_type_list

#plot
pdf(file = './res/step_67_fig_221021/human_fixed_SNC_cell_type_enrichment_overlap.pdf',width = 12,height = 7)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#overlap number
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- rtracklayer::as.data.frame(x = SNC_GRanges)
temp <- temp[,cell_type_list]
temp <- base::lapply(X = rownames(temp),FUN = function(x){
  return(sum(temp[x,] == 'YES'))
})
temp <- unlist(temp)

temp <- as.data.frame(table(temp))
colnames(temp) <- c('Intersect','percent')
temp$percent <- temp$percent/sum(temp$percent)*100
temp$Intersect <- factor(temp$Intersect,levels = temp$Intersect)

pdf(file = './res/step_67_fig_221021/human_fixed_SNC_cell_type_active_number_barplot.pdf',width = 5,height = 3)
ggplot(data = temp,aes(x = Intersect,y = percent)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8,color = 'grey') + 
  theme_bw() + 
  scale_y_break(breaks = c(8,80)) + 
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

# overlap with NPC active -------------------------------------------------
temp <- list(NPC_MPRA = names(SNC_GRanges)[SNC_GRanges$NPC_active == 'Yes'],
             RG_ATAC = names(SNC_GRanges)[SNC_GRanges$RG.1 == 'YES'|SNC_GRanges$RG.2 == 'YES'])

pdf(file = './res/step_67_fig_221021/human_fixed_SNC_NPC_activity_RG_ATAC_overlap.pdf',width = 4,height = 4)
ggvenn(data = temp)
dev.off()

#seems not that high

# validate the accessibility by enrich heatmap ----------------------------

## RG.1 --------------------------------------------------------------------
cell_type <- 'RG.1'

#load coverage
signal_coverage <- list.files(path = './res/step_67_fig_221021/GroupBigWigs/cell_type/')
signal_coverage <- signal_coverage[grep(pattern = cell_type,x = signal_coverage,fixed = TRUE)]
signal_coverage <- paste('./res/step_67_fig_221021/GroupBigWigs/cell_type',signal_coverage,sep = '/')
file.exists(signal_coverage)

#normalize to matrix
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_GRanges,extend = 250,w = 25,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = FALSE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_GRanges@elementMetadata[,cell_type],levels = c('YES','NO')),
                      use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(5,'inches'))

col_fun <- colorRamp2(breaks = c(0,0.002,0.004),colors = c('#FFC30F','#C70039','#581845'))

char <- paste0('./res/step_67_fig_221021/human_fixed_SNC_',cell_type,'_enrichHeatmap.pdf')
pdf(file = char,width = 3.5,height = 6.5)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(SNC_GRanges[rownames(p1@matrix)]@elementMetadata[,cell_type],levels = c('YES','NO')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE,
                width = unit(1,'inches'),height = unit(5,'inches'))
dev.off()


## for loop ----------------------------------------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

for (i in cell_type_list) {
  #cell type
  cell_type <- i
  
  #load coverage
  signal_coverage <- list.files(path = './res/step_67_fig_221021/GroupBigWigs/cell_type/')
  signal_coverage <- signal_coverage[grep(pattern = cell_type,x = signal_coverage,fixed = TRUE)]
  signal_coverage <- paste('./res/step_67_fig_221021/GroupBigWigs/cell_type',signal_coverage,sep = '/')
  file.exists(signal_coverage)
  
  #normalize to matrix
  signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
  mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_GRanges,extend = 250,w = 25,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = FALSE,verbose = TRUE)
  
  p1 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_GRanges@elementMetadata[,cell_type],levels = c('YES','NO')),
                        use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(4,'inches'))
  
  col_fun <- colorRamp2(breaks = c(0,0.002,0.004),colors = c('#FFC30F','#C70039','#581845'))
  
  char <- paste0('./res/step_67_fig_221021/human_fixed_SNC_',cell_type,'_enrichHeatmap.pdf')
  pdf(file = char,width = 3.5,height = 6.5)
  print(EnrichedHeatmap(mat = p1@matrix,row_split = factor(SNC_GRanges[rownames(p1@matrix)]@elementMetadata[,cell_type],levels = c('YES','NO')),
                        use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE,
                        width = unit(1,'inches'),height = unit(5,'inches')))
  dev.off()
  print(paste(i,'done!',sep = ' '))
  gc()
}
my_send_sms('enrichHeatmap done!')

#save data
saveRDS(object = SNC_GRanges,file = './res/step_67_fig_221021/SNC_GRanges.rds')

# try intersect with +-250bp ----------------------------------------------
#modify SNC_GRanges
temp <- rtracklayer::as.data.frame(x = SNC_GRanges)
temp <- temp[,c("seqnames","start","end")]
temp$start <- temp$start-250
temp$end <- temp$end+250
colnames(temp) <- c('chrom','start','end')
temp <- as(temp,'GRanges')
names(temp) <- names(SNC_GRanges)
temp@elementMetadata <- SNC_GRanges@elementMetadata
temp

SNC_GRanges_ori <- SNC_GRanges
SNC_GRanges <- temp

## get human cell type peaks -----------------------------------------------

#which cell type enrich these SNC
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

for (i in cell_type_list) {
  file_list <- paste0('./processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/PeakCalls/',i,'-reproduciblePeaks.gr.rds')
  print(file.exists(file_list))
  temp <- paste(i,'peak',sep = '.')
  assign(x = temp,value = readRDS(file_list))
  temp <- countOverlaps(query = SNC_GRanges,subject = get(temp))
  SNC_GRanges@elementMetadata[,i] <- 'NO'
  SNC_GRanges@elementMetadata[temp > 0,i] <- 'YES'
}

temp <- data.frame(cell_type = cell_type_list,
                   SNC = unlist(base::lapply(X = cell_type_list,FUN = function(x){
                     sum(SNC_GRanges@elementMetadata[,x] == 'YES')
                   })))
temp$cell_type <- gsub(pattern = '.',replacement = '-',x = temp$cell_type,fixed = TRUE)
temp$cell_type <- factor(temp$cell_type,levels = rev(temp$cell_type))

pdf(file = './res/step_67_fig_221021/extent_plot/human_fixed_SNC_cell_type_enrichment.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x = cell_type,y = SNC,fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8) + 
  scale_fill_manual(values = color_param$celltype[as.character(temp$cell_type)]) + 
  theme_bw() + 
  ylim(c(0,1000)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),
        aspect.ratio = 2.5) + 
  xlab('Cell type') + ylab('SNC number') + 
  coord_flip()
dev.off()

## modern human SNC cell type enrichment intersect -------------------------
#UPset plot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- base::lapply(X = cell_type_list,FUN = function(x){
  temp <- names(SNC_GRanges)[SNC_GRanges@elementMetadata[,x] == 'YES']
  return(temp)
})

cell_type_list <- gsub(pattern = '.',replacement = '-',x = cell_type_list,fixed = TRUE)
names(temp) <- cell_type_list

pdf(file = './res/step_67_fig_221021/extent_plot/human_fixed_SNC_cell_type_enrichment_overlap.pdf',width = 12,height = 7)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#overlap number
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- rtracklayer::as.data.frame(x = SNC_GRanges)
temp <- temp[,cell_type_list]
temp <- base::lapply(X = rownames(temp),FUN = function(x){
  return(sum(temp[x,] == 'YES'))
})
temp <- unlist(temp)

temp <- as.data.frame(table(temp))
colnames(temp) <- c('Intersect','percent')
temp$percent <- temp$percent/sum(temp$percent)*100
temp$Intersect <- factor(temp$Intersect,levels = temp$Intersect)

pdf(file = './res/step_67_fig_221021/extent_plot/human_fixed_SNC_cell_type_active_number_barplot.pdf',width = 5,height = 3)
ggplot(data = temp,aes(x = Intersect,y = percent)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8,color = 'grey') + 
  theme_bw() + 
  scale_y_break(breaks = c(8,74)) + 
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

## overlap with NPC active -------------------------------------------------
temp <- list(NPC_MPRA = names(SNC_GRanges)[SNC_GRanges$NPC_active == 'Yes'],
             NPC_ATAC = names(SNC_GRanges)[SNC_GRanges$RG.1 == 'YES'|SNC_GRanges$RG.2 == 'YES'])

pdf(file = './res/step_67_fig_221021/extent_plot/human_fixed_SNC_NPC_activity_RG_ATAC_overlap.pdf',width = 4,height = 4)
ggvenn(data = temp)
dev.off()

## validate the accessibility by enrich heatmap ----------------------------

cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

for (i in cell_type_list) {
  #cell type
  cell_type <- i
  
  #load coverage
  signal_coverage <- list.files(path = './res/step_67_fig_221021/GroupBigWigs/cell_type/')
  signal_coverage <- signal_coverage[grep(pattern = cell_type,x = signal_coverage,fixed = TRUE)]
  signal_coverage <- paste('./res/step_67_fig_221021/GroupBigWigs/cell_type',signal_coverage,sep = '/')
  file.exists(signal_coverage)
  
  #normalize to matrix
  signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
  mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_GRanges_ori,extend = 250,w = 25,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = FALSE,verbose = TRUE)
  
  p1 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_GRanges@elementMetadata[,cell_type],levels = c('YES','NO')),
                        use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(4,'inches'))
  
  col_fun <- colorRamp2(breaks = c(0,0.002,0.004),colors = c('#FFC30F','#C70039','#581845'))
  
  char <- paste0('./res/step_67_fig_221021/extent_plot/human_fixed_SNC_',cell_type,'_enrichHeatmap.pdf')
  pdf(file = char,width = 3.5,height = 6.5)
  print(EnrichedHeatmap(mat = p1@matrix,row_split = factor(SNC_GRanges[rownames(p1@matrix)]@elementMetadata[,cell_type],levels = c('YES','NO')),
                        use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE,
                        width = unit(1,'inches'),height = unit(5,'inches')))
  dev.off()
  print(paste(i,'done!',sep = ' '))
  gc()
}
my_send_sms('enrichHeatmap done!')

#save data
saveRDS(object = SNC_GRanges,file = './res/step_67_fig_221021/extent_plot/SNC_GRanges.rds')
