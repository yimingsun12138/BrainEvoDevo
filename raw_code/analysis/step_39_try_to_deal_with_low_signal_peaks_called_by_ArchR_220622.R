#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: try to deal with low signal peaks called by ArchR               ##
## Data: 2022.06.22                                                                ##
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
library(parallel)
library(Seurat)
library(ArchR)
library(rtracklayer)
library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(ComplexHeatmap)
library(parallel)
library(patchwork)
library(ggpubr)
library(ggvenn)
library(circlize)
library(ChIPpeakAnno)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnrichedHeatmap)
library(circlize)
library(scales)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(ggpointdensity)
library(networkD3)
library(htmlwidgets)
library(DESeq2)
library(Hmisc)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# regenerage macaque coverage file ----------------------------------------
#load peak matrix
macaque_peak_matrix <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')

#create ArchR project
temp
ArrowFiles
file.exists(ArrowFiles)
macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles,
  outputDirectory = './macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)
macaque_multiome_ArchR <- macaque_multiome_ArchR[colnames(macaque_peak_matrix)]

#add meta data
macaque_multiome_ArchR$nFrags_matrix <- macaque_peak_matrix$nFrags
macaque_multiome_ArchR$ReadsInTSS_matrix <- macaque_peak_matrix$ReadsInTSS
macaque_multiome_ArchR$cell_type <- macaque_peak_matrix$cell_type
macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks <- macaque_peak_matrix$Reads_In_Cell_Type_Peaks
macaque_multiome_ArchR$ReadsInPeaks <- macaque_peak_matrix$ReadsInPeaks

#export
table(macaque_multiome_ArchR$cell_type)
table(macaque_multiome_ArchR$nFrags_matrix == macaque_multiome_ArchR$nFrags)
table(macaque_multiome_ArchR$ReadsInTSS_matrix == macaque_multiome_ArchR$ReadsInTSS)

macaque_multiome_ArchR$ReadsInTSS <- macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
coverage_file <- getGroupBW(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',normMethod = 'ReadsInTSS',tileSize = 501,maxCells = NULL,verbose = TRUE)

# Ex-1 --------------------------------------------------------------------
cell_type <- 'Ex-1'

#load data
char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
coverage_file <- list.files('./res/step_39_fig_220622/macaque_TileSize-501-normMethod-Reads_In_Cell_Type_Peaks_coverage_file/')
coverage_file <- coverage_file[grep(pattern = paste0(char,'-'),x = coverage_file,fixed = TRUE)]
coverage_file <- paste('./res/step_39_fig_220622/macaque_TileSize-501-normMethod-Reads_In_Cell_Type_Peaks_coverage_file',coverage_file,sep = '/')
file.exists(coverage_file)
coverage_file <- rtracklayer::import.bw(con = coverage_file)
coverage_file$score <- (coverage_file$score / coverage_file@ranges@width) * 10^4

#add cell type peaks (run on cluster, ArchR is fucking shit)
macaque_multiome_ArchR <- loadArchRProject(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/')
macaque_multiome_ArchR <- macaque_multiome_ArchR[macaque_multiome_ArchR$cell_type == cell_type]

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
cell_type_peakset <- list.files('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/')
cell_type_peakset <- cell_type_peakset[grep(pattern = paste0(char,'-'),x = cell_type_peakset,fixed = TRUE)]
cell_type_peakset <- paste('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls',cell_type_peakset,sep = '/')
file.exists(cell_type_peakset)
cell_type_peakset <- readRDS(file = cell_type_peakset)

macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = cell_type_peakset,
                                     genomeAnnotation = getGenomeAnnotation(macaque_multiome_ArchR),force = TRUE)
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)
macaque_multiome_ArchR <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

meta_data <- as.data.frame(macaque_multiome_ArchR@colData)
ii <- macaque_multiome_ArchR@rowRanges
ii <- paste(ii@seqnames,as.character(ii@ranges),sep = '-')
macaque_multiome_ArchR <- macaque_multiome_ArchR@assays@data$PeakMatrix
rownames(macaque_multiome_ArchR) <- ii
macaque_multiome_ArchR <- (rowSums(macaque_multiome_ArchR)/sum(meta_data$ReadsInPeaks))*10^8

cell_type_peakset <- GRanges(seqnames = cell_type_peakset@seqnames,ranges = cell_type_peakset@ranges)
names(cell_type_peakset) <- paste(cell_type_peakset@seqnames,as.character(cell_type_peakset@ranges),sep = '-')
cell_type_peakset$score <- macaque_multiome_ArchR[names(cell_type_peakset)]
cell_type_peakset$score <- cell_type_peakset$score / cell_type_peakset@ranges@width

#plot cpm distribution
temp <- data.frame(tile = paste(coverage_file@seqnames,as.character(coverage_file@ranges),sep = '-'),
                   logcpm = log(coverage_file$score))
temp$overlap <- countOverlaps(query = coverage_file,subject = cell_type_peakset)
temp$group <- 'non-peak-tile'
temp[temp$overlap > 0,"group"] <- 'peak-tile'
temp <- rbind(temp,data.frame(tile = names(cell_type_peakset),logcpm = log(cell_type_peakset$score),overlap = 1,group = 'cell_type_peakset'))
cut_off <- quantile(temp$logcpm[temp$group == 'non-peak-tile'],probs = c(0.99))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_non_peak_tile_coverage_signal_density_plot.pdf')
pdf(file = char,width = 6,height = 4)
ggplot(data = temp,aes(x = logcpm,color = group)) + 
  geom_density() + 
  scale_color_manual(values = c('#272E6A','#208A42','#F47D2B')) + 
  geom_vline(xintercept = cut_off,color = 'grey',linetype = 'dashed') + 
  annotate(geom = 'text',label = round(x = cut_off,digits = 2),x = cut_off,y = 0,color = 'red',size = 6) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.6,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold')) + 
  labs(title = 'peak/non-peak tile coverage signal') + xlab('log(signal)') + ylab('density')
dev.off()

#quantile peak signal
cell_type_peakset$signal <- log(cell_type_peakset$score)
cell_type_peakset$group <- cut2(x = cell_type_peakset$signal,cuts = quantile(cell_type_peakset$signal,probs = seq(0,1,0.1),type = 5))

#generate enrichheatmap
char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
heatmap_coverage <- list.files('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque/')
heatmap_coverage <- heatmap_coverage[grep(pattern = paste0(char,'-'),x = heatmap_coverage,fixed = TRUE)]
heatmap_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque',heatmap_coverage,sep = '/')
file.exists(heatmap_coverage)
heatmap_coverage <- rtracklayer::import.bw(con = heatmap_coverage)
heatmap_coverage$score <- heatmap_coverage$score * 10^4

temp <- as.character(cell_type_peakset@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(cell_type_peakset@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(cell_type_peakset)

mat <- normalizeToMatrix(signal = heatmap_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

col_fun <- colorRamp2(breaks = c(0,100,200),colors = c('#f7fcfd','#8c96c6','#4d004b'))
p <- EnrichedHeatmap(mat = mat,pos_line = FALSE,row_split = factor(cell_type_peakset[rownames(mat)]$group,levels = rev(levels(cell_type_peakset$group))),
                     use_raster = TRUE,raster_resize_mat = mean,name = 'insertion',row_order = order(cell_type_peakset$score,decreasing = TRUE),
                     col = col_fun,width = unit(2,'inches'),height = unit(11,'inches'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_quantile_enriched_heatmap.pdf')
pdf(file = char,width = 4,height = 13)
p
dev.off()

#dot line plot
insertion_matrix <- base::do.call(what = rbind,args = base::lapply(X = levels(cell_type_peakset$group),FUN = function(x){
  temp <- names(cell_type_peakset)[cell_type_peakset$group == x]
  temp <- mat[temp,]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$group <- x
  return(temp)
}))
insertion_matrix$group <- factor(insertion_matrix$group,levels = levels(cell_type_peakset$group))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_quantile_enriched_dot_line_plot.pdf')
pdf(file = char,width = 10,height = 10)
ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=group)) + 
  geom_point(alpha = 0) + 
  geom_line() + 
  geom_hline(yintercept = 8,color = 'red',size = 0.5,linetype = 'dashed') + 
  geom_vline(xintercept = 5,color = 'red',size = 0.5,linetype = 'dashed') + 
  geom_vline(xintercept = -5,color = 'red',size = 0.5,linetype = 'dashed') + 
  scale_color_manual(values = as.character(ArchRPalettes$stallion)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA,colour = 'black'),
        axis.line = element_blank(),
        plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
  xlab('position') + ylab('normalized insertion')
dev.off()

#cut peakset
cell_type_peakset$peak <- 'YES'
cell_type_peakset$peak[cell_type_peakset$signal < cut_off] <- 'NO'

temp <- data.frame(tile = paste(coverage_file@seqnames,as.character(coverage_file@ranges),sep = '-'),
                   logcpm = log(coverage_file$score))
temp$overlap <- countOverlaps(query = coverage_file,subject = cell_type_peakset[cell_type_peakset$peak == 'YES'])
temp$group <- 'non-peak-tile'
temp[temp$overlap > 0,"group"] <- 'peak-tile'
temp <- rbind(temp,data.frame(tile = names(cell_type_peakset[cell_type_peakset$peak == 'YES']),logcpm = log(cell_type_peakset[cell_type_peakset$peak == 'YES']$score),overlap = 1,group = 'cell_type_peakset'))
cut_off <- quantile(temp$logcpm[temp$group == 'non-peak-tile'],probs = c(0.99))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_non_peak_tile_coverage_signal_density_plot_after_1_round_clean.pdf')
pdf(file = char,width = 6,height = 4)
ggplot(data = temp,aes(x = logcpm,color = group)) + 
  geom_density() + 
  scale_color_manual(values = c('#272E6A','#208A42','#F47D2B')) + 
  geom_vline(xintercept = cut_off,color = 'grey',linetype = 'dashed') + 
  annotate(geom = 'text',label = round(x = cut_off,digits = 2),x = cut_off,y = 0,color = 'red',size = 6) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.6,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold')) + 
  labs(title = 'peak/non-peak tile coverage signal') + xlab('log(signal)') + ylab('density')
dev.off()

# RG --------------------------------------------------------------------
cell_type <- 'RG'

#load data
char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
coverage_file <- list.files('./res/step_39_fig_220622/macaque_TileSize-501-normMethod-Reads_In_Cell_Type_Peaks_coverage_file/')
coverage_file <- coverage_file[grep(pattern = paste0(char,'-'),x = coverage_file,fixed = TRUE)]
coverage_file <- paste('./res/step_39_fig_220622/macaque_TileSize-501-normMethod-Reads_In_Cell_Type_Peaks_coverage_file',coverage_file,sep = '/')
file.exists(coverage_file)
coverage_file <- rtracklayer::import.bw(con = coverage_file)
coverage_file$score <- (coverage_file$score / coverage_file@ranges@width) * 10^4

#add cell type peaks (run on cluster, ArchR is fucking shit)
macaque_multiome_ArchR <- loadArchRProject(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/')
macaque_multiome_ArchR <- macaque_multiome_ArchR[macaque_multiome_ArchR$cell_type == cell_type]

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
cell_type_peakset <- list.files('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/')
cell_type_peakset <- cell_type_peakset[grep(pattern = paste0(char,'-'),x = cell_type_peakset,fixed = TRUE)]
cell_type_peakset <- paste('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls',cell_type_peakset,sep = '/')
file.exists(cell_type_peakset)
cell_type_peakset <- readRDS(file = cell_type_peakset)

macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = cell_type_peakset,
                                     genomeAnnotation = getGenomeAnnotation(macaque_multiome_ArchR),force = TRUE)
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR,verbose = TRUE,force = TRUE)
macaque_multiome_ArchR <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

meta_data <- as.data.frame(macaque_multiome_ArchR@colData)
ii <- macaque_multiome_ArchR@rowRanges
ii <- paste(ii@seqnames,as.character(ii@ranges),sep = '-')
macaque_multiome_ArchR <- macaque_multiome_ArchR@assays@data$PeakMatrix
rownames(macaque_multiome_ArchR) <- ii
macaque_multiome_ArchR <- (rowSums(macaque_multiome_ArchR)/sum(meta_data$ReadsInPeaks))*10^8

cell_type_peakset <- GRanges(seqnames = cell_type_peakset@seqnames,ranges = cell_type_peakset@ranges)
names(cell_type_peakset) <- paste(cell_type_peakset@seqnames,as.character(cell_type_peakset@ranges),sep = '-')
cell_type_peakset$score <- macaque_multiome_ArchR[names(cell_type_peakset)]
cell_type_peakset$score <- cell_type_peakset$score / cell_type_peakset@ranges@width

#plot cpm distribution
temp <- data.frame(tile = paste(coverage_file@seqnames,as.character(coverage_file@ranges),sep = '-'),
                   logcpm = log(coverage_file$score))
temp$overlap <- countOverlaps(query = coverage_file,subject = cell_type_peakset)
temp$group <- 'non-peak-tile'
temp[temp$overlap > 0,"group"] <- 'peak-tile'
temp <- rbind(temp,data.frame(tile = names(cell_type_peakset),logcpm = log(cell_type_peakset$score),overlap = 1,group = 'cell_type_peakset'))
cut_off <- quantile(temp$logcpm[temp$group == 'non-peak-tile'],probs = c(0.99))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_non_peak_tile_coverage_signal_density_plot.pdf')
pdf(file = char,width = 6,height = 4)
ggplot(data = temp,aes(x = logcpm,color = group)) + 
  geom_density() + 
  scale_color_manual(values = c('#272E6A','#208A42','#F47D2B')) + 
  geom_vline(xintercept = cut_off,color = 'grey',linetype = 'dashed') + 
  annotate(geom = 'text',label = round(x = cut_off,digits = 2),x = cut_off,y = 0,color = 'red',size = 6) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.6,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold')) + 
  labs(title = 'peak/non-peak tile coverage signal') + xlab('log(signal)') + ylab('density')
dev.off()

#quantile peak signal
cell_type_peakset$signal <- log(cell_type_peakset$score)
cell_type_peakset$group <- cut2(x = cell_type_peakset$signal,cuts = quantile(cell_type_peakset$signal,probs = seq(0,1,0.1),type = 5))

#generate enrichheatmap
char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
heatmap_coverage <- list.files('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque/')
heatmap_coverage <- heatmap_coverage[grep(pattern = paste0(char,'-'),x = heatmap_coverage,fixed = TRUE)]
heatmap_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque',heatmap_coverage,sep = '/')
file.exists(heatmap_coverage)
heatmap_coverage <- rtracklayer::import.bw(con = heatmap_coverage)
heatmap_coverage$score <- heatmap_coverage$score * 10^4

temp <- as.character(cell_type_peakset@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(cell_type_peakset@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(cell_type_peakset)

mat <- normalizeToMatrix(signal = heatmap_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

col_fun <- colorRamp2(breaks = c(0,100,200),colors = c('#f7fcfd','#8c96c6','#4d004b'))
p <- EnrichedHeatmap(mat = mat,pos_line = FALSE,row_split = factor(cell_type_peakset[rownames(mat)]$group,levels = rev(levels(cell_type_peakset$group))),
                     use_raster = TRUE,raster_resize_mat = mean,name = 'insertion',row_order = order(cell_type_peakset$score,decreasing = TRUE),
                     col = col_fun,width = unit(2,'inches'),height = unit(13,'inches'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_quantile_enriched_heatmap.pdf')
pdf(file = char,width = 4,height = 15)
p
dev.off()

#dot line plot
insertion_matrix <- base::do.call(what = rbind,args = base::lapply(X = levels(cell_type_peakset$group),FUN = function(x){
  temp <- names(cell_type_peakset)[cell_type_peakset$group == x]
  temp <- mat[temp,]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$group <- x
  return(temp)
}))
insertion_matrix$group <- factor(insertion_matrix$group,levels = levels(cell_type_peakset$group))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_quantile_enriched_dot_line_plot.pdf')
pdf(file = char,width = 10,height = 10)
ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=group)) + 
  geom_point(alpha = 0) + 
  geom_line() + 
  geom_hline(yintercept = 8,color = 'red',size = 0.5,linetype = 'dashed') + 
  geom_vline(xintercept = 5,color = 'red',size = 0.5,linetype = 'dashed') + 
  geom_vline(xintercept = -5,color = 'red',size = 0.5,linetype = 'dashed') + 
  scale_color_manual(values = as.character(ArchRPalettes$stallion)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        panel.background = element_rect(fill = NA,colour = 'black'),
        axis.line = element_blank(),
        plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
  xlab('position') + ylab('normalized insertion')
dev.off()

#cut peakset
cell_type_peakset$peak <- 'YES'
cell_type_peakset$peak[cell_type_peakset$signal < cut_off] <- 'NO'

temp <- data.frame(tile = paste(coverage_file@seqnames,as.character(coverage_file@ranges),sep = '-'),
                   logcpm = log(coverage_file$score))
temp$overlap <- countOverlaps(query = coverage_file,subject = cell_type_peakset[cell_type_peakset$peak == 'YES'])
temp$group <- 'non-peak-tile'
temp[temp$overlap > 0,"group"] <- 'peak-tile'
temp <- rbind(temp,data.frame(tile = names(cell_type_peakset[cell_type_peakset$peak == 'YES']),logcpm = log(cell_type_peakset[cell_type_peakset$peak == 'YES']$score),overlap = 1,group = 'cell_type_peakset'))
cut_off <- quantile(temp$logcpm[temp$group == 'non-peak-tile'],probs = c(0.99))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_39_fig_220622/',char,'_peak_non_peak_tile_coverage_signal_density_plot_after_1_round_clean.pdf')
pdf(file = char,width = 6,height = 4)
ggplot(data = temp,aes(x = logcpm,color = group)) + 
  geom_density() + 
  scale_color_manual(values = c('#272E6A','#208A42','#F47D2B')) + 
  geom_vline(xintercept = cut_off,color = 'grey',linetype = 'dashed') + 
  annotate(geom = 'text',label = round(x = cut_off,digits = 2),x = cut_off,y = 0,color = 'red',size = 6) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.6,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold')) + 
  labs(title = 'peak/non-peak tile coverage signal') + xlab('log(signal)') + ylab('density')
dev.off()