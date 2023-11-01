#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque multiome peak greenleaf ATAC peak comparision           ##
## Data: 2022.03.09                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

# macaque peakset liftover ------------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_peakset <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
names(macaque_peakset) <- 1:length(macaque_peakset)
temp <- rtracklayer::as.data.frame(macaque_peakset)
names(macaque_peakset) <- paste(temp$seqnames,temp$start,temp$end,sep = '-')
table(duplicated(names(macaque_peakset)))
macaque_peakset$peakname <- names(macaque_peakset)

#lift over
# rtracklayer::export.bed(object = macaque_peakset,con = './res/step_19_fig_220309/macaque_multiome_peakset.bed')
macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_19_fig_220309/macaque_multiome_peakset_lifted.bed')

table(duplicated(macaque_peakset_lifted$name))

temp <- data.frame(group=c('mapped','unmapped'),
                   num=c(sum(names(macaque_peakset) %in% macaque_peakset_lifted$name),sum(!(names(macaque_peakset) %in% macaque_peakset_lifted$name))))
temp$label <- paste0(temp$group,' ',round(temp$num/sum(temp$num)*100),'%')

pdf(file = './res/step_19_fig_220309/macaque_multiome_liftover_mapped_proportion_pieplot.pdf',width = 3,height = 3.5)
ggpie(data = temp,x='num',label = 'label',lab.pos = 'out',fill='group') + 
  theme(legend.position = 'none') + 
  labs(title = 'macaque multiome peakset liftover') + 
  theme(plot.title = element_text(size = 12,face = 'bold',hjust = 0.5))
dev.off()

summary(macaque_peakset_lifted@ranges@width)

temp <- data.frame(width=log(macaque_peakset_lifted@ranges@width))

pdf(file = './res/step_19_fig_220309/macaque_multiome_liftover_mapped_width_distribution_plot.pdf',width = 6,height = 4)
ggplot(data = temp,aes(x=width)) + geom_density() + xlab('log(width)') + 
  geom_vline(aes(xintercept = 6.125),color = 'red',size = 0.5) + 
  geom_vline(aes(xintercept = 6.312),color = 'red',size = 0.5)
dev.off()

#cut at 550bp
macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_19_fig_220309/macaque_multiome_peakset_lifted.bed')
macaque_peakset_lifted <- macaque_peakset_lifted[macaque_peakset_lifted@ranges@width < 550]

table(duplicated(macaque_peakset_lifted$name))

temp <- data.frame(group=c('mapped','unmapped'),
                   num=c(sum(names(macaque_peakset) %in% macaque_peakset_lifted$name),sum(!(names(macaque_peakset) %in% macaque_peakset_lifted$name))))
temp$label <- paste0(temp$group,' ',round(temp$num/sum(temp$num)*100),'%')

pdf(file = './res/step_19_fig_220309/macaque_multiome_liftover_after_filter_mapped_proportion_pieplot.pdf',width = 3,height = 3.5)
ggpie(data = temp,x='num',label = 'label',lab.pos = 'out',fill='group') + 
  theme(legend.position = 'none') + 
  labs(title = 'macaque multiome peakset liftover') + 
  theme(plot.title = element_text(size = 12,face = 'bold',hjust = 0.5))
dev.off()

pdf(file = './res/step_19_fig_220309/macaque_multiome_liftover_after_filter_width_histgram.pdf',width = 8,height = 5)
hist(macaque_peakset_lifted@ranges@width)
dev.off()

# overlap with Greenleaf peakset ------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_19_fig_220309/macaque_multiome_peakset_lifted.bed')
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220309/Save-ArchR-Project.rds')
human_peakset <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
human_peakset <- GRanges(seqnames = human_peakset@seqnames,ranges = human_peakset@ranges)
macaque_peakset_lifted <- GRanges(seqnames = macaque_peakset_lifted@seqnames,ranges = macaque_peakset_lifted@ranges)
temp <- append(human_peakset,macaque_peakset_lifted)
rtracklayer::export.bed(object = temp,con = './res/step_19_fig_220309/macaque_human_peakset.bed')

macaque_human_merged_peakset <- rtracklayer::import.bed(con = './res/step_19_fig_220309/macaque_human_merged_peakset_d0.bed')
macaque_human_merged_peakset <- macaque_human_merged_peakset[macaque_human_merged_peakset@seqnames %in% unique(human_peakset@seqnames)]

macaque_peakset_lifted <- rtracklayer::import.bed(con = './res/step_19_fig_220309/macaque_multiome_peakset_lifted.bed')
macaque_peakset_lifted <- macaque_peakset_lifted[macaque_peakset_lifted@ranges@width < 550]
macaque_peakset_lifted <- macaque_peakset_lifted[macaque_peakset_lifted@seqnames %in% unique(human_peakset@seqnames)]

table(countOverlaps(query = macaque_peakset_lifted,subject = macaque_human_merged_peakset))
temp <- list(macaque=unique(findOverlaps(query = macaque_peakset_lifted,subject = macaque_human_merged_peakset)@to),
             human=unique(findOverlaps(query = human_peakset,subject = macaque_human_merged_peakset)@to))

pdf(file = './res/step_19_fig_220309/macaque_human_peakset_overlap_vennplot.pdf',width = 6,height = 5)
ggvenn::ggvenn(data = temp,c('macaque','human'))
dev.off()

pdf(file = './res/step_19_fig_220309/macaque_human_merged_peakset_width_histgram.pdf',width = 8,height = 5)
hist(log(macaque_human_merged_peakset@ranges@width))
dev.off()

#try different d
d_matrix <- matrix(ncol = 2,nrow = 10)
rownames(d_matrix) <- c(0,10,50,80,100,200,400,600,800,1000)
colnames(d_matrix) <- c('d','overlap')
d_matrix <- as.data.frame(d_matrix)
d_matrix$d <- rownames(d_matrix)

for(i in d_matrix$d){
  char <- paste0('./res/step_19_fig_220309/macaque_human_merged_peakset_d',i,'.bed')
  temp <- rtracklayer::import.bed(con = char)
  temp <- list(macaque=unique(findOverlaps(query = macaque_peakset_lifted,subject = temp)@to),
               human=unique(findOverlaps(query = human_peakset,subject = temp)@to))
  temp <- length(dplyr::intersect(temp$human,temp$macaque))/length(dplyr::union(temp$human,temp$macaque))
  d_matrix[i,"overlap"] <- round(temp*100,digits = 1)
}
d_matrix$d <- as.numeric(d_matrix$d)

pdf(file = './res/step_19_fig_220309/macaque_human_peakset_overlap_vs_d_dotplot.pdf',width = 6,height = 3)
ggplot(data = d_matrix,aes(x=d,y=overlap)) + geom_point(size=1,color='black') + 
  geom_smooth(method='lm',se=FALSE,linetype=1,color='red',size = 0.5)
dev.off()

#conclusion use d = 0
macaque_human_merged_peakset <- rtracklayer::import.bed(con = './res/step_19_fig_220309/macaque_human_merged_peakset_d0.bed')
macaque_human_merged_peakset <- macaque_human_merged_peakset[macaque_human_merged_peakset@seqnames %in% unique(human_peakset@seqnames)]

#marker peak overlap
macaque_markersPeaks <- getMarkerFeatures(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment","log10(nFrags)"),
  testMethod = "wilcoxon"
)

macaque_markersPeaks <- getMarkers(macaque_markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE)

human_markersPeaks <- getMarkerFeatures(
  ArchRProj = Greenleaf_ATAC_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment","log10(nFrags)"),
  testMethod = "wilcoxon"
)

human_markersPeaks <- getMarkers(human_markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE)

marker_peak_overlap <- matrix(nrow = length(macaque_markersPeaks),ncol = length(human_markersPeaks))
rownames(marker_peak_overlap) <- names(macaque_markersPeaks)
colnames(marker_peak_overlap) <- names(human_markersPeaks)
marker_peak_overlap <- as.data.frame(marker_peak_overlap)

for (i in rownames(marker_peak_overlap)) {
  for (j in names(marker_peak_overlap)) {
    A <- macaque_markersPeaks[[i]]
    A <- rtracklayer::as.data.frame(A)
    B <- human_markersPeaks[[j]]
    
    A <- paste(A$seqnames,A$start,A$end,sep = '-')
    A <- A[A %in% macaque_peakset_lifted$name]
    A <- macaque_peakset_lifted[macaque_peakset_lifted$name %in% A]
    
    A <- GenomicRanges::findOverlaps(query = A,subject = macaque_human_merged_peakset)@to
    B <- GenomicRanges::findOverlaps(query = B,subject = macaque_human_merged_peakset)@to
    marker_peak_overlap[i,j] <- length(dplyr::intersect(A,B))/length(dplyr::union(A,B))
  }
}

marker_peak_overlap <- marker_peak_overlap[c('InMGE','InCGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per','Ependymal'),
                                           c('InMGE','InCGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG','OPC','Mic','End','Per')]

pdf(file = './res/step_19_fig_220309/macaque_human_marker_peak_overlap_by_cell_type_heatmap.pdf',width = 7,height = 7)
Heatmap(marker_peak_overlap,cluster_rows = FALSE,cluster_columns = FALSE,
        row_split = factor(c('In','In','Ex','Ex','Ex','Ex','NP','NP','NP','NP','NN','NN','NN','NN'),
                           levels = c('In','Ex','NP','NN')),
        column_split = factor(c('In','In','Ex','Ex','Ex','Ex','NP','NP','NN','NN','NN','NN'),
                              levels = c('In','Ex','NP','NN')),border = TRUE,
        height = unit(5.25,"inches"),width = unit(4.5,"inches"),name = 'Jaccard idx',
        row_title = 'macaque multiome marker peak',column_title = 'Greenleaf ATAC marker peak')
dev.off()


# create pseudo bulk peak matrix ------------------------------------------
temp <- rtracklayer::as.data.frame(macaque_human_merged_peakset)
macaque_human_merged_peakset$peakname <- paste(temp$seqnames,temp$start,temp$end,sep = '-')

temp <- rtracklayer::import.chain(con = './data/reference/hg38ToRheMac10.over.chain')
macaque_peakset <- rtracklayer::liftOver(x = macaque_human_merged_peakset,chain = temp)
macaque_peakset <- unlist(macaque_peakset)
macaque_peakset <- macaque_peakset[macaque_peakset@seqnames %in% unique(getPeakSet(macaque_multiome_ArchR)@seqnames)]
temp <- macaque_peakset
temp <- rtracklayer::as.data.frame(temp)
temp <- paste(temp$seqnames,temp$start,temp$end,sep = '-')
macaque_peakset$new_peakname <- temp
macaque_peakset <- macaque_peakset[!duplicated(macaque_peakset$new_peakname)]

saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "/data/User/sunym/temp/macaque_multiome_ArchR", load = TRUE, overwrite = TRUE)
macaque_multiome_ArchR <- readRDS(file = '/data/User/sunym/temp/macaque_multiome_ArchR/Save-ArchR-Project.rds')

#macaque
macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = macaque_peakset,force = TRUE)
macaque_multiome_ArchR <- addPeakMatrix(macaque_multiome_ArchR)
macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix')
temp <- macaque_peak_matrix@rowRanges
temp <- rtracklayer::as.data.frame(temp)
temp <- paste(temp$seqnames,temp$start,temp$end,sep = '-')
macaque_peak_matrix <- macaque_peak_matrix@assays@data@listData$PeakMatrix
rownames(macaque_peak_matrix) <- temp
macaque_peak_matrix <- CreateSeuratObject(counts = macaque_peak_matrix,project = 'temp',assay = 'RNA',meta.data = as.data.frame(macaque_multiome_ArchR@cellColData),min.cells = 0,min.features = 0)

