#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: evaluate DESeq2 reliability in species ATAC peak comparision    ##
## Data: 2022.06.01                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# human and macaque, human and mouse --------------------------------------

## IPC ---------------------------------------------------------------------
cell_type <- 'IP'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

macaque_peak <- list.files(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
macaque_peak <- macaque_peak[grep(pattern = paste0('^',temp,'-'),x = macaque_peak,fixed = FALSE)]
macaque_peak <- readRDS(file = paste('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls',macaque_peak,sep = '/'))
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- countOverlaps(query = Brain_ATAC_peakset$macaque,subject = macaque_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#macaque lfc
macaque_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
target_site <- Brain_ATAC_peakset$macaque[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

macaque_coverage <- as.numeric(mat[,1])
names(macaque_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(macaque_coverage/human_coverage)
names(temp) <- peak_list
macaque_lfc$coverage_lfc <- temp[rownames(macaque_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_macaque_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = macaque_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

mouse_peak <- list.files(path = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
mouse_peak <- mouse_peak[grep(pattern = paste0('^',temp,'-'),x = mouse_peak,fixed = FALSE)]
mouse_peak <- readRDS(file = paste('./data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls',mouse_peak,sep = '/'))
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- countOverlaps(query = Brain_ATAC_peakset$mouse,subject = mouse_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#mouse lfc
mouse_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
target_site <- Brain_ATAC_peakset$mouse[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

mouse_coverage <- as.numeric(mat[,1])
names(mouse_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(mouse_coverage/human_coverage)
names(temp) <- peak_list
mouse_lfc$coverage_lfc <- temp[rownames(mouse_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_mouse_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = mouse_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#define DAP
human_up <- macaque_lfc[macaque_lfc$avg_log2FC < -1 & macaque_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- macaque_lfc[abs(macaque_lfc$avg_log2FC) <= 1 & 
                           countOverlaps(query = Brain_ATAC_peakset$macaque[rownames(macaque_lfc)],subject = macaque_peak) > 0 & 
                           countOverlaps(query = Brain_ATAC_peakset$human[rownames(macaque_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

macaque_up <- macaque_lfc[macaque_lfc$avg_log2FC > 1 & macaque_lfc$sig > 2,]
macaque_up <- rownames(macaque_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('macaque_up',times = length(macaque_up))
names(ii) <- macaque_up
group <- append(group,ii)

peak_lfc <- macaque_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_macaque_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- as.character(Brain_ATAC_peakset$macaque[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$macaque[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.01,0.02),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','macaque_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','macaque_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()

#define DAP
human_up <- mouse_lfc[mouse_lfc$avg_log2FC < -1 & mouse_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- mouse_lfc[abs(mouse_lfc$avg_log2FC) <= 1 & 
                         countOverlaps(query = Brain_ATAC_peakset$mouse[rownames(mouse_lfc)],subject = mouse_peak) > 0 & 
                         countOverlaps(query = Brain_ATAC_peakset$human[rownames(mouse_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

mouse_up <- mouse_lfc[mouse_lfc$avg_log2FC > 1 & mouse_lfc$sig > 2,]
mouse_up <- rownames(mouse_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('mouse_up',times = length(mouse_up))
names(ii) <- mouse_up
group <- append(group,ii)

peak_lfc <- mouse_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_mouse_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- as.character(Brain_ATAC_peakset$mouse[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$mouse[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','mouse_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','mouse_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()

## RG ---------------------------------------------------------------------
cell_type <- 'RG'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

macaque_peak <- list.files(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
macaque_peak <- macaque_peak[grep(pattern = paste0('^',temp,'-'),x = macaque_peak,fixed = FALSE)]
macaque_peak <- readRDS(file = paste('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls',macaque_peak,sep = '/'))
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- countOverlaps(query = Brain_ATAC_peakset$macaque,subject = macaque_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#macaque lfc
macaque_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
target_site <- Brain_ATAC_peakset$macaque[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

macaque_coverage <- as.numeric(mat[,1])
names(macaque_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(macaque_coverage/human_coverage)
names(temp) <- peak_list
macaque_lfc$coverage_lfc <- temp[rownames(macaque_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_macaque_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = macaque_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

mouse_peak <- list.files(path = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
mouse_peak <- mouse_peak[grep(pattern = paste0('^',temp,'-'),x = mouse_peak,fixed = FALSE)]
mouse_peak <- readRDS(file = paste('./data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls',mouse_peak,sep = '/'))
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- countOverlaps(query = Brain_ATAC_peakset$mouse,subject = mouse_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#mouse lfc
mouse_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
target_site <- Brain_ATAC_peakset$mouse[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

mouse_coverage <- as.numeric(mat[,1])
names(mouse_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(mouse_coverage/human_coverage)
names(temp) <- peak_list
mouse_lfc$coverage_lfc <- temp[rownames(mouse_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_mouse_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = mouse_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#define DAP
human_up <- macaque_lfc[macaque_lfc$avg_log2FC < -1 & macaque_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- macaque_lfc[abs(macaque_lfc$avg_log2FC) <= 1 & 
                           countOverlaps(query = Brain_ATAC_peakset$macaque[rownames(macaque_lfc)],subject = macaque_peak) > 0 & 
                           countOverlaps(query = Brain_ATAC_peakset$human[rownames(macaque_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

macaque_up <- macaque_lfc[macaque_lfc$avg_log2FC > 1 & macaque_lfc$sig > 2,]
macaque_up <- rownames(macaque_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('macaque_up',times = length(macaque_up))
names(ii) <- macaque_up
group <- append(group,ii)

peak_lfc <- macaque_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_macaque_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- as.character(Brain_ATAC_peakset$macaque[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$macaque[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','macaque_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','macaque_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()

#define DAP
human_up <- mouse_lfc[mouse_lfc$avg_log2FC < -1 & mouse_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- mouse_lfc[abs(mouse_lfc$avg_log2FC) <= 1 & 
                         countOverlaps(query = Brain_ATAC_peakset$mouse[rownames(mouse_lfc)],subject = mouse_peak) > 0 & 
                         countOverlaps(query = Brain_ATAC_peakset$human[rownames(mouse_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

mouse_up <- mouse_lfc[mouse_lfc$avg_log2FC > 1 & mouse_lfc$sig > 2,]
mouse_up <- rownames(mouse_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('mouse_up',times = length(mouse_up))
names(ii) <- mouse_up
group <- append(group,ii)

peak_lfc <- mouse_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',cell_type,'_mouse_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- as.character(Brain_ATAC_peakset$mouse[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$mouse[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','mouse_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','mouse_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()

## Ex-1 ---------------------------------------------------------------------
cell_type <- 'Ex-1'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

macaque_peak <- list.files(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
macaque_peak <- macaque_peak[grep(pattern = paste0('^',temp,'-'),x = macaque_peak,fixed = FALSE)]
macaque_peak <- readRDS(file = paste('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls',macaque_peak,sep = '/'))
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- countOverlaps(query = Brain_ATAC_peakset$macaque,subject = macaque_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#macaque lfc
macaque_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
target_site <- Brain_ATAC_peakset$macaque[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

macaque_coverage <- as.numeric(mat[,1])
names(macaque_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(macaque_coverage/human_coverage)
names(temp) <- peak_list
macaque_lfc$coverage_lfc <- temp[rownames(macaque_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = macaque_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

mouse_peak <- list.files(path = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
mouse_peak <- mouse_peak[grep(pattern = paste0('^',temp,'-'),x = mouse_peak,fixed = FALSE)]
mouse_peak <- readRDS(file = paste('./data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls',mouse_peak,sep = '/'))
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- countOverlaps(query = Brain_ATAC_peakset$mouse,subject = mouse_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#mouse lfc
mouse_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
target_site <- Brain_ATAC_peakset$mouse[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

mouse_coverage <- as.numeric(mat[,1])
names(mouse_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(mouse_coverage/human_coverage)
names(temp) <- peak_list
mouse_lfc$coverage_lfc <- temp[rownames(mouse_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = mouse_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#define DAP
human_up <- macaque_lfc[macaque_lfc$avg_log2FC < -1 & macaque_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- macaque_lfc[abs(macaque_lfc$avg_log2FC) <= 1 & 
                           countOverlaps(query = Brain_ATAC_peakset$macaque[rownames(macaque_lfc)],subject = macaque_peak) > 0 & 
                           countOverlaps(query = Brain_ATAC_peakset$human[rownames(macaque_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

macaque_up <- macaque_lfc[macaque_lfc$avg_log2FC > 1 & macaque_lfc$sig > 2,]
macaque_up <- rownames(macaque_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('macaque_up',times = length(macaque_up))
names(ii) <- macaque_up
group <- append(group,ii)

peak_lfc <- macaque_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- as.character(Brain_ATAC_peakset$macaque[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$macaque[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.01,0.02),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','macaque_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','macaque_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()

#define DAP
human_up <- mouse_lfc[mouse_lfc$avg_log2FC < -1 & mouse_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- mouse_lfc[abs(mouse_lfc$avg_log2FC) <= 1 & 
                         countOverlaps(query = Brain_ATAC_peakset$mouse[rownames(mouse_lfc)],subject = mouse_peak) > 0 & 
                         countOverlaps(query = Brain_ATAC_peakset$human[rownames(mouse_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

mouse_up <- mouse_lfc[mouse_lfc$avg_log2FC > 1 & mouse_lfc$sig > 2,]
mouse_up <- rownames(mouse_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('mouse_up',times = length(mouse_up))
names(ii) <- mouse_up
group <- append(group,ii)

peak_lfc <- mouse_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- as.character(Brain_ATAC_peakset$mouse[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$mouse[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.005,0.01),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','mouse_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','mouse_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()

## Ex-3 ---------------------------------------------------------------------
cell_type <- 'Ex-3'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

macaque_peak <- list.files(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
macaque_peak <- macaque_peak[grep(pattern = paste0('^',temp,'-'),x = macaque_peak,fixed = FALSE)]
macaque_peak <- readRDS(file = paste('./ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls',macaque_peak,sep = '/'))
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- countOverlaps(query = Brain_ATAC_peakset$macaque,subject = macaque_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#macaque lfc
macaque_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
target_site <- Brain_ATAC_peakset$macaque[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

macaque_coverage <- as.numeric(mat[,1])
names(macaque_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(macaque_coverage/human_coverage)
names(temp) <- peak_list
macaque_lfc$coverage_lfc <- temp[rownames(macaque_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = macaque_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#get peak set
human_peak <- list.files(path = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_peak <- human_peak[grep(pattern = paste0('^',temp,'-'),x = human_peak,fixed = FALSE)]
human_peak <- readRDS(file = paste('./data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls',human_peak,sep = '/'))
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

mouse_peak <- list.files(path = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/')
temp <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
mouse_peak <- mouse_peak[grep(pattern = paste0('^',temp,'-'),x = mouse_peak,fixed = FALSE)]
mouse_peak <- readRDS(file = paste('./data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls',mouse_peak,sep = '/'))
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- countOverlaps(query = Brain_ATAC_peakset$mouse,subject = mouse_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(subset_peak_matrix))

#mouse lfc
mouse_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#coverage matrix
#human
target_site <- Brain_ATAC_peakset$human[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

human_coverage <- as.numeric(mat[,1])
names(human_coverage) <- peak_list

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
target_site <- Brain_ATAC_peakset$mouse[peak_list]
names(target_site) <- peak_list

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 0,k = 1,limit = NA,
                         value_column = 'score',background = 0,target_ratio = 1,mean_mode = 'w0',
                         smooth = FALSE,verbose = TRUE)

mouse_coverage <- as.numeric(mat[,1])
names(mouse_coverage) <- peak_list

#compare DESeq2 lfc and coverage lfc
temp <- log2(mouse_coverage/human_coverage)
names(temp) <- peak_list
mouse_lfc$coverage_lfc <- temp[rownames(mouse_lfc)]

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_peak_list_lfc_DESeq2_vs_coverage_dotplot.pdf')
pdf(file = char,width = 6,height = 6)
ggplot(data = mouse_lfc,aes(x=avg_log2FC,y=coverage_lfc)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',color = 'red',formula = y ~ x,size = 1) + 
  geom_abline(slope = 1,intercept = 0,color = 'red',size = 1,linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1)
dev.off()

#define DAP
human_up <- macaque_lfc[macaque_lfc$avg_log2FC < -1 & macaque_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- macaque_lfc[abs(macaque_lfc$avg_log2FC) <= 1 & 
                           countOverlaps(query = Brain_ATAC_peakset$macaque[rownames(macaque_lfc)],subject = macaque_peak) > 0 & 
                           countOverlaps(query = Brain_ATAC_peakset$human[rownames(macaque_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

macaque_up <- macaque_lfc[macaque_lfc$avg_log2FC > 1 & macaque_lfc$sig > 2,]
macaque_up <- rownames(macaque_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('macaque_up',times = length(macaque_up))
names(ii) <- macaque_up
group <- append(group,ii)

peak_lfc <- macaque_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- as.character(Brain_ATAC_peakset$macaque[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$macaque[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','macaque_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.005,0.01),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_macaque_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','macaque_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','macaque_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()

#define DAP
human_up <- mouse_lfc[mouse_lfc$avg_log2FC < -1 & mouse_lfc$sig > 2,]
human_up <- rownames(human_up)

conserved <- mouse_lfc[abs(mouse_lfc$avg_log2FC) <= 1 & 
                         countOverlaps(query = Brain_ATAC_peakset$mouse[rownames(mouse_lfc)],subject = mouse_peak) > 0 & 
                         countOverlaps(query = Brain_ATAC_peakset$human[rownames(mouse_lfc)],subject = human_peak) > 0,]
conserved <- rownames(conserved)

mouse_up <- mouse_lfc[mouse_lfc$avg_log2FC > 1 & mouse_lfc$sig > 2,]
mouse_up <- rownames(mouse_up)

ii <- rep('human_up',times = length(human_up))
names(ii) <- human_up
group <- ii
ii <- rep('conserved',times = length(conserved))
names(ii) <- conserved
group <- append(group,ii)
ii <- rep('mouse_up',times = length(mouse_up))
names(ii) <- mouse_up
group <- append(group,ii)

peak_lfc <- mouse_lfc[names(group),]
peak_lfc$group <- group

p1 <- ggplot(data = peak_lfc,aes(x=group,y=avg_log2FC,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'DESeq2')

p2 <- ggplot(data = peak_lfc,aes(x=group,y=coverage_lfc,fill=group)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 14,face = 'bold',hjust = 0.5)) + 
  labs(title = 'coverage')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_grouped_peak_list_lfc_boxplot.pdf')
pdf(file = char,width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#coverage plot
#human
temp <- as.character(Brain_ATAC_peakset$human[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$human[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- as.character(Brain_ATAC_peakset$mouse[names(group)]@ranges)
temp <- base::lapply(temp,function(x){
  a <- as.numeric(strsplit(x = x,split = '-')[[1]][1])
  b <- as.numeric(strsplit(x = x,split = '-')[[1]][2])
  c <- round(abs(b+a)/2,digits = 0)
  return(c)
})
temp <- unlist(temp)
temp <- paste0(Brain_ATAC_peakset$mouse[names(group)]@seqnames,':',temp,'-',temp)
target_site <- as(temp,'GRanges')
names(target_site) <- names(group)

char <- sub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_up','mouse_up','conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

col_fun <- colorRamp2(breaks = c(0,0.005,0.01),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_34_fig_220601/',char,'_mouse_human_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_up','mouse_up','conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_up','mouse_up','conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  row_order = p1@row_order)
dev.off()