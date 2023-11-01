#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: compare DESeq2 and wilcoxon calculating DAP                     ##
## Data: 2022.06.14                                                                ##
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

# Ex-3 --------------------------------------------------------------------
cell_type <- 'Ex-3'

#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
human_peak_matrix <- readRDS(file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')
macaque_peak_matrix <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
mouse_peak_matrix <- readRDS(file = './res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

char <- sub(pattern = '-',replacement = '_',fixed = TRUE,x = cell_type)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
DESeq2_group <- readRDS(file = char)

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]
human_peak_matrix <- human_peak_matrix[,human_peak_matrix$cell_type == cell_type]
macaque_peak_matrix <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type]
mouse_peak_matrix <- mouse_peak_matrix[,mouse_peak_matrix$cell_type == cell_type]
gc()

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

#rename peakset
names(Brain_ATAC_peakset$macaque) <- paste(Brain_ATAC_peakset$macaque@seqnames,as.character(Brain_ATAC_peakset$macaque@ranges),sep = '-')
table(rownames(macaque_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$macaque))

names(Brain_ATAC_peakset$mouse) <- paste(Brain_ATAC_peakset$mouse@seqnames,as.character(Brain_ATAC_peakset$mouse@ranges),sep = '-')
table(rownames(mouse_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$mouse))

#human
temp <- human_peak_matrix@assays$RNA@counts[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/human_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
human_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
human_peak_matrix@assays$converted@data <- temp

#macaque
temp <- macaque_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$macaque[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/macaque_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
macaque_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
macaque_peak_matrix@assays$converted@data <- temp

#mouse
temp <- mouse_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$mouse[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/mouse_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
mouse_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
mouse_peak_matrix@assays$converted@data <- temp

gc()

#DF analysis using wilcoxon
HR_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = macaque_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

HM_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

RM_df <- my_DF_wilcox_test(mat1 = macaque_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

#rename peakset
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name

#human specific
res_HR <- HR_df[HR_df$log2FC > 1 & HR_df$fdr < 0.01,]
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]

human_specific_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]
res_RM <- RM_df[RM_df$log2FC > 1 & RM_df$fdr < 0.01,]
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(rownames(res_HM),rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(human_macaque_conserved_peak,rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]
res_HM <- HM_df[abs(HM_df$log2FC) <= 1,]
res_RM <- RM_df[abs(RM_df$log2FC) <= 1,]

species_conserved_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
species_conserved_peak <- dplyr::intersect(species_conserved_peak,rownames(res_RM))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_group.rds')
saveRDS(group,file = char)

#intersect with DESeq2 group

##human specific
temp <- list(wilcoxon = names(group)[group == 'human_specific'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_specific'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##human_macaque_conserved
temp <- list(wilcoxon = names(group)[group == 'human_macaque_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##species_conserved
temp <- list(wilcoxon = names(group)[group == 'species_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'species_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

#Reads in cell type peaks normalized lfc
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('human_specific','human_macaque_conserved','species_conserved'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
pdf(file = char,width = 8,height = 4)
ggplot(data = peak_lfc,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  facet_wrap(~group,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = 'peak matrix lfc normalized by Reads In Cell Type Peaks')
dev.off()

#enriched heatmap
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order)
dev.off()

#enrichment distribution
my_insertion_plot <- function(group=group,p1 = p1,p2 = p2,p3 = p3,group.by=group.by){
  temp <- p1@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'mouse'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5','mouse' = '#FDB31A')) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group.by) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_specific')
p_human_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

#compare lfc between wilcoxon with DESeq2
wilcoxon_group <- group
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

##human_specific
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_specific'],names(DESeq2_group)[DESeq2_group == 'human_specific']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human specific peak lfc')
dev.off()

##human_macaque_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_macaque_conserved'],names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human macaque conserved peak lfc')
dev.off()

##species_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'species_conserved'],names(DESeq2_group)[DESeq2_group == 'species_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'species conserved peak lfc')
dev.off()

# Ex-1 --------------------------------------------------------------------
cell_type <- 'Ex-1'

#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
human_peak_matrix <- readRDS(file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')
macaque_peak_matrix <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
mouse_peak_matrix <- readRDS(file = './res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

char <- sub(pattern = '-',replacement = '_',fixed = TRUE,x = cell_type)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
DESeq2_group <- readRDS(file = char)

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]
human_peak_matrix <- human_peak_matrix[,human_peak_matrix$cell_type == cell_type]
macaque_peak_matrix <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type]
mouse_peak_matrix <- mouse_peak_matrix[,mouse_peak_matrix$cell_type == cell_type]
gc()

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

#rename peakset
names(Brain_ATAC_peakset$macaque) <- paste(Brain_ATAC_peakset$macaque@seqnames,as.character(Brain_ATAC_peakset$macaque@ranges),sep = '-')
table(rownames(macaque_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$macaque))

names(Brain_ATAC_peakset$mouse) <- paste(Brain_ATAC_peakset$mouse@seqnames,as.character(Brain_ATAC_peakset$mouse@ranges),sep = '-')
table(rownames(mouse_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$mouse))

#human
temp <- human_peak_matrix@assays$RNA@counts[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/human_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
human_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
human_peak_matrix@assays$converted@data <- temp

#macaque
temp <- macaque_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$macaque[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/macaque_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
macaque_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
macaque_peak_matrix@assays$converted@data <- temp

#mouse
temp <- mouse_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$mouse[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/mouse_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
mouse_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
mouse_peak_matrix@assays$converted@data <- temp

gc()

#DF analysis using wilcoxon
HR_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = macaque_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 4,
                           future.globals.maxSize = 20*(1024^3))

HM_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 4,
                           future.globals.maxSize = 20*(1024^3))

RM_df <- my_DF_wilcox_test(mat1 = macaque_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 4,
                           future.globals.maxSize = 20*(1024^3))

#rename peakset
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name

#human specific
res_HR <- HR_df[HR_df$log2FC > 1 & HR_df$fdr < 0.01,]
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]

human_specific_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]
res_RM <- RM_df[RM_df$log2FC > 1 & RM_df$fdr < 0.01,]
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(rownames(res_HM),rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(human_macaque_conserved_peak,rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]
res_HM <- HM_df[abs(HM_df$log2FC) <= 1,]
res_RM <- RM_df[abs(RM_df$log2FC) <= 1,]

species_conserved_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
species_conserved_peak <- dplyr::intersect(species_conserved_peak,rownames(res_RM))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_group.rds')
saveRDS(group,file = char)

#intersect with DESeq2 group

##human specific
temp <- list(wilcoxon = names(group)[group == 'human_specific'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_specific'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##human_macaque_conserved
temp <- list(wilcoxon = names(group)[group == 'human_macaque_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##species_conserved
temp <- list(wilcoxon = names(group)[group == 'species_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'species_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

#Reads in cell type peaks normalized lfc
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('human_specific','human_macaque_conserved','species_conserved'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
pdf(file = char,width = 8,height = 4)
ggplot(data = peak_lfc,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  facet_wrap(~group,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = 'peak matrix lfc normalized by Reads In Cell Type Peaks')
dev.off()

#enriched heatmap
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order)
dev.off()

#enrichment distribution
my_insertion_plot <- function(group=group,p1 = p1,p2 = p2,p3 = p3,group.by=group.by){
  temp <- p1@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'mouse'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5','mouse' = '#FDB31A')) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group.by) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_specific')
p_human_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

#compare lfc between wilcoxon with DESeq2
wilcoxon_group <- group
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

##human_specific
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_specific'],names(DESeq2_group)[DESeq2_group == 'human_specific']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human specific peak lfc')
dev.off()

##human_macaque_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_macaque_conserved'],names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human macaque conserved peak lfc')
dev.off()

##species_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'species_conserved'],names(DESeq2_group)[DESeq2_group == 'species_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'species conserved peak lfc')
dev.off()

# IP --------------------------------------------------------------------
cell_type <- 'IP'

#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
human_peak_matrix <- readRDS(file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')
macaque_peak_matrix <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
mouse_peak_matrix <- readRDS(file = './res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

char <- sub(pattern = '-',replacement = '_',fixed = TRUE,x = cell_type)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
DESeq2_group <- readRDS(file = char)

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]
human_peak_matrix <- human_peak_matrix[,human_peak_matrix$cell_type == cell_type]
macaque_peak_matrix <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type]
mouse_peak_matrix <- mouse_peak_matrix[,mouse_peak_matrix$cell_type == cell_type]
gc()

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

#rename peakset
names(Brain_ATAC_peakset$macaque) <- paste(Brain_ATAC_peakset$macaque@seqnames,as.character(Brain_ATAC_peakset$macaque@ranges),sep = '-')
table(rownames(macaque_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$macaque))

names(Brain_ATAC_peakset$mouse) <- paste(Brain_ATAC_peakset$mouse@seqnames,as.character(Brain_ATAC_peakset$mouse@ranges),sep = '-')
table(rownames(mouse_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$mouse))

#human
temp <- human_peak_matrix@assays$RNA@counts[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/human_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
human_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
human_peak_matrix@assays$converted@data <- temp

#macaque
temp <- macaque_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$macaque[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/macaque_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
macaque_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
macaque_peak_matrix@assays$converted@data <- temp

#mouse
temp <- mouse_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$mouse[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/mouse_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
mouse_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
mouse_peak_matrix@assays$converted@data <- temp

gc()

#DF analysis using wilcoxon
HR_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = macaque_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

HM_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

RM_df <- my_DF_wilcox_test(mat1 = macaque_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

#rename peakset
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name

#human specific
res_HR <- HR_df[HR_df$log2FC > 1 & HR_df$fdr < 0.01,]
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]

human_specific_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]
res_RM <- RM_df[RM_df$log2FC > 1 & RM_df$fdr < 0.01,]
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(rownames(res_HM),rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(human_macaque_conserved_peak,rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]
res_HM <- HM_df[abs(HM_df$log2FC) <= 1,]
res_RM <- RM_df[abs(RM_df$log2FC) <= 1,]

species_conserved_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
species_conserved_peak <- dplyr::intersect(species_conserved_peak,rownames(res_RM))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_group.rds')
saveRDS(group,file = char)

#intersect with DESeq2 group

##human specific
temp <- list(wilcoxon = names(group)[group == 'human_specific'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_specific'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##human_macaque_conserved
temp <- list(wilcoxon = names(group)[group == 'human_macaque_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##species_conserved
temp <- list(wilcoxon = names(group)[group == 'species_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'species_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

#Reads in cell type peaks normalized lfc
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('human_specific','human_macaque_conserved','species_conserved'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
pdf(file = char,width = 8,height = 4)
ggplot(data = peak_lfc,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  facet_wrap(~group,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = 'peak matrix lfc normalized by Reads In Cell Type Peaks')
dev.off()

#enriched heatmap
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order)
dev.off()

#enrichment distribution
my_insertion_plot <- function(group=group,p1 = p1,p2 = p2,p3 = p3,group.by=group.by){
  temp <- p1@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'mouse'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5','mouse' = '#FDB31A')) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group.by) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_specific')
p_human_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

#compare lfc between wilcoxon with DESeq2
wilcoxon_group <- group
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

##human_specific
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_specific'],names(DESeq2_group)[DESeq2_group == 'human_specific']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human specific peak lfc')
dev.off()

##human_macaque_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_macaque_conserved'],names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human macaque conserved peak lfc')
dev.off()

##species_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'species_conserved'],names(DESeq2_group)[DESeq2_group == 'species_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'species conserved peak lfc')
dev.off()

# RG --------------------------------------------------------------------
cell_type <- 'RG'

#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
human_peak_matrix <- readRDS(file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')
macaque_peak_matrix <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
mouse_peak_matrix <- readRDS(file = './res/step_35_fig_220606/mouse_peak_matrix_Seurat.rds')

Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

char <- sub(pattern = '-',replacement = '_',fixed = TRUE,x = cell_type)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
DESeq2_group <- readRDS(file = char)

#subset peak matrix
subset_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == cell_type]
human_peak_matrix <- human_peak_matrix[,human_peak_matrix$cell_type == cell_type]
macaque_peak_matrix <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type]
mouse_peak_matrix <- mouse_peak_matrix[,mouse_peak_matrix$cell_type == cell_type]
gc()

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

#rename peakset
names(Brain_ATAC_peakset$macaque) <- paste(Brain_ATAC_peakset$macaque@seqnames,as.character(Brain_ATAC_peakset$macaque@ranges),sep = '-')
table(rownames(macaque_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$macaque))

names(Brain_ATAC_peakset$mouse) <- paste(Brain_ATAC_peakset$mouse@seqnames,as.character(Brain_ATAC_peakset$mouse@ranges),sep = '-')
table(rownames(mouse_peak_matrix@assays$RNA@counts) %in% names(Brain_ATAC_peakset$mouse))

#human
temp <- human_peak_matrix@assays$RNA@counts[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/human_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
human_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
human_peak_matrix@assays$converted@data <- temp

#macaque
temp <- macaque_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$macaque[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/macaque_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
macaque_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
macaque_peak_matrix@assays$converted@data <- temp

#mouse
temp <- mouse_peak_matrix@assays$RNA@counts
rownames(temp) <- Brain_ATAC_peakset$mouse[rownames(temp)]$name
temp <- temp[peak_list,]
temp_count <- temp
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(x){
  normed <- temp[,x]/mouse_peak_matrix@meta.data[x,"Reads_In_Cell_Type_Peaks"]
  normed <- normed*median(c(human_peak_matrix$Reads_In_Cell_Type_Peaks,macaque_peak_matrix$Reads_In_Cell_Type_Peaks,mouse_peak_matrix$Reads_In_Cell_Type_Peaks))
  return(normed)
}))
rownames(temp) <- rownames(temp_count)
colnames(temp) <- colnames(temp_count)
mouse_peak_matrix[['converted']] <- CreateAssayObject(counts = temp_count,min.cells = 0,min.features = 0)
mouse_peak_matrix@assays$converted@data <- temp

gc()

#DF analysis using wilcoxon
HR_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = macaque_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

HM_df <- my_DF_wilcox_test(mat1 = human_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

RM_df <- my_DF_wilcox_test(mat1 = macaque_peak_matrix@assays$converted@data,
                           mat2 = mouse_peak_matrix@assays$converted@data,
                           alternative = 'two.sided',paired = FALSE,workers = 6,
                           future.globals.maxSize = 20*(1024^3))

#rename peakset
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name

#human specific
res_HR <- HR_df[HR_df$log2FC > 1 & HR_df$fdr < 0.01,]
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]

human_specific_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- HM_df[HM_df$log2FC > 1 & HM_df$fdr < 0.01,]
res_RM <- RM_df[RM_df$log2FC > 1 & RM_df$fdr < 0.01,]
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(rownames(res_HM),rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(human_macaque_conserved_peak,rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HR <- HR_df[abs(HR_df$log2FC) <= 1,]
res_HM <- HM_df[abs(HM_df$log2FC) <= 1,]
res_RM <- RM_df[abs(RM_df$log2FC) <= 1,]

species_conserved_peak <- dplyr::intersect(rownames(res_HR),rownames(res_HM))
species_conserved_peak <- dplyr::intersect(species_conserved_peak,rownames(res_RM))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_group.rds')
saveRDS(group,file = char)

#intersect with DESeq2 group

##human specific
temp <- list(wilcoxon = names(group)[group == 'human_specific'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_specific'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##human_macaque_conserved
temp <- list(wilcoxon = names(group)[group == 'human_macaque_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

##species_conserved
temp <- list(wilcoxon = names(group)[group == 'species_conserved'],DESeq2 = names(DESeq2_group)[DESeq2_group == 'species_conserved'])
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_wilcoxon_DESeq2_vennplot.pdf')
pdf(file = char,width = 4.5,height = 4.5)
ggvenn(data = temp,c('wilcoxon','DESeq2'))
dev.off()

#Reads in cell type peaks normalized lfc
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('human_specific','human_macaque_conserved','species_conserved'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
pdf(file = char,width = 8,height = 4)
ggplot(data = peak_lfc,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  facet_wrap(~group,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = 'peak matrix lfc normalized by Reads In Cell Type Peaks')
dev.off()

#enriched heatmap
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/human',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/macaque',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2

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
signal_coverage <- list.files(path = './res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_35_fig_220606/Reads_In_Cell_Type_Peaks_coverage/mouse',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.0075,0.015),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,row_order = p1@row_order)
dev.off()

#enrichment distribution
my_insertion_plot <- function(group=group,p1 = p1,p2 = p2,p3 = p3,group.by=group.by){
  temp <- p1@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'mouse'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5','mouse' = '#FDB31A')) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group.by) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_specific')
p_human_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'human_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

#compare lfc between wilcoxon with DESeq2
wilcoxon_group <- group
temp_human <- Seurat::AggregateExpression(object = human_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_human <- (temp_human$converted/sum(human_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_macaque <- Seurat::AggregateExpression(object = macaque_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_macaque <- (temp_macaque$converted/sum(macaque_peak_matrix$Reads_In_Cell_Type_Peaks))
temp_mouse <- Seurat::AggregateExpression(object = mouse_peak_matrix,assays = 'converted',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp_mouse <- (temp_mouse$converted/sum(mouse_peak_matrix$Reads_In_Cell_Type_Peaks))

mat <- cbind(temp_human,temp_macaque[rownames(temp_human),],temp_mouse[rownames(temp_human),])
colnames(mat) <- c('human','macaque','mouse')

##human_specific
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_specific'],names(DESeq2_group)[DESeq2_group == 'human_specific']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_specific_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human specific peak lfc')
dev.off()

##human_macaque_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'human_macaque_conserved'],names(DESeq2_group)[DESeq2_group == 'human_macaque_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_human_macaque_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'human macaque conserved peak lfc')
dev.off()

##species_conserved
group <- unique(c(names(wilcoxon_group)[wilcoxon_group == 'species_conserved'],names(DESeq2_group)[DESeq2_group == 'species_conserved']))
temp <- base::lapply(X = group,FUN = function(x){
  if((x %in% names(wilcoxon_group)) & !(x %in% names(DESeq2_group))){
    return('wilcoxon_specific')
  }else if((x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('both')
  }else if(!(x %in% names(wilcoxon_group)) & (x %in% names(DESeq2_group))){
    return('DESeq2_specific')
  }else{
    return('error')
  }
})
temp <- unlist(temp)
peak_lfc <- data.frame(peak = group,group = temp)
rownames(peak_lfc) <- peak_lfc$peak

peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  macaque_exp <- mat[x,"macaque"]
  return(log2(macaque_exp/human_exp))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat[x,"human"]
  mouse_exp <- mat[x,"mouse"]
  return(log2(mouse_exp/human_exp))
}))

temp <- peak_lfc[,c("peak","group","mouse_lfc")]
colnames(temp) <- c("peak","group","lfc")
peak_lfc <- peak_lfc[,c("peak","group","macaque_lfc")]
colnames(peak_lfc) <- c("peak","group","lfc")
peak_lfc$species <- 'macaque'
temp$species <- 'mouse'
peak_lfc <- rbind(peak_lfc,temp)
peak_lfc$group <- factor(peak_lfc$group,levels = c('wilcoxon_specific','DESeq2_specific','both'))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_37_fig_220614/',char,'_species_conserved_peak_lfc_wilcoxon_vs_DESeq2_boxplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = peak_lfc,aes(x=group,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-10,2)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = 1,linetype = 'dashed',color = 'red') + 
  geom_hline(yintercept = -1,linetype = 'dashed',color = 'red') + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1)) + 
  labs(title = 'species conserved peak lfc')
dev.off()