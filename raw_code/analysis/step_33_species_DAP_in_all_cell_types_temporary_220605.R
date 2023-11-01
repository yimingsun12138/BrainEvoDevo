#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: species DAP in all cell types temporary                         ##
## Data: 2022.06.05                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#notice:
#This is a temporary script,
#Later we are going to adjust the peaking calling parameter and DAP calculation method.

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

# IP DAP ------------------------------------------------------------------
cell_type <- 'IP'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#human specific
res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[res_HR$log2FoldChange > 1 & res_HR$padj < 0.01,]

res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

human_specific_peak <- dplyr::intersect(x = rownames(res_HR),y = rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[res_RM$log2FoldChange > 1 & res_RM$padj < 0.01,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(x = human_macaque_conserved_peak,y = rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[abs(res_HM$log2FoldChange) <= 1,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[abs(res_RM$log2FoldChange) <= 1,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

species_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
species_conserved_peak <- dplyr::intersect(x = species_conserved_peak,y = rownames(res_HR))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
saveRDS(group,file = char)

#nFrags normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"nFrags"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_nFrags_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by nFrags')
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.005,0.01),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_CPM_center_besed_heatmap.pdf')
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
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
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
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

# RG DAP ------------------------------------------------------------------
cell_type <- 'RG'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#human specific
res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[res_HR$log2FoldChange > 1 & res_HR$padj < 0.01,]

res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

human_specific_peak <- dplyr::intersect(x = rownames(res_HR),y = rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[res_RM$log2FoldChange > 1 & res_RM$padj < 0.01,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(x = human_macaque_conserved_peak,y = rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[abs(res_HM$log2FoldChange) <= 1,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[abs(res_RM$log2FoldChange) <= 1,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

species_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
species_conserved_peak <- dplyr::intersect(x = species_conserved_peak,y = rownames(res_HR))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
saveRDS(group,file = char)

#nFrags normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"nFrags"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_nFrags_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by nFrags')
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.005,0.01),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 10)
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
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
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
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

# Ex-1 DAP ------------------------------------------------------------------
cell_type <- 'Ex-1'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#human specific
res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[res_HR$log2FoldChange > 1 & res_HR$padj < 0.01,]

res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

human_specific_peak <- dplyr::intersect(x = rownames(res_HR),y = rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[res_RM$log2FoldChange > 1 & res_RM$padj < 0.01,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(x = human_macaque_conserved_peak,y = rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[abs(res_HM$log2FoldChange) <= 1,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[abs(res_RM$log2FoldChange) <= 1,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

species_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
species_conserved_peak <- dplyr::intersect(x = species_conserved_peak,y = rownames(res_HR))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
saveRDS(group,file = char)

#nFrags normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"nFrags"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_nFrags_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by nFrags')
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.005,0.01),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_CPM_center_besed_heatmap.pdf')
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
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
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
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

# Ex-3 DAP ------------------------------------------------------------------
cell_type <- 'Ex-3'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#human specific
res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[res_HR$log2FoldChange > 1 & res_HR$padj < 0.01,]

res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

human_specific_peak <- dplyr::intersect(x = rownames(res_HR),y = rownames(res_HM))
human_specific_peak <- human_specific_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_specific_peak],subject = human_peak) > 0]

ii <- rep('human_specific',times = length(human_specific_peak))
names(ii) <- human_specific_peak
group <- ii

#human macaque conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[res_HM$log2FoldChange > 1 & res_HM$padj < 0.01,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[res_RM$log2FoldChange > 1 & res_RM$padj < 0.01,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

human_macaque_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
human_macaque_conserved_peak <- dplyr::intersect(x = human_macaque_conserved_peak,y = rownames(res_HR))
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[human_macaque_conserved_peak],subject = human_peak) > 0]
human_macaque_conserved_peak <- human_macaque_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[human_macaque_conserved_peak],subject = macaque_peak) > 0]

ii <- rep('human_macaque_conserved',times = length(human_macaque_conserved_peak))
names(ii) <- human_macaque_conserved_peak
group <- append(group,ii)

#species conserved
res_HM <- results(object = dds_HM,contrast = c('species','human','mouse'))
res_HM <- na.omit(res_HM)
res_HM <- res_HM[abs(res_HM$log2FoldChange) <= 1,]

res_RM <- results(object = dds_RM,contrast = c('species','macaque','mouse'))
res_RM <- na.omit(res_RM)
res_RM <- res_RM[abs(res_RM$log2FoldChange) <= 1,]

res_HR <- results(object = dds_HR,contrast = c('species','human','macaque'))
res_HR <- na.omit(res_HR)
res_HR <- res_HR[abs(res_HR$log2FoldChange) <= 1,]

species_conserved_peak <- dplyr::intersect(x = rownames(res_HM),y = rownames(res_RM))
species_conserved_peak <- dplyr::intersect(x = species_conserved_peak,y = rownames(res_HR))
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$human[species_conserved_peak],subject = human_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$macaque[species_conserved_peak],subject = macaque_peak) > 0]
species_conserved_peak <- species_conserved_peak[countOverlaps(query = Brain_ATAC_peakset$mouse[species_conserved_peak],subject = mouse_peak) > 0]

ii <- rep('species_conserved',times = length(species_conserved_peak))
names(ii) <- species_conserved_peak
group <- append(group,ii)

#save group
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
saveRDS(group,file = char)

#nFrags normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"nFrags"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_nFrags_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by nFrags')
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/human_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/human_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/macaque_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/macaque_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
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
signal_coverage <- list.files(path = './res/step_28_fig_220511/mouse_cell_type_coverage/')
signal_coverage <- signal_coverage[grep(pattern = paste0('^',char,'-'),x = signal_coverage,fixed = FALSE)]
signal_coverage <- paste('./res/step_28_fig_220511/mouse_cell_type_coverage',signal_coverage,sep = '/')
signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','human_macaque_conserved','species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.005,0.01),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_CPM_center_besed_heatmap.pdf')
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
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  temp <- p3@matrix[names(group)[group == group.by],]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-100:-1,1:100)
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
char <- paste0('./res/step_33_fig_220605/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

# addition ----------------------------------------------------------------
#notice that macaque peak strength is systematically weaker in human_macaque_conserved_peakset
#evaluate the size factor effect.

## IP ----------------------------------------------------------------------
cell_type <- 'IP'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#ReadsInPeaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"ReadsInPeaks"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_ReadsInPeaks_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by ReadsInPeaks')
dev.off()

#evaluate SF and nFrags
dds_HR$depth <- (dds_HR$nFrags*nrow(dds_HR@colData))/sum(dds_HR$nFrags)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs nFrags depth',sep = ' '))
dev.off()

#evaluate SF and ReadsInPeaks
dds_HR$Peak_depth <- (dds_HR$ReadsInPeaks*nrow(dds_HR@colData))/sum(dds_HR$ReadsInPeaks)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_Peak_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = Peak_depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs ReadsInPeaks depth',sep = ' '))
dev.off()

## RG ----------------------------------------------------------------------
cell_type <- 'RG'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#ReadsInPeaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"ReadsInPeaks"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_ReadsInPeaks_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by ReadsInPeaks')
dev.off()

#evaluate SF and nFrags
dds_HR$depth <- (dds_HR$nFrags*nrow(dds_HR@colData))/sum(dds_HR$nFrags)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs nFrags depth',sep = ' '))
dev.off()

#evaluate SF and ReadsInPeaks
dds_HR$Peak_depth <- (dds_HR$ReadsInPeaks*nrow(dds_HR@colData))/sum(dds_HR$ReadsInPeaks)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_Peak_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = Peak_depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs ReadsInPeaks depth',sep = ' '))
dev.off()

## Ex-1 ----------------------------------------------------------------------
cell_type <- 'Ex-1'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#ReadsInPeaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"ReadsInPeaks"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_ReadsInPeaks_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by ReadsInPeaks')
dev.off()

#evaluate SF and nFrags
dds_HR$depth <- (dds_HR$nFrags*nrow(dds_HR@colData))/sum(dds_HR$nFrags)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs nFrags depth',sep = ' '))
dev.off()

#evaluate SF and ReadsInPeaks
dds_HR$Peak_depth <- (dds_HR$ReadsInPeaks*nrow(dds_HR@colData))/sum(dds_HR$ReadsInPeaks)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_Peak_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = Peak_depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs ReadsInPeaks depth',sep = ' '))
dev.off()

## Ex-3 ----------------------------------------------------------------------
cell_type <- 'Ex-3'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220603/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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

#DESeq2
dds_HR <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','macaque')]
meta_data <- dds_HR@meta.data
dds_HR <- dds_HR@assays$RNA@counts[peak_list,]
dds_HR <- DESeqDataSetFromMatrix(countData = dds_HR,colData = meta_data,design = ~ species)
dds_HR <- DESeq(object = dds_HR)

dds_HM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('human','mouse')]
meta_data <- dds_HM@meta.data
dds_HM <- dds_HM@assays$RNA@counts[peak_list,]
dds_HM <- DESeqDataSetFromMatrix(countData = dds_HM,colData = meta_data,design = ~ species)
dds_HM <- DESeq(object = dds_HM)

dds_RM <- subset_peak_matrix[,subset_peak_matrix$species %in% c('macaque','mouse')]
meta_data <- dds_RM@meta.data
dds_RM <- dds_RM@assays$RNA@counts[peak_list,]
dds_RM <- DESeqDataSetFromMatrix(countData = dds_RM,colData = meta_data,design = ~ species)
dds_RM <- DESeq(object = dds_RM)

#ReadsInPeaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"ReadsInPeaks"])*(10^8)
}
mat@assays$RNA@data <- temp

peak_lfc <- data.frame(peak = names(group),group = as.character(group))
rownames(peak_lfc) <- peak_lfc$peak
peak_lfc$macaque_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  macaque_exp <- mat@assays$RNA@data[x,mat$species == 'macaque']
  return(log2(mean(macaque_exp)/mean(human_exp)))
}))
peak_lfc$mouse_lfc <- unlist(base::lapply(X = peak_lfc$peak,FUN = function(x){
  human_exp <- mat@assays$RNA@data[x,mat$species == 'human']
  mouse_exp <- mat@assays$RNA@data[x,mat$species == 'mouse']
  return(log2(mean(mouse_exp)/mean(human_exp)))
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
char <- paste0('./res/step_33_fig_220605/',char,'_peak_matrix_lfc_normalized_by_ReadsInPeaks_boxplot.pdf')
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
  labs(title = 'peak matrix lfc normalized by ReadsInPeaks')
dev.off()

#evaluate SF and nFrags
dds_HR$depth <- (dds_HR$nFrags*nrow(dds_HR@colData))/sum(dds_HR$nFrags)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs nFrags depth',sep = ' '))
dev.off()

#evaluate SF and ReadsInPeaks
dds_HR$Peak_depth <- (dds_HR$ReadsInPeaks*nrow(dds_HR@colData))/sum(dds_HR$ReadsInPeaks)

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_human_macaque_SF_vs_Peak_depth_dotplot.pdf')
pdf(file = char,width = 5,height = 5)
ggplot(data = as.data.frame(dds_HR@colData),aes(x = Peak_depth,y = sizeFactor,color = species)) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5')) + 
  geom_smooth(method = 'lm',size = 0.1,color = 'blue') + 
  geom_abline(intercept = 0,slope = 1,size = 0.1,color = 'red',linetype = 'dashed') + 
  theme_cowplot() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5)) + 
  labs(title = paste(cell_type,'size factor vs ReadsInPeaks depth',sep = ' '))
dev.off()