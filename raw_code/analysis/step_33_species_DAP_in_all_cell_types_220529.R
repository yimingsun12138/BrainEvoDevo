#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: species DAP in all cell types                                   ##
## Data: 2022.05.29                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#notice:
#the coverage file is wrong, while doing IP enrichedHeatmap.
#now the code has been restored but the res are still out of date.

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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# RG DAP ------------------------------------------------------------------
cell_type <- 'RG'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')

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

#macaque lfc
macaque_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#mouse lfc
mouse_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#peak lfc matrix
peak_lfc_matrix <- data.frame(peak = peak_list,macaque_lfc = 0,macaque_sig = 0,mouse_lfc = 0,mouse_sig = 0)
rownames(peak_lfc_matrix) <- peak_lfc_matrix$peak
peak_lfc_matrix[rownames(macaque_lfc),"macaque_lfc"] <- macaque_lfc$avg_log2FC
peak_lfc_matrix[rownames(macaque_lfc),"macaque_sig"] <- macaque_lfc$sig
peak_lfc_matrix[rownames(mouse_lfc),"mouse_lfc"] <- mouse_lfc$avg_log2FC
peak_lfc_matrix[rownames(mouse_lfc),"mouse_sig"] <- mouse_lfc$sig

#try lfc=1, sig=2
peak_lfc_matrix$group <- 'mixed'

temp <- rownames(macaque_lfc)[abs(macaque_lfc$avg_log2FC) > 1 & macaque_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'
temp <- rownames(mouse_lfc)[abs(mouse_lfc$avg_log2FC) > 1 & mouse_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'

peak_lfc_matrix[abs(peak_lfc_matrix$macaque_lfc) <= 1 & 
                  abs(peak_lfc_matrix$mouse_lfc) <= 1 & 
                  countOverlaps(query = Brain_ATAC_peakset$macaque[peak_lfc_matrix$peak],subject = macaque_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$mouse[peak_lfc_matrix$peak],subject = mouse_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$human[peak_lfc_matrix$peak],subject = human_peak) > 0,
                "group"] <- 'conserved'

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peak_lfc_matrix.rds')
saveRDS(object = peak_lfc_matrix,file = char)

#kmeans
mat <- subset_peak_matrix
mat <- NormalizeData(object = mat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
VariableFeatures(mat) <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'differential']
mat <- ScaleData(object = mat,features = VariableFeatures(mat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)

group <- kmeans(mat@assays$RNA@scale.data[VariableFeatures(mat),],centers = 6)$cluster
col_anno <- HeatmapAnnotation(species = mat$species,which = 'column')

p1 <- Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
              show_column_names = FALSE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(mat@assays$RNA@scale.data)],
                                                      levels = unique(group)),
              cluster_columns = FALSE,
              column_split = factor(mat$species,levels = c('human','macaque','mouse')),
              height = unit(10,'inches'),width = unit(6,'inches'),
              name = 'z-score',top_annotation = col_anno,border = TRUE)

temp <- peak_lfc_matrix[rownames(p1@matrix)[p1@row_order],c("macaque_lfc","mouse_lfc")]
p2 <- Heatmap(matrix = as.matrix(temp),
              show_column_names = TRUE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(temp)],levels = unique(group)),
              cluster_columns = FALSE,
              height = unit(10,'inches'),width = unit(1,'inches'),
              name = 'lfc',border = TRUE)

p1+p2

#modify groups
temp <- peak_lfc_matrix[peak_lfc_matrix$group == 'conserved',]
ii <- rep('conserved',times = nrow(temp))
names(ii) <- temp$peak
group <- append(group,ii)
group <- data.frame(peak_name = names(group),raw_cluster = as.character(group))

group$cluster <- NA
group[group$raw_cluster == '6',"cluster"] <- 'macaque_human_conserved'
group[group$raw_cluster == '3',"cluster"] <- 'macaque_specific'
group[group$raw_cluster == '4',"cluster"] <- 'mouse_macaque_conserved'
group[group$raw_cluster == '2',"cluster"] <- 'mouse_specific'
group[group$raw_cluster == '5',"cluster"] <- 'mouse_human_conserved'
group[group$raw_cluster == '1',"cluster"] <- 'human_specific'
group[group$raw_cluster == 'conserved',"cluster"] <- 'species_conserved'

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peak_kmeans_group.rds')
saveRDS(group,file = char)

#plot
mat <- subset_peak_matrix
mat <- NormalizeData(object = mat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
VariableFeatures(mat) <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'differential']
mat <- ScaleData(object = mat,features = VariableFeatures(mat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)

group <- readRDS(file = char)
ii <- group$peak_name
group <- group$cluster
names(group) <- ii
col_anno <- HeatmapAnnotation(species = mat$species,which = 'column')

p1 <- Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
              show_column_names = FALSE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(mat@assays$RNA@scale.data)],
                                                      levels = c('human_specific','macaque_human_conserved',
                                                                 'mouse_human_conserved','macaque_specific',
                                                                 'mouse_macaque_conserved','mouse_specific',
                                                                 'species_conserved')),
              cluster_columns = FALSE,
              column_split = factor(mat$species,levels = c('human','macaque','mouse')),
              height = unit(10,'inches'),width = unit(6,'inches'),
              name = 'z-score',top_annotation = col_anno,border = TRUE,
              row_title = 'merged peakset')

temp <- peak_lfc_matrix[rownames(p1@matrix)[p1@row_order],c("macaque_lfc","mouse_lfc")]
p2 <- Heatmap(matrix = as.matrix(temp),
              show_column_names = TRUE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(temp)],levels = c('human_specific','macaque_human_conserved',
                                                                                       'mouse_human_conserved','macaque_specific',
                                                                                       'mouse_macaque_conserved','mouse_specific',
                                                                                       'species_conserved')),
              cluster_columns = FALSE,
              height = unit(10,'inches'),width = unit(1,'inches'),
              name = 'lfc',border = TRUE,row_title = 'merged peakset')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_differential_peakset_z_score_heatmap.pdf')
pdf(file = char,width = 10,height = 12)
p1+p2
dev.off()

#boxplot show lfc
peak_lfc_matrix <- peak_lfc_matrix[peak_lfc_matrix$group == 'differential',]
peak_lfc_matrix$cluster <- group[rownames(peak_lfc_matrix)]
temp <- data.frame(lfc = peak_lfc_matrix$macaque_lfc,species = 'macaque',cluster = peak_lfc_matrix$cluster)
temp <- rbind(temp,data.frame(lfc = peak_lfc_matrix$mouse_lfc,species = 'mouse',cluster = peak_lfc_matrix$cluster))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_differential_peakset_lfc_boxplot.pdf')
pdf(file = char,width = 7,height = 5)
ggplot(data = temp,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-5,5)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  facet_wrap(~cluster,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1))
dev.off()

# enrichheatmap
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

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
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

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
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

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.004,0.008),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                         'mouse_human_conserved','macaque_specific',
                                                                                         'mouse_macaque_conserved','mouse_specific',
                                                                                         'species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE)
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
p_macaque_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_specific')
p_mouse_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_specific')
p_macaque_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_human_conserved')
p_mouse_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_human_conserved')
p_mouse_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 18,height = 7)
p_human_specific+p_macaque_specific+p_mouse_specific+p_macaque_human_conserved+
  p_mouse_human_conserved+p_mouse_macaque_conserved+p_species_conserved + plot_layout(ncol = 4)
dev.off()

#peak contribution
#human
human_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$human[peak_list],subject = human_peak)
human_contributed_peak <- human_contributed_peak[human_contributed_peak > 0]
human_contributed_peak <- names(human_contributed_peak)

#macaque
macaque_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$macaque[peak_list],subject = macaque_peak)
macaque_contributed_peak <- macaque_contributed_peak[macaque_contributed_peak > 0]
macaque_contributed_peak <- names(macaque_contributed_peak)

#mouse
mouse_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$mouse[peak_list],subject = mouse_peak)
mouse_contributed_peak <- mouse_contributed_peak[mouse_contributed_peak > 0]
mouse_contributed_peak <- names(mouse_contributed_peak)

#ggvenn
temp <- list(human = human_contributed_peak,macaque = macaque_contributed_peak,mouse = mouse_contributed_peak)
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/species_contribution_to_final_',char,'_peakset.pdf')
pdf(file = char,width = 6.5,height = 6.5)
ggvenn(data = temp,c('human','macaque','mouse'))
dev.off()

#sankey plot
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peak_lfc_matrix.rds')
peak_lfc_matrix <- readRDS(file = char)
peak_lfc_matrix$kmeans_cluster <- 'dropped'
peak_lfc_matrix[names(group),"kmeans_cluster"] <- as.character(group)

temp <- base::lapply(X = peak_lfc_matrix$peak,FUN = function(x){
  contribution_score <- 1*sum(x %in% mouse_contributed_peak) + 10*sum(x %in% macaque_contributed_peak) + 100*sum(x %in% human_contributed_peak)
  if(contribution_score == 1){
    return('mouse_specific')
  }else if(contribution_score == 10){
    return('macaque_specific')
  }else if(contribution_score == 100){
    return('human_specific')
  }else if(contribution_score == 11){
    return('mouse_macaque_conserved')
  }else if(contribution_score == 101){
    return('mouse_human_conserved')
  }else if(contribution_score == 110){
    return('macaque_human_conserved')
  }else if(contribution_score == 111){
    return('species_conserved')
  }else{
    return('error!')
  }
})
temp <- unlist(temp)
peak_lfc_matrix$peak_cluster <- temp

#scibet
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_kmeans_cluster_and_peak_only_cluster_confusion_heatmap.pdf')
pdf(file = char,width = 6,height = 6)
scibet::Confusion_heatmap(ori = peak_lfc_matrix$peak_cluster,prd = peak_lfc_matrix$kmeans_cluster) + 
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 16,face = 'bold')) + 
  ylab('kmeans cluster') + xlab('peak only cluster')
dev.off()

# IP DAP ------------------------------------------------------------------
cell_type <- 'IP'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_32_fig_220528/Brain_peak_matrix_Seurat.rds')
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220528/Brain_ATAC_peakset.rds')

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

#macaque lfc
macaque_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

#mouse lfc
mouse_lfc <- FindMarkers(object = subset_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

#peak lfc matrix
peak_lfc_matrix <- data.frame(peak = peak_list,macaque_lfc = 0,macaque_sig = 0,mouse_lfc = 0,mouse_sig = 0)
rownames(peak_lfc_matrix) <- peak_lfc_matrix$peak
peak_lfc_matrix[rownames(macaque_lfc),"macaque_lfc"] <- macaque_lfc$avg_log2FC
peak_lfc_matrix[rownames(macaque_lfc),"macaque_sig"] <- macaque_lfc$sig
peak_lfc_matrix[rownames(mouse_lfc),"mouse_lfc"] <- mouse_lfc$avg_log2FC
peak_lfc_matrix[rownames(mouse_lfc),"mouse_sig"] <- mouse_lfc$sig

#try lfc=1, sig=2
peak_lfc_matrix$group <- 'mixed'

temp <- rownames(macaque_lfc)[abs(macaque_lfc$avg_log2FC) > 1 & macaque_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'
temp <- rownames(mouse_lfc)[abs(mouse_lfc$avg_log2FC) > 1 & mouse_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'

peak_lfc_matrix[abs(peak_lfc_matrix$macaque_lfc) <= 1 & 
                  abs(peak_lfc_matrix$mouse_lfc) <= 1 & 
                  countOverlaps(query = Brain_ATAC_peakset$macaque[peak_lfc_matrix$peak],subject = macaque_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$mouse[peak_lfc_matrix$peak],subject = mouse_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$human[peak_lfc_matrix$peak],subject = human_peak) > 0,
                "group"] <- 'conserved'

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peak_lfc_matrix.rds')
saveRDS(object = peak_lfc_matrix,file = char)

#kmeans
mat <- subset_peak_matrix
mat <- NormalizeData(object = mat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
VariableFeatures(mat) <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'differential']
mat <- ScaleData(object = mat,features = VariableFeatures(mat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)

group <- kmeans(mat@assays$RNA@scale.data[VariableFeatures(mat),],centers = 7,nstart = 7)$cluster
col_anno <- HeatmapAnnotation(species = mat$species,which = 'column')

p1 <- Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
              show_column_names = FALSE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(mat@assays$RNA@scale.data)],
                                                      levels = unique(group)),
              cluster_columns = FALSE,
              column_split = factor(mat$species,levels = c('human','macaque','mouse')),
              height = unit(10,'inches'),width = unit(6,'inches'),
              name = 'z-score',top_annotation = col_anno,border = TRUE)

temp <- peak_lfc_matrix[rownames(p1@matrix)[p1@row_order],c("macaque_lfc","mouse_lfc")]
p2 <- Heatmap(matrix = as.matrix(temp),
              show_column_names = TRUE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(temp)],levels = unique(group)),
              cluster_columns = FALSE,
              height = unit(10,'inches'),width = unit(1,'inches'),
              name = 'lfc',border = TRUE)

p1+p2

#modify groups
temp <- peak_lfc_matrix[peak_lfc_matrix$group == 'conserved',]
ii <- rep('conserved',times = nrow(temp))
names(ii) <- temp$peak
group <- append(group,ii)
group <- data.frame(peak_name = names(group),raw_cluster = as.character(group))

group$cluster <- NA
group[group$raw_cluster %in% c('3','4'),"cluster"] <- 'macaque_human_conserved'
group[group$raw_cluster == '6',"cluster"] <- 'macaque_specific'
group[group$raw_cluster == '5',"cluster"] <- 'mouse_macaque_conserved'
group[group$raw_cluster == '2',"cluster"] <- 'mouse_specific'
group[group$raw_cluster == '1',"cluster"] <- 'mouse_human_conserved'
group[group$raw_cluster == '7',"cluster"] <- 'human_specific'
group[group$raw_cluster == 'conserved',"cluster"] <- 'species_conserved'

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peak_kmeans_group.rds')
saveRDS(group,file = char)

#plot
mat <- subset_peak_matrix
mat <- NormalizeData(object = mat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
VariableFeatures(mat) <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'differential']
mat <- ScaleData(object = mat,features = VariableFeatures(mat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)

group <- readRDS(file = char)
ii <- group$peak_name
group <- group$cluster
names(group) <- ii
col_anno <- HeatmapAnnotation(species = mat$species,which = 'column')

p1 <- Heatmap(matrix = as.matrix(mat@assays$RNA@scale.data),
              show_column_names = FALSE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(mat@assays$RNA@scale.data)],
                                                      levels = c('human_specific','macaque_human_conserved',
                                                                 'mouse_human_conserved','macaque_specific',
                                                                 'mouse_macaque_conserved','mouse_specific',
                                                                 'species_conserved')),
              cluster_columns = FALSE,
              column_split = factor(mat$species,levels = c('human','macaque','mouse')),
              height = unit(10,'inches'),width = unit(6,'inches'),
              name = 'z-score',top_annotation = col_anno,border = TRUE,
              row_title = 'merged peakset')

temp <- peak_lfc_matrix[rownames(p1@matrix)[p1@row_order],c("macaque_lfc","mouse_lfc")]
p2 <- Heatmap(matrix = as.matrix(temp),
              show_column_names = TRUE,show_row_names = FALSE,
              cluster_rows = FALSE,row_split = factor(group[rownames(temp)],levels = c('human_specific','macaque_human_conserved',
                                                                                       'mouse_human_conserved','macaque_specific',
                                                                                       'mouse_macaque_conserved','mouse_specific',
                                                                                       'species_conserved')),
              cluster_columns = FALSE,
              height = unit(10,'inches'),width = unit(1,'inches'),
              name = 'lfc',border = TRUE,row_title = 'merged peakset')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_differential_peakset_z_score_heatmap.pdf')
pdf(file = char,width = 10,height = 12)
p1+p2
dev.off()

#boxplot show lfc
peak_lfc_matrix <- peak_lfc_matrix[peak_lfc_matrix$group == 'differential',]
peak_lfc_matrix$cluster <- group[rownames(peak_lfc_matrix)]
temp <- data.frame(lfc = peak_lfc_matrix$macaque_lfc,species = 'macaque',cluster = peak_lfc_matrix$cluster)
temp <- rbind(temp,data.frame(lfc = peak_lfc_matrix$mouse_lfc,species = 'mouse',cluster = peak_lfc_matrix$cluster))

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_differential_peakset_lfc_boxplot.pdf')
pdf(file = char,width = 7,height = 5)
ggplot(data = temp,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  ylim(c(-5,5)) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  facet_wrap(~cluster,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1))
dev.off()

# enrichheatmap
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

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
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

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean,row_order = p1@row_order)

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

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE,raster_resize_mat = mean)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.004,0.008),colors = c('#FFC30F','#C70039','#581845'))
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peakset_CPM_center_besed_heatmap.pdf')
pdf(file = char,width = 4,height = 9)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                         'mouse_human_conserved','macaque_specific',
                                                                                         'mouse_macaque_conserved','mouse_specific',
                                                                                         'species_conserved')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE)
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
p_macaque_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_specific')
p_mouse_specific <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_specific')
p_macaque_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'macaque_human_conserved')
p_mouse_human_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_human_conserved')
p_mouse_macaque_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'mouse_macaque_conserved')
p_species_conserved <- my_insertion_plot(group = group,p1 = p1,p2 = p2,p3 = p3,group.by = 'species_conserved')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 18,height = 7)
p_human_specific+p_macaque_specific+p_mouse_specific+p_macaque_human_conserved+
  p_mouse_human_conserved+p_mouse_macaque_conserved+p_species_conserved + plot_layout(ncol = 4)
dev.off()

#peak contribution
#human
human_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$human[peak_list],subject = human_peak)
human_contributed_peak <- human_contributed_peak[human_contributed_peak > 0]
human_contributed_peak <- names(human_contributed_peak)

#macaque
macaque_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$macaque[peak_list],subject = macaque_peak)
macaque_contributed_peak <- macaque_contributed_peak[macaque_contributed_peak > 0]
macaque_contributed_peak <- names(macaque_contributed_peak)

#mouse
mouse_contributed_peak <- countOverlaps(query = Brain_ATAC_peakset$mouse[peak_list],subject = mouse_peak)
mouse_contributed_peak <- mouse_contributed_peak[mouse_contributed_peak > 0]
mouse_contributed_peak <- names(mouse_contributed_peak)

#ggvenn
temp <- list(human = human_contributed_peak,macaque = macaque_contributed_peak,mouse = mouse_contributed_peak)
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/species_contribution_to_final_',char,'_peakset.pdf')
pdf(file = char,width = 6.5,height = 6.5)
ggvenn(data = temp,c('human','macaque','mouse'))
dev.off()

#sankey plot
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_peak_lfc_matrix.rds')
peak_lfc_matrix <- readRDS(file = char)
peak_lfc_matrix$kmeans_cluster <- 'dropped'
peak_lfc_matrix[names(group),"kmeans_cluster"] <- as.character(group)

temp <- base::lapply(X = peak_lfc_matrix$peak,FUN = function(x){
  contribution_score <- 1*sum(x %in% mouse_contributed_peak) + 10*sum(x %in% macaque_contributed_peak) + 100*sum(x %in% human_contributed_peak)
  if(contribution_score == 1){
    return('mouse_specific')
  }else if(contribution_score == 10){
    return('macaque_specific')
  }else if(contribution_score == 100){
    return('human_specific')
  }else if(contribution_score == 11){
    return('mouse_macaque_conserved')
  }else if(contribution_score == 101){
    return('mouse_human_conserved')
  }else if(contribution_score == 110){
    return('macaque_human_conserved')
  }else if(contribution_score == 111){
    return('species_conserved')
  }else{
    return('error!')
  }
})
temp <- unlist(temp)
peak_lfc_matrix$peak_cluster <- temp

#scibet
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220529/',char,'_kmeans_cluster_and_peak_only_cluster_confusion_heatmap.pdf')
pdf(file = char,width = 6,height = 6)
scibet::Confusion_heatmap(ori = peak_lfc_matrix$peak_cluster,prd = peak_lfc_matrix$kmeans_cluster) + 
  theme(aspect.ratio = 1,
        axis.title = element_text(size = 16,face = 'bold')) + 
  ylab('kmeans cluster') + xlab('peak only cluster')
dev.off()