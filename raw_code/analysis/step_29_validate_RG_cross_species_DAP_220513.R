#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: validate RG cross species DAP                                   ##
## Data: 2022.05.13                                                                ##
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
library(EnrichedHeatmap)
library(circlize)
library(scales)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
Brain_ATAC_peakset <- readRDS(file = './res/step_27_fig_220510/Brain_ATAC_peakset.rds')
Brain_peak_matrix <- readRDS(file = './res/step_27_fig_220510/Brain_peak_matrix_Seurat.rds')

macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')

# data distribution of RG_peak_matrix ----------------------------------
RG_peak_matrix <- Brain_peak_matrix[,Brain_peak_matrix$cell_type == 'RG']

#get peak_list
human_RG_peak <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/RG-reproduciblePeaks.gr.rds')
temp <- countOverlaps(query = Brain_ATAC_peakset$human,subject = human_RG_peak)
temp <- temp[temp > 0]
peak_list <- names(temp)

macaque_RG_peak <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/RG-reproduciblePeaks.gr.rds')
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
temp <- countOverlaps(query = Brain_ATAC_peakset$macaque,subject = macaque_RG_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

mouse_RG_peak <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/RG-reproduciblePeaks.gr.rds')
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
temp <- countOverlaps(query = Brain_ATAC_peakset$mouse,subject = mouse_RG_peak)
temp <- temp[temp > 0]
temp <- names(temp)
peak_list <- append(peak_list,temp)

peak_list <- unique(peak_list)
table(peak_list %in% rownames(RG_peak_matrix))

# macaque lfc -------------------------------------------------------------
macaque_lfc <- FindMarkers(object = RG_peak_matrix,ident.1 = 'macaque',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
macaque_lfc$sig <- -log10(macaque_lfc$p_val_adj)
macaque_lfc[macaque_lfc$sig > 300,"sig"] <- 300

# mouse lfc ---------------------------------------------------------------
mouse_lfc <- FindMarkers(object = RG_peak_matrix,ident.1 = 'mouse',ident.2 = 'human',group.by = 'species',assay = 'RNA',slot = 'counts',features = peak_list,test.use = 'DESeq2',verbose = TRUE,only.pos = FALSE)
mouse_lfc$sig <- -log10(mouse_lfc$p_val_adj)
mouse_lfc[mouse_lfc$sig > 300,"sig"] <- 300

# peak lfc matrix ---------------------------------------------------------
peak_lfc_matrix <- data.frame(peak = peak_list,macaque_lfc = 0,macaque_sig = 0,mouse_lfc = 0,mouse_sig = 0)
rownames(peak_lfc_matrix) <- peak_lfc_matrix$peak
peak_lfc_matrix[rownames(macaque_lfc),"macaque_lfc"] <- macaque_lfc$avg_log2FC
peak_lfc_matrix[rownames(macaque_lfc),"macaque_sig"] <- macaque_lfc$sig
peak_lfc_matrix[rownames(mouse_lfc),"mouse_lfc"] <- mouse_lfc$avg_log2FC
peak_lfc_matrix[rownames(mouse_lfc),"mouse_sig"] <- mouse_lfc$sig

#try lfc=1/0.7, sig=2
peak_lfc_matrix$group <- 'mixed'

temp <- rownames(macaque_lfc)[abs(macaque_lfc$avg_log2FC) > 1 & macaque_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'
temp <- rownames(mouse_lfc)[abs(mouse_lfc$avg_log2FC) > 1 & mouse_lfc$sig > 2]
peak_lfc_matrix[temp,"group"] <- 'differential'

peak_lfc_matrix[abs(peak_lfc_matrix$macaque_lfc) < 0.7 & 
                  abs(peak_lfc_matrix$mouse_lfc) < 0.7 & 
                  countOverlaps(query = Brain_ATAC_peakset$macaque[peak_lfc_matrix$peak],subject = macaque_RG_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$mouse[peak_lfc_matrix$peak],subject = mouse_RG_peak) > 0 & 
                  countOverlaps(query = Brain_ATAC_peakset$human[peak_lfc_matrix$peak],subject = human_RG_peak) > 0,
                "group"] <- 'conserved'

#kmeans
mat <- RG_peak_matrix
mat <- NormalizeData(object = mat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000,verbose = TRUE)
VariableFeatures(mat) <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'differential']
mat <- ScaleData(object = mat,features = VariableFeatures(mat),assay = 'RNA',vars.to.regress = NULL,do.scale = TRUE,do.center = TRUE,verbose = TRUE)


# group <- kmeans(mat@assays$RNA@scale.data[VariableFeatures(mat),],centers = 6)$cluster
group <- readRDS(file = './res/step_29_fig_220513/RG_peak_kmeans_group.rds')
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
              height = unit(8,'inches'),width = unit(6,'inches'),
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
              height = unit(8,'inches'),width = unit(1,'inches'),
              name = 'lfc',border = TRUE,row_title = 'merged peakset')

pdf(file = './res/step_29_fig_220513/RG_differential_peakset_z_score_heatmap.pdf',width = 10,height = 10)
p1+p2
dev.off()

#boxplot show lfc
peak_lfc_matrix <- peak_lfc_matrix[peak_lfc_matrix$group == 'differential',]
peak_lfc_matrix$cluster <- group[rownames(peak_lfc_matrix)]
temp <- data.frame(lfc = peak_lfc_matrix$macaque_lfc,species = 'macaque',cluster = peak_lfc_matrix$cluster)
temp <- rbind(temp,data.frame(lfc = peak_lfc_matrix$mouse_lfc,species = 'mouse',cluster = peak_lfc_matrix$cluster))

pdf(file = './res/step_29_fig_220513/RG_differential_peakset_lfc_boxplot.pdf',width = 7,height = 5)
ggplot(data = temp,aes(x=species,y=lfc,fill=species)) + 
  geom_boxplot(outlier.alpha = 0,width = 0.5) + 
  geom_hline(yintercept = 0,linetype = 'dashed',color = 'red') + 
  facet_wrap(~cluster,ncol = 3) + 
  theme_classic() + 
  scale_fill_manual(values = c('macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme(aspect.ratio = 1.2,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1))
dev.off()

#cluster 1 mouse_macaque_conserved
#cluster 2 mouse_human_conserved
#cluster 3 macaque_human_conserved
#cluster 4 macaque_specific
#cluster 5 mouse_specific
#cluster 6 human_specific

# ii <- rownames(peak_lfc_matrix)[peak_lfc_matrix$group == 'conserved']
# temp <- rep(x = 'conserved',times = length(ii))
# names(temp) <- ii
# group <- append(group,temp)
# 
# group <- data.frame(peak_name = names(group),raw_cluster = as.character(group))
# rownames(group) <- group$peak_name
# group$cluster <- NA
# group[group$raw_cluster == '1',"cluster"] <- 'mouse_macaque_conserved'
# group[group$raw_cluster == '2',"cluster"] <- 'mouse_human_conserved'
# group[group$raw_cluster == '3',"cluster"] <- 'macaque_human_conserved'
# group[group$raw_cluster == '4',"cluster"] <- 'macaque_specific'
# group[group$raw_cluster == '5',"cluster"] <- 'mouse_specific'
# group[group$raw_cluster == '6',"cluster"] <- 'human_specific'
# group[group$raw_cluster == 'conserved',"cluster"] <- 'species_conserved'
# saveRDS(object = group,file = './res/step_29_fig_220513/RG_peak_kmeans_group.rds')

# enriched heatmap --------------------------------------------------------

#human
target_site <- Brain_ATAC_peakset$human[names(group)]
signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/human_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE)

#macaque
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
target_site <- Brain_ATAC_peakset$macaque[names(group)]
signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/macaque_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE)

p1+p2

#mouse
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
target_site <- Brain_ATAC_peakset$mouse[names(group)]
signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/mouse_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.003,0.006),colors = c('#FFC30F','#C70039','#581845'))
pdf(file = './res/step_29_fig_220513/RG_peakset_CPM_heatmap.pdf',width = 4,height = 8)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                         'mouse_human_conserved','macaque_specific',
                                                                                         'mouse_macaque_conserved','mouse_specific',
                                                                                         'species_conserved')),
                use_raster = TRUE,col = col_fun,name = 'insertion') + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE)
dev.off()


# center based enrich heatmap ---------------------------------------------

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

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/human_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE)

#macaque
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

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/macaque_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE)

p1+p2

#mouse
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

signal_coverage <- rtracklayer::import.bw(con = './res/step_28_fig_220511/mouse_cell_type_coverage/RG-TileSize-25-normMethod-nFrags-ArchR.bw')
mat <- normalizeToMatrix(signal = signal_coverage,target = target_site,extend = 5000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p3 <- EnrichedHeatmap(mat = mat,row_split = factor(group[rownames(mat)],levels = c('human_specific','macaque_human_conserved',
                                                                                   'mouse_human_conserved','macaque_specific',
                                                                                   'mouse_macaque_conserved','mouse_specific',
                                                                                   'species_conserved')),
                      use_raster = TRUE)

p1+p2+p3

#plot
col_fun <- colorRamp2(breaks = c(0,0.003,0.006),colors = c('#FFC30F','#C70039','#581845'))
pdf(file = './res/step_29_fig_220513/RG_peakset_CPM_center_besed_heatmap.pdf',width = 4,height = 8)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(group[rownames(p1@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                         'mouse_human_conserved','macaque_specific',
                                                                                         'mouse_macaque_conserved','mouse_specific',
                                                                                         'species_conserved')),
                use_raster = TRUE,col = col_fun,name = 'insertion',pos_line = FALSE) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(group[rownames(p2@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE) + 
  EnrichedHeatmap(mat = p3@matrix,row_split = factor(group[rownames(p3@matrix)],levels = c('human_specific','macaque_human_conserved',
                                                                                           'mouse_human_conserved','macaque_specific',
                                                                                           'mouse_macaque_conserved','mouse_specific',
                                                                                           'species_conserved')),
                  use_raster = TRUE,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE)
dev.off()

# enrichment distribution -------------------------------------------------
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

pdf(file = './res/step_29_fig_220513/RG_peakset_mean_insertion_dot_line_plot.pdf',width = 18,height = 7)
p_human_specific+p_macaque_specific+p_mouse_specific+p_macaque_human_conserved+
  p_mouse_human_conserved+p_mouse_macaque_conserved+p_species_conserved + plot_layout(ncol = 4)
dev.off()

# putative downstream genes -----------------------------------------------
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_27_fig_220510/Brain_ATAC_peakset.rds')
names(Brain_ATAC_peakset$macaque) <- Brain_ATAC_peakset$macaque$name
names(Brain_ATAC_peakset$mouse) <- Brain_ATAC_peakset$mouse$name
Brain_peak_matrix <- readRDS(file = './res/step_27_fig_220510/Brain_peak_matrix_Seurat.rds')

group <- readRDS(file = './res/step_29_fig_220513/RG_peak_kmeans_group.rds')
ii <- group$peak_name
group <- group$cluster
names(group) <- ii

human_p2g <- readRDS(file = './res/step_28_fig_220511/human_p2g.rds')
macaque_p2g <- readRDS(file = './res/step_28_fig_220511/macaque_p2g.rds')
mouse_p2g <- readRDS(file = './res/step_28_fig_220511/mouse_p2g.rds')

macaque_integration_Seurat <- readRDS(file = './processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/greenleaf_human_cortex_RNA_seurat_reanno_210728.rds')
mouse_RNA_Seurat <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scRNA_seq/processed_data/mouse_RNA_Seurat_pcaproject_label_transfer_Seuart_220321.rds')
#create function
count_target <- function(peak_list,p2g){
  #select peak anchor
  p2g@metadata$peakSet$idx <- 1:length(p2g@metadata$peakSet)
  p2g@metadata$peakSet$name <- paste(p2g@metadata$peakSet@seqnames,as.character(p2g@metadata$peakSet@ranges),sep = '-')
  temp <- p2g@metadata$peakSet[p2g@metadata$peakSet$name %in% peak_list]
  #count gene
  gene_list <- p2g[p2g$idxATAC %in% temp$idx,]
  gene_list <- p2g@metadata$geneSet$name[gene_list$idxRNA]
  #create data.frame
  temp <- data.frame(table(gene_list))
  colnames(temp) <- c('gene','counts')
  temp$counts <- temp$counts/length(peak_list)
  rownames(temp) <- temp$gene
  temp$p <- NA
  for (i in rownames(temp)) {
    temp[i,'p'] <- dhyper(x = sum(gene_list == i),
                          m = sum(p2g@metadata$geneSet$name[p2g$idxRNA] == i),
                          n = sum(p2g@metadata$geneSet$name[p2g$idxRNA] != i),
                          k = length(gene_list))
  }
  temp$p <- -log10(temp$p)
  #return
  return(temp)
}

## human_specific_target ---------------------------------------------------
human_specific_target <- count_target(peak_list = names(group)[group == 'human_specific'],p2g = human_p2g)

#GO analysis
gene_list <- rownames(Greenleaf_RNA_Seurat@assays$converted@counts)
gene_list <- c(gene_list %in% human_specific_target$gene)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(Greenleaf_RNA_Seurat@assays$converted@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human_specific",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = allRes$Term)

pdf(file = './res/step_29_fig_220513/RG_human_specific_peak_target_gene_BP.pdf',width = 6,height = 4)
ggplot(allRes[c(1,2,3,5,7,14,15,16,24,26),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#87D2DB',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.8,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('human_specific BP')
dev.off()

#dotplot
human_specific_target$group <- 'human_specific'
for (i in 1:10) {
  temp <- sample(x = paste(human_p2g@metadata$peakSet@seqnames,as.character(human_p2g@metadata$peakSet@ranges),sep = '-'),size = sum(group == 'human_specific'))
  temp <- count_target(peak_list = temp,p2g = human_p2g)
  temp$group <- 'background'
  human_specific_target <- rbind(human_specific_target,temp)
}

human_specific_target$label <- human_specific_target$gene
human_specific_target[human_specific_target$group == 'background',"label"] <- NA
human_specific_target[human_specific_target$p < 2,"label"] <- NA

pdf(file = './res/step_29_fig_220513/RG_human_specific_peak_target_gene_dotplot.pdf',width = 10,height = 10)
ggplot(data = human_specific_target,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','human_specific' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')
dev.off()


## macaque_specific_target --------------------------------------------------------
macaque_specific_target <- names(group)[group == 'macaque_specific']
macaque_specific_target <- Brain_ATAC_peakset$macaque[macaque_specific_target]
macaque_specific_target <- paste(macaque_specific_target@seqnames,as.character(macaque_specific_target@ranges),sep = '-')
macaque_specific_target <- count_target(peak_list = macaque_specific_target,p2g = macaque_p2g)

#GO analysis
gene_list <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- c(gene_list %in% macaque_specific_target$gene)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "macaque_specific",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = unique(allRes$Term))

pdf(file = './res/step_29_fig_220513/RG_macaque_specific_peak_target_gene_BP.pdf',width = 6,height = 4)
ggplot(allRes[c(1,2,3,4,6,8,9,12,13,15),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#87D2DB',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.8,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('macaque_specific BP')
dev.off()

#dotplot
macaque_specific_target$group <- 'macaque_specific'
for (i in 1:10) {
  temp <- sample(x = paste(macaque_p2g@metadata$peakSet@seqnames,as.character(macaque_p2g@metadata$peakSet@ranges),sep = '-'),size = sum(group == 'macaque_specific'))
  temp <- count_target(peak_list = temp,p2g = macaque_p2g)
  temp$group <- 'background'
  macaque_specific_target <- rbind(macaque_specific_target,temp)
}

macaque_specific_target$label <- macaque_specific_target$gene
macaque_specific_target[macaque_specific_target$group == 'background',"label"] <- NA
macaque_specific_target[macaque_specific_target$p < 2,"label"] <- NA

pdf(file = './res/step_29_fig_220513/RG_macaque_specific_peak_target_gene_dotplot.pdf',width = 10,height = 10)
ggplot(data = macaque_specific_target,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','macaque_specific' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')
dev.off()

## macaque_human_conserved_target ------------------------------------------

#in human
macaque_human_conserved_target_in_human <- count_target(peak_list = names(group)[group == 'macaque_human_conserved'],p2g = human_p2g)

#GO analysis
gene_list <- rownames(Greenleaf_RNA_Seurat@assays$converted@counts)
gene_list <- c(gene_list %in% macaque_human_conserved_target_in_human$gene)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(Greenleaf_RNA_Seurat@assays$converted@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "macaque_human_conserved_in_human",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = unique(allRes$Term))

pdf(file = './res/step_29_fig_220513/RG_macaque_human_conserved_in_human_peak_target_gene_BP.pdf',width = 6,height = 4)
ggplot(allRes[c(1,3,4,5,6,7,8,9,10,12),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#87D2DB',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.8,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('macaque_human_conserved_in_human BP')
dev.off()

#dotplot
macaque_human_conserved_target_in_human$group <- 'macaque_human_conserved'
for (i in 1:10) {
  temp <- sample(x = paste(human_p2g@metadata$peakSet@seqnames,as.character(human_p2g@metadata$peakSet@ranges),sep = '-'),size = sum(group == 'macaque_human_conserved'))
  temp <- count_target(peak_list = temp,p2g = human_p2g)
  temp$group <- 'background'
  macaque_human_conserved_target_in_human <- rbind(macaque_human_conserved_target_in_human,temp)
}

macaque_human_conserved_target_in_human$label <- macaque_human_conserved_target_in_human$gene
macaque_human_conserved_target_in_human[macaque_human_conserved_target_in_human$group == 'background',"label"] <- NA
macaque_human_conserved_target_in_human[macaque_human_conserved_target_in_human$p < 2,"label"] <- NA

pdf(file = './res/step_29_fig_220513/RG_macaque_human_conserved_in_human_peak_target_gene_dotplot.pdf',width = 10,height = 10)
ggplot(data = macaque_human_conserved_target_in_human,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','macaque_human_conserved' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')
dev.off()

#in macaque
macaque_human_conserved_target_in_macaque <- names(group)[group == 'macaque_human_conserved']
macaque_human_conserved_target_in_macaque <- Brain_ATAC_peakset$macaque[macaque_human_conserved_target_in_macaque]
macaque_human_conserved_target_in_macaque <- paste(macaque_human_conserved_target_in_macaque@seqnames,as.character(macaque_human_conserved_target_in_macaque@ranges),sep = '-')
macaque_human_conserved_target_in_macaque <- count_target(peak_list = macaque_human_conserved_target_in_macaque,p2g = macaque_p2g)

#GO analysis
gene_list <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- c(gene_list %in% macaque_specific_target$gene)
gene_list[gene_list == TRUE] <- 1
gene_list[gene_list == FALSE] <- 0
names(gene_list) <- rownames(macaque_integration_Seurat@assays$RNA@counts)
gene_list <- factor(gene_list,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "macaque_human_conserved_in_macaque",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes$Term <- factor(allRes$Term,levels = unique(allRes$Term))

pdf(file = './res/step_29_fig_220513/RG_macaque_human_conserved_in_macaque_peak_target_gene_BP.pdf',width = 6,height = 4)
ggplot(allRes[c(1,2,3,6,8,9,11,13,14,16),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = '#87D2DB',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.8,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',size = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('macaque_human_conserved_in_macaque BP')
dev.off()

#dotplot
macaque_human_conserved_target_in_macaque$group <- 'macaque_human_conserved'
for (i in 1:10) {
  temp <- sample(x = paste(macaque_p2g@metadata$peakSet@seqnames,as.character(macaque_p2g@metadata$peakSet@ranges),sep = '-'),size = sum(group == 'macaque_human_conserved'))
  temp <- count_target(peak_list = temp,p2g = macaque_p2g)
  temp$group <- 'background'
  macaque_human_conserved_target_in_macaque <- rbind(macaque_human_conserved_target_in_macaque,temp)
}

macaque_human_conserved_target_in_macaque$label <- macaque_human_conserved_target_in_macaque$gene
macaque_human_conserved_target_in_macaque[macaque_human_conserved_target_in_macaque$group == 'background',"label"] <- NA
macaque_human_conserved_target_in_macaque[macaque_human_conserved_target_in_macaque$p < 2,"label"] <- NA

pdf(file = './res/step_29_fig_220513/RG_macaque_human_conserved_in_macaque_peak_target_gene_dotplot.pdf',width = 10,height = 10)
ggplot(data = macaque_human_conserved_target_in_macaque,aes(x=counts,y=p,color=group)) + 
  theme_cowplot() + 
  geom_point(size = 2) + 
  scale_color_manual(values = c('background' = 'blue','macaque_human_conserved' = 'red')) + 
  geom_text_repel(aes(label=label),color='black') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 14,face = 'bold'),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,colour = 'black',size = 1)) + 
  labs(title = 'target gene counts') + 
  xlab('count times') + ylab('-log10(p-value)')
dev.off()