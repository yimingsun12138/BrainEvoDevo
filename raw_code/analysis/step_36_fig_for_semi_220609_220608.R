#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: fig for semi 220609                                             ##
## Data: 2022.06.08                                                                ##
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

# UMAP for scATAC data ----------------------------------------------------
#load color
color_param <- read.csv(file = './data/parameter/scATACseq_color_param.csv')
col_value <- color_param$color
names(col_value) <- color_param$term

color_param <- read.csv(file = './data/parameter/scRNAseq_color_param.csv')
col_value <- append(col_value,c('End' = '#e0598b','Per' = '#342739','Cycling' = '#A9A911FF'))

#human
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')

pdf(file = './res/step_36_fig_220608/human_scATAC_seq_umap_dimplot.pdf',width = 5,height = 5)
my_dimplot(embedding = Greenleaf_ATAC_ArchR@embeddings$UMAP$df,
           meta_data = as.data.frame(Greenleaf_ATAC_ArchR@cellColData),
           group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_cowplot() + 
  scale_color_manual(values = col_value[unique(Greenleaf_ATAC_ArchR$cell_type)]) + 
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'human scATAC-seq')
dev.off()

#macaque
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')

pdf(file = './res/step_36_fig_220608/macaque_scATAC_seq_umap_dimplot.pdf',width = 5,height = 5)
my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,
           meta_data = as.data.frame(macaque_multiome_ArchR@cellColData),
           group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_cowplot() + 
  scale_color_manual(values = col_value[unique(macaque_multiome_ArchR$cell_type)]) + 
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'macaque scATAC-seq')
dev.off()

#mouse
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')

pdf(file = './res/step_36_fig_220608/mouse_scATAC_seq_umap_dimplot.pdf',width = 5,height = 5)
my_dimplot(embedding = mouse_ATAC_ArchR@embeddings$UMAP$df,
           meta_data = as.data.frame(mouse_ATAC_ArchR@cellColData),
           group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme_cowplot() + 
  scale_color_manual(values = col_value[unique(mouse_ATAC_ArchR$cell_type)]) + 
  theme(aspect.ratio = 1,
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1),
        axis.text = element_blank(),
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'mouse scATAC-seq')
dev.off()


# cell type peak number ---------------------------------------------------

#human
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
peak_path <- './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/PeakCalls/'
peak_file <- list.files(peak_path)
peak_file <- peak_file[grep(pattern = '.gr.rds$',x = peak_file,fixed = FALSE)]

temp <- unlist(base::lapply(X = peak_file,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][1])
}))
temp <- sub(pattern = '.',replacement = '-',x = temp,fixed = TRUE)

peak_file <- data.frame(cell_type = c(temp,'merged'),path = c(paste0(peak_path,peak_file),NA),num = 0)


for (i in 1:(nrow(peak_file)-1)) {
  peak_path <- peak_file[i,"path"]
  temp_peakset <- readRDS(peak_path)
  peak_file[i,"num"] <- length(temp_peakset)
}

peak_file[peak_file$cell_type == 'merged',"num"] <- length(getPeakSet(Greenleaf_ATAC_ArchR))
peak_file$num <- peak_file$num/1000
peak_file$cell_type <- factor(peak_file$cell_type,levels = c('End','Per','Mic','OPC','RG','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE','merged'))

pdf(file = './res/step_36_fig_220608/human_scATAC_peaks_barplot.pdf',width = 5,height = 4)
ggplot(data = peak_file,aes(x = cell_type,y = num)) + 
  geom_bar(stat="identity",width = 0.7,fill = 'steelblue') + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_text(size = 14,face = 'bold'),
        panel.grid.minor = element_line(size = 0.5,color = 'grey',linetype = 'dashed'),
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Cell Type') + ylab('peak num (k)') + labs(title = 'human scATAC-seq peaks')
dev.off()

human_peak_file <- peak_file
human_peak_file$species <- 'human'

#macaque
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
peak_path <- './ArchR/processed_data/macaque_multiome_ArchR_220411/PeakCalls/'
peak_file <- list.files(peak_path)
peak_file <- peak_file[grep(pattern = '.gr.rds$',x = peak_file,fixed = FALSE)]

temp <- unlist(base::lapply(X = peak_file,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][1])
}))
temp <- sub(pattern = '.',replacement = '-',x = temp,fixed = TRUE)
temp[which(temp == 'End-Per')] <- 'End/Per'

peak_file <- data.frame(cell_type = c(temp,'merged'),path = c(paste0(peak_path,peak_file),NA),num = 0)


for (i in 1:(nrow(peak_file)-1)) {
  peak_path <- peak_file[i,"path"]
  temp_peakset <- readRDS(peak_path)
  peak_file[i,"num"] <- length(temp_peakset)
}

peak_file[peak_file$cell_type == 'merged',"num"] <- length(getPeakSet(macaque_multiome_ArchR))
peak_file$num <- peak_file$num/1000
peak_file$cell_type <- factor(peak_file$cell_type,levels = c('End/Per','Mic','OPC','Ependymal','RG','Cycling','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE','merged'))

pdf(file = './res/step_36_fig_220608/macaque_scATAC_peaks_barplot.pdf',width = 5,height = 4)
ggplot(data = peak_file,aes(x = cell_type,y = num)) + 
  geom_bar(stat="identity",width = 0.7,fill = 'steelblue') + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_text(size = 14,face = 'bold'),
        panel.grid.minor = element_line(size = 0.5,color = 'grey',linetype = 'dashed'),
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Cell Type') + ylab('peak num (k)') + labs(title = 'macaque scATAC-seq peaks')
dev.off()

macaque_peak_file <- peak_file
macaque_peak_file$species <- 'macaque'

#mouse
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')
peak_path <- './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/PeakCalls/'
peak_file <- list.files(peak_path)
peak_file <- peak_file[grep(pattern = '.gr.rds$',x = peak_file,fixed = FALSE)]

temp <- unlist(base::lapply(X = peak_file,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][1])
}))
temp <- sub(pattern = '.',replacement = '-',x = temp,fixed = TRUE)

peak_file <- data.frame(cell_type = c(temp,'merged'),path = c(paste0(peak_path,peak_file),NA),num = 0)


for (i in 1:(nrow(peak_file)-1)) {
  peak_path <- peak_file[i,"path"]
  temp_peakset <- readRDS(peak_path)
  peak_file[i,"num"] <- length(temp_peakset)
}

peak_file[peak_file$cell_type == 'merged',"num"] <- length(getPeakSet(mouse_ATAC_ArchR))
peak_file$num <- peak_file$num/1000
peak_file$cell_type <- factor(peak_file$cell_type,levels = c('RG','IP','Ex-1','Ex-3','merged'))

pdf(file = './res/step_36_fig_220608/mouse_scATAC_peaks_barplot.pdf',width = 5,height = 4)
ggplot(data = peak_file,aes(x = cell_type,y = num)) + 
  geom_bar(stat="identity",width = 0.7,fill = 'steelblue') + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_text(size = 14,face = 'bold'),
        panel.grid.minor = element_line(size = 0.5,color = 'grey',linetype = 'dashed'),
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Cell Type') + ylab('peak num (k)') + labs(title = 'mouse scATAC-seq peaks')
dev.off()

mouse_peak_file <- peak_file
mouse_peak_file$species <- 'mouse'

#merged plot
human_peak_file <- human_peak_file[!(human_peak_file$cell_type %in% c('End','Per')),]
macaque_peak_file <- macaque_peak_file[macaque_peak_file$cell_type %in% human_peak_file$cell_type,]
mouse_peak_file <- mouse_peak_file[mouse_peak_file$cell_type %in% human_peak_file$cell_type,]
mouse_peak_file <- rbind(mouse_peak_file,data.frame(cell_type = c('Ex-2','InCGE','InMGE','Mic','OPC'),path = NA,num = 0,species = 'mouse'))

peak_file <- rbind(human_peak_file,macaque_peak_file,mouse_peak_file)
peak_file$cell_type <- factor(peak_file$cell_type,levels = c('Mic','OPC','RG','IP','Ex-1','Ex-2','Ex-3','InMGE','InCGE','merged'))

pdf(file = './res/step_36_fig_220608/species_scATAC_peak_num_barplot.pdf',width = 7,height = 4)
ggplot(data = peak_file,aes(x = cell_type,y = num,fill = species)) + 
  geom_bar(stat="identity",position=position_dodge(),width = 0.8) + 
  scale_fill_manual(values = c('human' = '#A31D1D','macaque' = '#3361A5','mouse' = '#FDB31A')) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.title = element_text(size = 14,face = 'bold'),
        panel.grid.minor = element_line(size = 0.5,color = 'grey',linetype = 'dashed'),
        plot.title = element_text(hjust = 0.5)) + 
  xlab('Cell Type') + ylab('peak num (k)') + labs(title = 'scATAC-seq peaks')
dev.off()

# venn plot of species peaks per cell -------------------------------------

## IP DAP ------------------------------------------------------------------
cell_type <- 'IP'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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
table(peak_list %in% names(Brain_ATAC_peakset$human))

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
char <- paste0('./res/step_36_fig_220608/species_contribution_to_final_',char,'_peakset.pdf')
pdf(file = char,width = 6.5,height = 6.5)
ggvenn(data = temp,c('human','macaque','mouse'))
dev.off()

## Ex-1 DAP ------------------------------------------------------------------
cell_type <- 'Ex-1'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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
table(peak_list %in% names(Brain_ATAC_peakset$human))

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
char <- paste0('./res/step_36_fig_220608/species_contribution_to_final_',char,'_peakset.pdf')
pdf(file = char,width = 6.5,height = 6.5)
ggvenn(data = temp,c('human','macaque','mouse'))
dev.off()

## Ex-3 DAP ------------------------------------------------------------------
cell_type <- 'Ex-3'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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
table(peak_list %in% names(Brain_ATAC_peakset$human))

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
char <- paste0('./res/step_36_fig_220608/species_contribution_to_final_',char,'_peakset.pdf')
pdf(file = char,width = 6.5,height = 6.5)
ggvenn(data = temp,c('human','macaque','mouse'))
dev.off()

## RG DAP ------------------------------------------------------------------
cell_type <- 'RG'
#load data
Brain_ATAC_peakset <- readRDS(file = './res/step_32_fig_220603/Brain_ATAC_peakset.rds')

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
table(peak_list %in% names(Brain_ATAC_peakset$human))

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
char <- paste0('./res/step_36_fig_220608/species_contribution_to_final_',char,'_peakset.pdf')
pdf(file = char,width = 6.5,height = 6.5)
ggvenn(data = temp,c('human','macaque','mouse'))
dev.off()

# validate DAP ------------------------------------------------------------

## Ex-3 DAP ------------------------------------------------------------------
cell_type <- 'Ex-3'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
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

#Reads in cell type peaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"Reads_In_Cell_Type_Peaks"])*median(mat$Reads_In_Cell_Type_Peaks)
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
char <- paste0('./res/step_36_fig_220608/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_CPM_center_besed_heatmap.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

## Ex-1 DAP ------------------------------------------------------------------
cell_type <- 'Ex-1'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
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

#Reads in cell type peaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"Reads_In_Cell_Type_Peaks"])*median(mat$Reads_In_Cell_Type_Peaks)
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
char <- paste0('./res/step_36_fig_220608/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_CPM_center_besed_heatmap.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

## IP DAP ------------------------------------------------------------------
cell_type <- 'IP'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
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

#Reads in cell type peaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"Reads_In_Cell_Type_Peaks"])*median(mat$Reads_In_Cell_Type_Peaks)
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
char <- paste0('./res/step_36_fig_220608/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_CPM_center_besed_heatmap.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()

## RG DAP ------------------------------------------------------------------
cell_type <- 'RG'
#load data
Brain_peak_matrix <- readRDS(file = './res/step_35_fig_220606/Brain_peak_matrix_Seurat.rds')
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

#Reads in cell type peaks normalized lfc
char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_33_fig_220605/',char,'_peak_group.rds')
group <- readRDS(file = char)

mat <- subset_peak_matrix
temp <- mat@assays$RNA@counts[names(group),]
for (i in colnames(temp)) {
  temp[,i] <- (temp[,i]/mat@meta.data[i,"Reads_In_Cell_Type_Peaks"])*median(mat$Reads_In_Cell_Type_Peaks)
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
char <- paste0('./res/step_36_fig_220608/',char,'_peak_matrix_lfc_normalized_by_Reads_In_Cell_Type_Peaks_boxplot.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_CPM_center_besed_heatmap.pdf')
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
char <- paste0('./res/step_36_fig_220608/',char,'_peakset_mean_insertion_dot_line_plot.pdf')
pdf(file = char,width = 12,height = 3)
p_human_specific+p_human_macaque_conserved+p_species_conserved+plot_layout(ncol = 3)
dev.off()