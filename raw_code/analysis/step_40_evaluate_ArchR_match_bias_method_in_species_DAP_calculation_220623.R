#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: evaluate ArchR match bias method in species DAP calculation     ##
## Data: 2022.06.23                                                                ##
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

# load data ---------------------------------------------------------------
macaque_multiome_ArchR <- loadArchRProject(path = '/data/User/sunym/project/Brain/ArchR/processed_data/macaque_multiome_ArchR_220411')

macaque_peak_matrix <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
human_peak_matrix <- readRDS(file = './res/step_35_fig_220606/human_peak_matrix_Seurat.rds')
macaque_peak_matrix$species <- 'macaque'
human_peak_matrix$species <- 'human'

# intra-species comparison -----------------------------------------------

## RG vs Ex-1 ----------------------------------------------------------
cell_type_1 <- 'RG'
cell_type_2 <- 'Ex-1'

matchobj <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                         groups = macaque_multiome_ArchR$cell_type,
                                         useGroups = cell_type_1,
                                         bgdGroups = cell_type_2,
                                         bias = c("TSSEnrichment", "log10(nFrags)"),
                                         k = 100,n = 1000,seed = 1,bufferRatio = 0.8)

#TSSenrichment distribution
p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c(cell_type_1,cell_type_2),]),
             aes(x = TSSEnrichment,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution')

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),]),
             aes(x = TSSEnrichment,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution')

char <- sub(pattern = '-',replacement = '_',x = cell_type_1,fixed = TRUE)
char <- paste(char,'vs',sub(pattern = '-',replacement = '_',x = cell_type_2,fixed = TRUE),sep = '_')
char <- paste0('./res/step_40_fig_220623/',char,'_intra_species_TSSEnrichment_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

#nFrags distribution
macaque_multiome_ArchR$log_nFrags <- log10(macaque_multiome_ArchR$nFrags)

p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c(cell_type_1,cell_type_2),]),
             aes(x = log_nFrags,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution') + xlab('log10(nFrags)')

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),]),
             aes(x = log_nFrags,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution') + xlab('log10(nFrags)')

char <- sub(pattern = '-',replacement = '_',x = cell_type_1,fixed = TRUE)
char <- paste(char,'vs',sub(pattern = '-',replacement = '_',x = cell_type_2,fixed = TRUE),sep = '_')
char <- paste0('./res/step_40_fig_220623/',char,'_intra_species_nFrags_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

## InCGE vs InMGE ----------------------------------------------------------
cell_type_1 <- 'InCGE'
cell_type_2 <- 'InMGE'

matchobj <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                         groups = macaque_multiome_ArchR$cell_type,
                                         useGroups = cell_type_1,
                                         bgdGroups = cell_type_2,
                                         bias = c("TSSEnrichment", "log10(nFrags)"),
                                         k = 100,n = 1000,seed = 1,bufferRatio = 0.8)

#TSSenrichment distribution
p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c(cell_type_1,cell_type_2),]),
             aes(x = TSSEnrichment,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution')

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),]),
             aes(x = TSSEnrichment,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution')

char <- sub(pattern = '-',replacement = '_',x = cell_type_1,fixed = TRUE)
char <- paste(char,'vs',sub(pattern = '-',replacement = '_',x = cell_type_2,fixed = TRUE),sep = '_')
char <- paste0('./res/step_40_fig_220623/',char,'_intra_species_TSSEnrichment_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

#nFrags distribution
macaque_multiome_ArchR$log_nFrags <- log10(macaque_multiome_ArchR$nFrags)

p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c(cell_type_1,cell_type_2),]),
             aes(x = log_nFrags,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution') + xlab('log10(nFrags)')

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),]),
             aes(x = log_nFrags,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution') + xlab('log10(nFrags)')

char <- sub(pattern = '-',replacement = '_',x = cell_type_1,fixed = TRUE)
char <- paste(char,'vs',sub(pattern = '-',replacement = '_',x = cell_type_2,fixed = TRUE),sep = '_')
char <- paste0('./res/step_40_fig_220623/',char,'_intra_species_nFrags_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

## Ex-2 vs Ex-3 ----------------------------------------------------------
cell_type_1 <- 'Ex-2'
cell_type_2 <- 'Ex-3'

matchobj <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                         groups = macaque_multiome_ArchR$cell_type,
                                         useGroups = cell_type_1,
                                         bgdGroups = cell_type_2,
                                         bias = c("TSSEnrichment", "log10(nFrags)"),
                                         k = 100,n = 1000,seed = 1,bufferRatio = 0.8)

#TSSenrichment distribution
p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c(cell_type_1,cell_type_2),]),
             aes(x = TSSEnrichment,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution')

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),]),
             aes(x = TSSEnrichment,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution')

char <- sub(pattern = '-',replacement = '_',x = cell_type_1,fixed = TRUE)
char <- paste(char,'vs',sub(pattern = '-',replacement = '_',x = cell_type_2,fixed = TRUE),sep = '_')
char <- paste0('./res/step_40_fig_220623/',char,'_intra_species_TSSEnrichment_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

#nFrags distribution
macaque_multiome_ArchR$log_nFrags <- log10(macaque_multiome_ArchR$nFrags)

p1 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c(cell_type_1,cell_type_2),]),
             aes(x = log_nFrags,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution') + xlab('log10(nFrags)')

p2 <- ggplot(data = as.data.frame(macaque_multiome_ArchR@cellColData[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),]),
             aes(x = log_nFrags,color = cell_type)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution') + xlab('log10(nFrags)')

char <- sub(pattern = '-',replacement = '_',x = cell_type_1,fixed = TRUE)
char <- paste(char,'vs',sub(pattern = '-',replacement = '_',x = cell_type_2,fixed = TRUE),sep = '_')
char <- paste0('./res/step_40_fig_220623/',char,'_intra_species_nFrags_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

# inter-species comparison ------------------------------------------------

## RG ----------------------------------------------------------------------
cell_type <- 'RG'

meta_data <- rbind(macaque_peak_matrix@meta.data[macaque_peak_matrix$cell_type == cell_type,c("TSSEnrichment","nFrags","cell_type","species")],
                   human_peak_matrix@meta.data[human_peak_matrix$cell_type == cell_type,c("TSSEnrichment","nFrags","cell_type","species")])

matchobj <- ArchR:::.matchBiasCellGroups(input = meta_data,
                                         groups = meta_data$species,
                                         useGroups = 'macaque',
                                         bgdGroups = 'human',
                                         bias = c("TSSEnrichment", "log10(nFrags)"),
                                         k = 100,n = 1000,seed = 1,bufferRatio = 0.8)

#TSSenrichment distribution
p1 <- ggplot(data = meta_data,
             aes(x = TSSEnrichment,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution')

p2 <- ggplot(data = meta_data[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),],
             aes(x = TSSEnrichment,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_40_fig_220623/',char,'_macaque_human_TSSEnrichment_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

#nFrags distribution
meta_data$log_nFrags <- log10(meta_data$nFrags)

p1 <- ggplot(data = meta_data,
             aes(x = log_nFrags,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution') + xlab('log10(nFrags)')

p2 <- ggplot(data = meta_data[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),],
             aes(x = log_nFrags,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution') + xlab('log10(nFrags)')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_40_fig_220623/',char,'_macaque_human_nFrags_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

## Ex-1 ----------------------------------------------------------------------
cell_type <- 'Ex-1'

meta_data <- rbind(macaque_peak_matrix@meta.data[macaque_peak_matrix$cell_type == cell_type,c("TSSEnrichment","nFrags","cell_type","species")],
                   human_peak_matrix@meta.data[human_peak_matrix$cell_type == cell_type,c("TSSEnrichment","nFrags","cell_type","species")])

matchobj <- ArchR:::.matchBiasCellGroups(input = meta_data,
                                         groups = meta_data$species,
                                         useGroups = 'macaque',
                                         bgdGroups = 'human',
                                         bias = c("TSSEnrichment", "log10(nFrags)"),
                                         k = 100,n = 1000,seed = 1,bufferRatio = 0.8)

#TSSenrichment distribution
p1 <- ggplot(data = meta_data,
             aes(x = TSSEnrichment,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution')

p2 <- ggplot(data = meta_data[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),],
             aes(x = TSSEnrichment,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_40_fig_220623/',char,'_macaque_human_TSSEnrichment_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

#nFrags distribution
meta_data$log_nFrags <- log10(meta_data$nFrags)

p1 <- ggplot(data = meta_data,
             aes(x = log_nFrags,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution') + xlab('log10(nFrags)')

p2 <- ggplot(data = meta_data[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),],
             aes(x = log_nFrags,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution') + xlab('log10(nFrags)')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_40_fig_220623/',char,'_macaque_human_nFrags_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

## Ex-3 ----------------------------------------------------------------------
cell_type <- 'Ex-3'

meta_data <- rbind(macaque_peak_matrix@meta.data[macaque_peak_matrix$cell_type == cell_type,c("TSSEnrichment","nFrags","cell_type","species")],
                   human_peak_matrix@meta.data[human_peak_matrix$cell_type == cell_type,c("TSSEnrichment","nFrags","cell_type","species")])

matchobj <- ArchR:::.matchBiasCellGroups(input = meta_data,
                                         groups = meta_data$species,
                                         useGroups = 'macaque',
                                         bgdGroups = 'human',
                                         bias = c("TSSEnrichment", "log10(nFrags)"),
                                         k = 100,n = 1000,seed = 1,bufferRatio = 0.8)

#TSSenrichment distribution
p1 <- ggplot(data = meta_data,
             aes(x = TSSEnrichment,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution')

p2 <- ggplot(data = meta_data[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),],
             aes(x = TSSEnrichment,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_40_fig_220623/',char,'_macaque_human_TSSEnrichment_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()

#nFrags distribution
meta_data$log_nFrags <- log10(meta_data$nFrags)

p1 <- ggplot(data = meta_data,
             aes(x = log_nFrags,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'raw distribution') + xlab('log10(nFrags)')

p2 <- ggplot(data = meta_data[c(matchobj$matchbgd[[1]]$cells,matchobj$matchbgd[[1]]$bgd),],
             aes(x = log_nFrags,color = species)) + 
  geom_density() + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.5,
        axis.title = element_text(size = 14,face = 'bold'),
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        legend.position = c(0.9,0.9)) + 
  labs(title = 'matched distribution') + xlab('log10(nFrags)')

char <- sub(pattern = '-',replacement = '_',x = cell_type,fixed = TRUE)
char <- paste0('./res/step_40_fig_220623/',char,'_macaque_human_nFrags_distribution.pdf')
pdf(file = char,width = 6,height = 6)
p1+p2+plot_layout(ncol = 1)
dev.off()