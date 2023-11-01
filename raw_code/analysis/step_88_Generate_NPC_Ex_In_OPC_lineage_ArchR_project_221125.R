#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Generate NPC Ex In OPC lineage ArchR project                    ##
## Data: 2022.11.25                                                                ##
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
library(ChIPseeker)
library(org.Hs.eg.db)
library(org.Mmu.eg.db)
library(DESeq2)
library(universalmotif)
library(topGO)
library(future.apply)
library(transPlotR)
library(aplot)
library(Biostrings)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# filter genes regulated by human SNC -------------------------------------

#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

affinity_difference <- readRDS(file = './res/step_88_fig_221124/affinity_difference.rds')

#load p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#modify affinity difference
affinity_difference$seq_name <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))

#cell type list
cell_type_list <- c('NPC','Ex','In','OPC','others')
candidate_gene <- c()

for (i in cell_type_list) {
  
  #cell type
  cell_type <- i
  
  #load gene list
  gene_list <- list.files(path = './res/step_88_fig_221124/gene_regulated_by_human_specific_SNC')
  gene_list <- gene_list[grep(pattern = cell_type,x = gene_list,fixed = TRUE)]
  gene_list <- paste('./res/step_88_fig_221124/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
  gene_list <- readRDS(file = gene_list)
  
  #load DF list
  DF_list <- list.files(path = './res/step_88_fig_221124/DF_list')
  DF_list <- DF_list[grep(pattern = cell_type,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_88_fig_221124/DF_list',DF_list,sep = '/')
  DF_list <- readRDS(DF_list)
  
  DF_list <- DF_list[DF_list$group == 'human_specific',"human_SNC"]
  
  #for loop
  for (j in gene_list) {
    temp <- which(p2g@metadata$geneSet$name == j)
    temp <- which(p2g$idxRNA %in% temp)
    temp <- p2g$idxATAC[temp]
    peak_list <- p2g@metadata$peakSet[temp]
    names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')
    
    #human specific SNC intersect with peak list
    SNC_list <- human_SNC
    names(SNC_list) <- human_SNC$human_SNC
    SNC_list <- countOverlaps(query = SNC_list[DF_list],subject = peak_list)
    SNC_list <- names(SNC_list)[SNC_list > 0]
    
    #if SNC list in affinity difference?
    if(sum(SNC_list %in% affinity_difference$seq_name) >= 1){
      candidate_gene <- c(candidate_gene,j)
    }
  }
}

candidate_gene <- unique(candidate_gene)

#save gene list
saveRDS(object = candidate_gene,file = './res/step_88_fig_221125/candidate_gene_list.rds')