#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: reviewed genes affected by human specific SNC                   ##
## Data: 2022.11.27                                                                ##
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

# load gene list ----------------------------------------------------------
gene_list <- c('ANXA2','ADK','SFXN3','BIRC5','CSTB','ELP6','CSPG5','SST','HES1',
               'RTN4RL2','SNCB','UNC5A','STK33','DENND2B','TMEM9B','ATP1A3','BEGAIN',
               'FGF14','RBFOX2','MAPKAPK3','UCHL3','LMO7','IGSF9B','GRIP1','CAND1',
               'NXPH1','GFRA2','DMTN','LGI3')

# plot function -----------------------------------------------------------
plot_track <- function(gene_list = gene_list,species = 'human'){
  #load data
  human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
  macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')
  macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
  Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
  color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')
  human_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
  human_anno <- rtracklayer::as.data.frame(x = human_anno)
  macaque_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
  macaque_anno <- rtracklayer::as.data.frame(x = macaque_anno)
  
  affinity_difference <- readRDS(file = './res/step_88_fig_221124/affinity_difference.rds')
  affinity_difference$seq_name <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
    temp <- strsplit(x = x,split = '#')
    return(temp[[1]][2])
  }))
  
  #set parameter
  cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','InMGE','InCGE','OPC')
  cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)
  extend_size <- 100000
  resolution <- 1000
  ylim_track <- 0.03
  
  #load bw
  if(species == 'human'){
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type/')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }else if(species == 'macaque'){
    cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type')
    cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
      temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
      return(temp)
    }))
    cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/macaque/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
    cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
    cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
    cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
    gc()
  }else{
    stop('species error!')
  }
  
  #for loop
  for (i in gene_list) {
    #set gene name
    gene_name <- i
    
    #get p2g
    p2g <- getPeak2GeneLinks(
      ArchRProj = Greenleaf_ATAC_ArchR,
      corCutOff = 0.45,
      resolution = 1,
      returnLoops = FALSE
    )
    
    #get SNC list that may affect the gene express
    temp <- which(p2g@metadata$geneSet$name == gene_name)
    temp <- which(p2g$idxRNA %in% temp)
    temp <- p2g$idxATAC[temp]
    temp <- p2g@metadata$peakSet[temp]
    
    names(human_SNC) <- human_SNC$human_SNC
    names(macaque_SNC) <- macaque_SNC$human_SNC
    SNC_list <- countOverlaps(query = human_SNC,subject = temp)
    SNC_list <- names(SNC_list)[SNC_list > 0]
    
    temp <- affinity_difference
    temp <- temp[temp$seq_name %in% SNC_list,]
    SNC_list <- unique(temp$seq_name)
    if(species == 'human'){
      SNC_list <- human_SNC[SNC_list]
    }else if(species == 'macaque'){
      SNC_list <- macaque_SNC[SNC_list]
    }
    
    #get target site
    if(species == 'human'){
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = Greenleaf_ATAC_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }else if(species == 'macaque'){
      TSS_site <- getPeak2GeneLinks(
        ArchRProj = macaque_multiome_ArchR,
        corCutOff = 0.45,
        resolution = 1,
        returnLoops = FALSE
      )
    }else{
      stop('species error!')
    }
    
    TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]
    if(length(TSS_site) > 1){
      stop('TSS number error!')
    }
    chrom <- as.character(TSS_site@seqnames)
    start_site <- TSS_site@ranges@start - extend_size
    end_site <- TSS_site@ranges@start + extend_size
    
    #modify extension
    temp <- rtracklayer::as.data.frame(x = SNC_list)
    start_range <- min(temp$start,temp$end)
    end_range <- max(temp$start,temp$end)
    
    if(start_site >= start_range){
      start_site <- start_range - 10000
    }
    if(end_site <= end_range){
      end_site <- end_range + 10000
    }
    
    #plot gene annotation
    if(species == 'human'){
      p_gene_anno <- trancriptVis(gtfFile = human_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }else if(species == 'macaque'){
      p_gene_anno <- trancriptVis(gtfFile = macaque_anno,collapse = TRUE,gene = gene_name,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                                  posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',
                                  newStyleArrow = FALSE,xAxis.info = TRUE)
    }else{
      stop('species error!')
    }
    
    #plot peak
    if(species == 'human'){
      temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
      temp <- human_SNC
      names(temp) <- temp$human_SNC
      rtracklayer::export(object = temp,con = '/home/sunym/temp/SNC.bed',format = 'bed')
      rtracklayer::export(object = SNC_list,con = '/home/sunym/temp/SNC_list.bed',format = 'bed')
    }else if(species == 'macaque'){
      temp <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
      rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
      temp <- macaque_SNC
      names(temp) <- temp$human_SNC
      rtracklayer::export(object = temp,con = '/home/sunym/temp/SNC.bed',format = 'bed')
      rtracklayer::export(object = SNC_list,con = '/home/sunym/temp/SNC_list.bed',format = 'bed')
    }else{
      stop('species error!')
    }
    
    p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed','/home/sunym/temp/SNC.bed','/home/sunym/temp/SNC_list.bed'),
                          chr = chrom,region.min = start_site - 10000,region.max = end_site + 10000,
                          show.legend = FALSE,fill = ggsci::pal_d3()(3),track.width = 0.3)
    
    #plot track
    if(species == 'human'){
      temp_anno <- human_anno
    }else if(species == 'macaque'){
      temp_anno <- macaque_anno
    }else{
      stop('species error!')
    }
    
    p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = temp_anno,
                             chr = chrom,region.min = start_site,region.max = end_site,
                             sample.order = cell_type_list,space.y = 0,
                             y.max = ylim_track,color = color_param$celltype[cell_type_list],
                             theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)
    
    #plot link
    if(species == 'human'){
      p2g <- getPeak2GeneLinks(
        ArchRProj = Greenleaf_ATAC_ArchR,
        corCutOff = 0.45,
        resolution = resolution,
        returnLoops = FALSE
      )
    }else if(species == 'macaque'){
      p2g <- getPeak2GeneLinks(
        ArchRProj = macaque_multiome_ArchR,
        corCutOff = 0.45,
        resolution = resolution,
        returnLoops = FALSE
      )
    }else{
      stop('species error!')
    }
    
    temp <- which(p2g@metadata$geneSet$name == gene_name)
    temp <- which(p2g$idxRNA %in% temp)
    if(length(temp) == 0){
      links_meta_data <- data.frame(chr = chrom,
                                    start = c(start_site - 3000,start_site - 6000),
                                    end = c(end_site + 3000,end_site + 6000),
                                    cor = c(1,0),group = 'group1')
    }else{
      idx_ATAC <- p2g$idxATAC[temp]
      idx_ATAC <- p2g@metadata$peakSet[idx_ATAC]
      idx_ATAC <- rtracklayer::as.data.frame(idx_ATAC)
      idx_ATAC <- round((idx_ATAC$start + idx_ATAC$end) / 2)
      links_meta_data <- data.frame(chr = chrom,start = NA,end = NA,cor = p2g$Correlation[temp],group = 'group1')
      for (i in 1:length(idx_ATAC)) {
        if(TSS_site@ranges@start <= idx_ATAC[i]){
          links_meta_data[i,'start'] <- TSS_site@ranges@start
          links_meta_data[i,"end"] <- idx_ATAC[i]
        }else{
          links_meta_data[i,"start"] <- idx_ATAC[i]
          links_meta_data[i,'end'] <- TSS_site@ranges@start
        }
      }
      temp <- data.frame(chr = chrom,
                         start = c(start_site - 3000,start_site - 6000),
                         end = c(end_site + 3000,end_site + 6000),
                         cor = c(1,0),group = 'group1')
      links_meta_data <- rbind(links_meta_data,temp)
    }
    p_p2g <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                     facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),
                     xAixs.info = FALSE,curvature = 0.3)
    
    #save plot
    char <- paste0('./res/step_90_fig_221127/',species,'_',gene_name,'_trackplot.pdf')
    pdf(file = char,width = 8,height = 8)
    print(p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.05) %>% insert_bottom(plot = p_p2g,height = 0.1) %>% insert_bottom(plot = p_gene_anno,height = 0.05))
    dev.off()
  }
  return(NULL)
}

# track plot --------------------------------------------------------------
plot_track(gene_list = gene_list,species = 'human')
plot_track(gene_list = gene_list,species = 'macaque')
