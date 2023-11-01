#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: human chimp macaque RG_neu homer motif enrichment               ##
## Data: 2023.07.10                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.1/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(dplyr)
library(Seurat)
library(ggplot2)
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
library(BSgenome.Mmusculus.UCSC.mm10)
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
library(riverplot)
library(CellChat)
library(OpenAI4R)
library(paletteer)
library(ggpattern)
library(ggrepel)
library(openxlsx)
library(readxl)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0613')

# human specific RG-neu DAP convert to chimp ------------------------------
#load data
consensus_peak <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')
RG_neu_DAP <- readRDS(file = './res/step_118_fig_230422/Wilcox/RG.1_df_wilcox.rds')
RG_neu_DAP <- RG_neu_DAP$peak[RG_neu_DAP$group == 'HGARs']
RG_neu_DAP <- unique(RG_neu_DAP)
table(RG_neu_DAP %in% consensus_peak$human_coord)

#get chimp peak coordinate
chr_name <- names(rtracklayer::import.chain(con = './data/reference/UCSC_chain_file_for_liftOver/panTro6ToHg38.over.chain'))
chimp_DAP <- my_unique_peakset_liftover(ori_GRanges = as(object = RG_neu_DAP,Class = 'GRanges'),
                                        UCSC_liftOver_path = '/content/data/sunym/software/UCSC_liftOver/liftOver',
                                        chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToPanTro6.over.chain',
                                        liftOver_mismatch = 0.1,length_filter = TRUE,length_mismatch = 0.1,chr_filter = TRUE,mapped_chr = chr_name,
                                        overlap_filter = TRUE,tmp_path = '/home/sunym/temp')
RG_neu_DAP <- as.character(chimp_DAP$ori)
names(RG_neu_DAP) <- NULL
table(RG_neu_DAP %in% consensus_peak$human_coord)
chimp_DAP <- chimp_DAP$mapped
mcols(chimp_DAP) <- NULL
rtracklayer::export.bed(object = chimp_DAP,con = './res/step_124_fig_230710/chimp_RG_neu_HGARs.bed')

#get human peak coordinate
human_DAP <- as(object = RG_neu_DAP,Class = 'GRanges')
rtracklayer::export.bed(object = human_DAP,con = './res/step_124_fig_230710/human_RG_neu_HGARs.bed')

#get macaque peak coordinate
macaque_DAP <- consensus_peak$macaque_coord[which(consensus_peak$human_coord %in% RG_neu_DAP)]
macaque_DAP <- as(object = macaque_DAP,Class = 'GRanges')
rtracklayer::export(object = macaque_DAP,con = './res/step_124_fig_230710/macaque_RG_neu_HGARs.bed')

# homer -------------------------------------------------------------------
#human
system('nohup findMotifsGenome.pl /content/data/sunym/project/Brain/res/step_124_fig_230710/human_RG_neu_HGARs.bed hg38 /home/sunym/temp/homer_res/human_RG_neu_HGARs_wilcox_homer -size given -mask > /home/sunym/temp/human.log 2>&1 &')

#macaque
system('nohup findMotifsGenome.pl /content/data/sunym/project/Brain/res/step_124_fig_230710/macaque_RG_neu_HGARs.bed /home/sunym/temp_genome/rheMac10.fa /home/sunym/temp/homer_res/macaque_RG_neu_HGARs_wilcox_homer -size given -mask > /home/sunym/temp/macaque.log 2>&1 &')

#chimp
system('nohup findMotifsGenome.pl /content/data/sunym/project/Brain/res/step_124_fig_230710/chimp_RG_neu_HGARs.bed panTro6 /home/sunym/temp/homer_res/chimp_RG_neu_HGARs_wilcox_homer -size given -mask > /home/sunym/temp/chimp.log 2>&1 &')