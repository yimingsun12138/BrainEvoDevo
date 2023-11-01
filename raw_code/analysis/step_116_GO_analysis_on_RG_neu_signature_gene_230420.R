#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: GO analysis on RG_neu signature gene                            ##
## Data: 2023.04.20                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.2.3/')
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session()

# human-specific RG-neu signature gene ------------------------------------
#set gene list
category <- 'human-specific'
gene_list <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')

#pipeline
background_gene_list <- rownames(readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds'))
gene_list <- c(background_gene_list %in% gene_list)
names(gene_list) <- background_gene_list
gene_list[which(gene_list == TRUE)] <- '1'
gene_list[which(gene_list == FALSE)] <- '0'
gene_list <- factor(gene_list,levels = c('0','1'))

for (GO_ontology in c('BP','CC','MF')) {
  GO_enrich <- new("topGOdata",
                   description = paste(category,GO_ontology,sep = ' '),
                   ontology = GO_ontology,
                   allGenes = gene_list,
                   nodeSize = 10,annotationFun = annFUN.org,
                   mapping = 'org.Hs.eg.db',ID = 'symbol')
  resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GO_enrich, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = length(resultFisher@score), numChar = 10^3)
  allRes$classicFisher <- as.numeric(allRes$classicFisher)
  allRes[which(allRes$classicFisher < 10^-30),"classicFisher"] <- 10^-30
  allRes$sig <- -log10(allRes$classicFisher)
  assign(x = paste('allRes',GO_ontology,sep = '_'),value = allRes)
}

sig_BP <- allRes_BP$GO.ID[1:100]
sig_BP <- filter_child_GO_term(sig_BP,'BP')
sig_BP <- allRes_BP[allRes_BP$GO.ID %in% sig_BP,]
sig_CC <- allRes_CC$GO.ID[1:100]
sig_CC <- filter_child_GO_term(sig_CC,'CC')
sig_CC <- allRes_CC[allRes_CC$GO.ID %in% sig_CC,]
sig_MF <- allRes_MF$GO.ID[1:100]
sig_MF <- filter_child_GO_term(sig_MF,'MF')
sig_MF <- allRes_MF[allRes_MF$GO.ID %in% sig_MF,]
