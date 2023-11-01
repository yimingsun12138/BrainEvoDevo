#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: species cellchat analysis focused on ANXA2                      ##
## Data: 2023.08.31                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0314')

# load cell-cell communication database -----------------------------------
Interaction_table <- read.table(file = './res/step_132_fig_230831/receptor_ligand_interactions_mitab_v1.0_April2017.txt',header = TRUE,sep = '\t',quote = '')
ligand_table <- read.table(file = './res/step_132_fig_230831/ligand.txt',header = TRUE,sep = '\t',quote = '')
receptor_table <- read.table(file = './res/step_132_fig_230831/receptor.txt',header = TRUE,sep = '\t',quote = '')
ecm_table <- read.table(file = './res/step_132_fig_230831/ECM.txt',header = TRUE,sep = '\t',quote = '')

#check data
'LGALS3' %in% ligand_table$Hgnc.Symbol
'LGALS3' %in% receptor_table$Hgnc.Symbol
'LGALS3' %in% ecm_table$Hgnc.Symbol

'LGALS3BP' %in% ligand_table$Hgnc.Symbol
'LGALS3BP' %in% receptor_table$Hgnc.Symbol
'LGALS3BP' %in% ecm_table$Hgnc.Symbol

temp <- Interaction_table[which(Interaction_table$AliasA == 'LGALS3' | Interaction_table$AliasB == 'LGALS3'),]

#check data
'PDGFD' %in% ligand_table$Hgnc.Symbol
'PDGFD' %in% receptor_table$Hgnc.Symbol
'PDGFD' %in% ecm_table$Hgnc.Symbol

temp <- Interaction_table[which(Interaction_table$AliasA == 'PDGFD' | Interaction_table$AliasB == 'PDGFD'),]

#check data
'ANXA2' %in% ligand_table$Hgnc.Symbol
'ANXA2' %in% receptor_table$Hgnc.Symbol
'ANXA2' %in% ecm_table$Hgnc.Symbol

temp <- Interaction_table[which(Interaction_table$AliasA == 'ANXA2' | Interaction_table$AliasB == 'ANXA2'),]
temp <- temp[which(temp$AliasA %in% ligand_table$Hgnc.Symbol | temp$AliasA %in% ecm_table$Hgnc.Symbol),]
temp <- temp[which(temp$AliasB %in% receptor_table$Hgnc.Symbol),]
table(temp$interactionType)
temp <- temp[which(temp$interactionType %in% c('interacts-with','psi-mi:MI:0407(direct interaction)')),]
