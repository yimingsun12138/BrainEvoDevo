#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: generate GO enrichment table for supplementary                  ##
## Data: 2023.05.29                                                                ##
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

# load data ---------------------------------------------------------------
Human_specific_signature <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
Primate_conserved_signature <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
Species_conserved_signature <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')

# human specific gene -----------------------------------------------------
gene_list <- Human_specific_signature
group_label <- 'human-specific RG-neu marker genes'
res_name <- 'HS_GO'

background_gene_list <- rownames(readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds'))
gene_list <- c(background_gene_list %in% gene_list)
names(gene_list) <- background_gene_list
gene_list[which(gene_list == TRUE)] <- '1'
gene_list[which(gene_list == FALSE)] <- '0'
gene_list <- factor(gene_list,levels = c('0','1'))

GO_Res <- base::do.call(what = rbind,args = base::lapply(X = c('BP','CC','MF'),FUN = function(ontology_source){
  GO_ontology <- ontology_source
  GO_enrich <- new("topGOdata",
                   description = paste(group_label,GO_ontology,sep = ' '),
                   ontology = GO_ontology,
                   allGenes = gene_list,
                   nodeSize = 10,annotationFun = annFUN.org,
                   mapping = 'org.Hs.eg.db',ID = 'symbol')
  resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GO_enrich, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 100, numChar = 10^3)
  
  #get significant genes
  Significant_genes <- base::lapply(X = 1:nrow(allRes),FUN = function(x){
    #get GO term genes
    GO_genes <- getGOgeneSet(x = allRes$GO.ID[x],OrgDb = 'org.Hs.eg.db',ont = GO_ontology,keytype = 'SYMBOL')
    GO_genes <- dplyr::intersect(x = names(gene_list)[gene_list == '1'],y = GO_genes)
    if(length(GO_genes) != as.numeric(allRes$Significant[x])){
      stop('Significant gene number wrong!')
    }
    return(paste(GO_genes,collapse = ','))
  })
  
  allRes$Ontology <- GO_ontology
  allRes$Significant_genes <- unlist(Significant_genes)
  allRes <- allRes[,c("GO.ID","Term","Ontology","Annotated","Significant","Significant_genes","classicFisher")]
  return(allRes)
}))
assign(x = res_name,value = GO_Res)

# primate conserved gene -----------------------------------------------------
gene_list <- Primate_conserved_signature
group_label <- 'primate-conserved RG-neu marker genes'
res_name <- 'PC_GO'

background_gene_list <- rownames(readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds'))
gene_list <- c(background_gene_list %in% gene_list)
names(gene_list) <- background_gene_list
gene_list[which(gene_list == TRUE)] <- '1'
gene_list[which(gene_list == FALSE)] <- '0'
gene_list <- factor(gene_list,levels = c('0','1'))

GO_Res <- base::do.call(what = rbind,args = base::lapply(X = c('BP','CC','MF'),FUN = function(ontology_source){
  GO_ontology <- ontology_source
  GO_enrich <- new("topGOdata",
                   description = paste(group_label,GO_ontology,sep = ' '),
                   ontology = GO_ontology,
                   allGenes = gene_list,
                   nodeSize = 10,annotationFun = annFUN.org,
                   mapping = 'org.Hs.eg.db',ID = 'symbol')
  resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GO_enrich, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 100, numChar = 10^3)
  
  #get significant genes
  Significant_genes <- base::lapply(X = 1:nrow(allRes),FUN = function(x){
    #get GO term genes
    GO_genes <- getGOgeneSet(x = allRes$GO.ID[x],OrgDb = 'org.Hs.eg.db',ont = GO_ontology,keytype = 'SYMBOL')
    GO_genes <- dplyr::intersect(x = names(gene_list)[gene_list == '1'],y = GO_genes)
    if(length(GO_genes) != as.numeric(allRes$Significant[x])){
      stop('Significant gene number wrong!')
    }
    return(paste(GO_genes,collapse = ','))
  })
  
  allRes$Ontology <- GO_ontology
  allRes$Significant_genes <- unlist(Significant_genes)
  allRes <- allRes[,c("GO.ID","Term","Ontology","Annotated","Significant","Significant_genes","classicFisher")]
  return(allRes)
}))
assign(x = res_name,value = GO_Res)

# species conserved gene -----------------------------------------------------
gene_list <- Species_conserved_signature
group_label <- 'species-conserved RG-neu marker genes'
res_name <- 'SC_GO'

background_gene_list <- rownames(readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds'))
gene_list <- c(background_gene_list %in% gene_list)
names(gene_list) <- background_gene_list
gene_list[which(gene_list == TRUE)] <- '1'
gene_list[which(gene_list == FALSE)] <- '0'
gene_list <- factor(gene_list,levels = c('0','1'))

GO_Res <- base::do.call(what = rbind,args = base::lapply(X = c('BP','CC','MF'),FUN = function(ontology_source){
  GO_ontology <- ontology_source
  GO_enrich <- new("topGOdata",
                   description = paste(group_label,GO_ontology,sep = ' '),
                   ontology = GO_ontology,
                   allGenes = gene_list,
                   nodeSize = 10,annotationFun = annFUN.org,
                   mapping = 'org.Hs.eg.db',ID = 'symbol')
  resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GO_enrich, classicFisher = resultFisher, orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 100, numChar = 10^3)
  
  #get significant genes
  Significant_genes <- base::lapply(X = 1:nrow(allRes),FUN = function(x){
    #get GO term genes
    GO_genes <- getGOgeneSet(x = allRes$GO.ID[x],OrgDb = 'org.Hs.eg.db',ont = GO_ontology,keytype = 'SYMBOL')
    GO_genes <- dplyr::intersect(x = names(gene_list)[gene_list == '1'],y = GO_genes)
    if(length(GO_genes) != as.numeric(allRes$Significant[x])){
      stop('Significant gene number wrong!')
    }
    return(paste(GO_genes,collapse = ','))
  })
  
  allRes$Ontology <- GO_ontology
  allRes$Significant_genes <- unlist(Significant_genes)
  allRes <- allRes[,c("GO.ID","Term","Ontology","Annotated","Significant","Significant_genes","classicFisher")]
  return(allRes)
}))
assign(x = res_name,value = GO_Res)

# combine data into excel -------------------------------------------------
all_GO_Res <- createWorkbook()
addWorksheet(wb = all_GO_Res,sheetName = 'human-specific marker')
addWorksheet(wb = all_GO_Res,sheetName = 'primate-conserved marker')
addWorksheet(wb = all_GO_Res,sheetName = 'species-conserved marker')

writeDataTable(wb = all_GO_Res,sheet = 'human-specific marker',x = HS_GO)
writeDataTable(wb = all_GO_Res,sheet = 'primate-conserved marker',x = PC_GO)
writeDataTable(wb = all_GO_Res,sheet = 'species-conserved marker',x = SC_GO)

saveWorkbook(wb = all_GO_Res,file = './res/step_119_fig_230529/RG_neu_marker_gene_GO_enrichment.xlsx',overwrite = TRUE)
