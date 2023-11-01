#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Generate NPC Ex In OPC lineage ArchR project                    ##
## Data: 2022.11.24                                                                ##
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

# change work space -------------------------------------------------------
setwd('/home/sunym/temp')

# Greenleaf ATAC ArchR coverage ----------------------------------------------------

# copy Greenleaf_ATAC_ArchR
temp <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

#load Arrow file
ArrowFiles <- getArrowFiles(ArchRProj = temp)
file.exists(ArrowFiles)

#create ArchR project
Greenleaf_ATAC_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/home/sunym/temp/Greenleaf_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

#add meta data
Greenleaf_ATAC_ArchR <- Greenleaf_ATAC_ArchR[rownames(temp@cellColData)]
Greenleaf_ATAC_ArchR$cell_type <- temp$cell_type

Greenleaf_ATAC_ArchR$bulk <- 'others'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('RG-1','IP'),"bulk"] <- 'NPC'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"bulk"] <- 'Ex'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('RG-2','OPC'),"bulk"] <- 'OPC'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('InCGE','InMGE'),"bulk"] <- 'In'
table(Greenleaf_ATAC_ArchR$bulk)

#add peakset
Greenleaf_ATAC_ArchR <- addPeakSet(ArchRProj = Greenleaf_ATAC_ArchR,peakSet = getPeakSet(ArchRProj = temp))
Greenleaf_ATAC_ArchR <- addPeakMatrix(ArchRProj = Greenleaf_ATAC_ArchR)

#export bigwid
Greenleaf_ATAC_ArchR$ReadsInPromoter <- Greenleaf_ATAC_ArchR$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'bulk',normMethod = 'ReadsInPromoter',tileSize = 25,maxCells = NULL,verbose = TRUE)
ii <- gsub(pattern = 'ReadsInPromoter',replacement = 'ReadsInPeaks',fixed = TRUE,x = coverage_file)
for (i in 1:length(coverage_file)) {
  char <- paste('mv',coverage_file[i],ii[i],sep = ' ')
  system(char)
}

#notice point:
#maxpeaks = 150k
#macs2 FDR = 0.05
#modify mincells and max cells while doing pseudobulk
#mincells: 50

#add pseudo bulk
temp <- paste(Greenleaf_ATAC_ArchR$Sample,Greenleaf_ATAC_ArchR$bulk,sep = '_')
table(temp) %>% max()
table(temp) %>% min()
Greenleaf_ATAC_ArchR <- addGroupCoverages(ArchRProj = Greenleaf_ATAC_ArchR, 
                                          groupBy = "bulk", 
                                          minCells = 50, 
                                          maxCells = max(table(temp)), 
                                          maxFragments = 10^10, 
                                          minReplicates = 3, 
                                          maxReplicates = 3, 
                                          sampleRatio = 0.8)

temp <- as.data.frame(Greenleaf_ATAC_ArchR@projectMetadata$GroupCoverages$bulk$coverageMetadata)

#call peak
pathToMacs2 <- '/home/sunym/env/MACS2/bin/macs2'

Greenleaf_ATAC_ArchR <- addReproduciblePeakSet(
  ArchRProj = Greenleaf_ATAC_ArchR, 
  groupBy = "bulk", 
  pathToMacs2 = pathToMacs2,
  maxPeaks = 150000,
  genomeSize = 2.7e+09,
  cutOff = 0.05,
  force = TRUE
)

my_send_sms('call peak done!')

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#save data
peak_set <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
saveRDS(object = peak_set,file = '/content/data/sunym/project/Brain/res/step_88_fig_221124/human_merged_peakset.rds')

# macaque multiome ArchR coverage --------------------------------------------------

# copy macaque multiome ArchR
temp <- loadArchRProject(path = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_ArchR_221011/')

#load Arrow file
ArrowFiles <- getArrowFiles(ArchRProj = temp)
file.exists(ArrowFiles)

#create ArchR project
macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = '/home/sunym/temp/macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = temp),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = temp)
)

#add meta data
macaque_multiome_ArchR <- macaque_multiome_ArchR[rownames(temp@cellColData)]
macaque_multiome_ArchR$cell_type <- temp$cell_type
macaque_multiome_ArchR <- macaque_multiome_ArchR[macaque_multiome_ArchR$cell_type != 'Cycling']

macaque_multiome_ArchR$bulk <- 'others'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('RG-1','IP'),"bulk"] <- 'NPC'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"bulk"] <- 'Ex'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('RG-2','OPC'),"bulk"] <- 'OPC'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('InCGE','InMGE'),"bulk"] <- 'In'
table(macaque_multiome_ArchR$bulk)

#add peakset
macaque_multiome_ArchR <- addPeakSet(ArchRProj = macaque_multiome_ArchR,peakSet = getPeakSet(ArchRProj = temp))
macaque_multiome_ArchR <- addPeakMatrix(ArchRProj = macaque_multiome_ArchR)

#export bigwid
macaque_multiome_ArchR$ReadsInPromoter <- macaque_multiome_ArchR$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = macaque_multiome_ArchR,groupBy = 'bulk',normMethod = 'ReadsInPromoter',tileSize = 25,maxCells = NULL,verbose = TRUE)
ii <- gsub(pattern = 'ReadsInPromoter',replacement = 'ReadsInPeaks',fixed = TRUE,x = coverage_file)
for (i in 1:length(coverage_file)) {
  char <- paste('mv',coverage_file[i],ii[i],sep = ' ')
  system(char)
}

#notice point:
#maxpeaks = 150k
#macs2 FDR = 0.05
#modify mincells and max cells while doing pseudobulk
#mincells: 200

#add pseudo bulk
temp <- paste(macaque_multiome_ArchR$Sample,macaque_multiome_ArchR$bulk,sep = '_')
table(temp) %>% max()
table(temp) %>% min()
macaque_multiome_ArchR <- addGroupCoverages(ArchRProj = macaque_multiome_ArchR, 
                                            groupBy = "bulk", 
                                            minCells = 200, 
                                            maxCells = max(table(temp)), 
                                            maxFragments = 10^10, 
                                            minReplicates = 4, 
                                            maxReplicates = 4, 
                                            sampleRatio = 0.8)

temp <- as.data.frame(macaque_multiome_ArchR@projectMetadata$GroupCoverages$bulk$coverageMetadata)

#call peak
pathToMacs2 <- '/home/sunym/env/MACS2/bin/macs2'

macaque_multiome_ArchR <- addReproduciblePeakSet(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "bulk", 
  pathToMacs2 = pathToMacs2,
  maxPeaks = 150000,
  genomeSize = 2.7e+09,
  cutOff = 0.05,
  force = TRUE
)

my_send_sms('call peak done!')

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#save data
peak_set <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
saveRDS(object = peak_set,file = '/content/data/sunym/project/Brain/res/step_88_fig_221124/macaque_merged_peakset.rds')

# DF human SNC list -------------------------------------------------------

#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

human_peak_matrix <- readRDS(file = './res/step_78_fig_221111/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_78_fig_221111/macaque_peak_matrix.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#add meta data
Greenleaf_ATAC_ArchR$bulk <- 'others'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('RG-1','IP'),"bulk"] <- 'NPC'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"bulk"] <- 'Ex'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('RG-2','OPC'),"bulk"] <- 'OPC'
Greenleaf_ATAC_ArchR@cellColData[Greenleaf_ATAC_ArchR$cell_type %in% c('InCGE','InMGE'),"bulk"] <- 'In'
table(Greenleaf_ATAC_ArchR$bulk)

macaque_multiome_ArchR$bulk <- 'others'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('RG-1','IP'),"bulk"] <- 'NPC'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('Ex-1','Ex-2','Ex-3','Ex-4'),"bulk"] <- 'Ex'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('RG-2','OPC'),"bulk"] <- 'OPC'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('InCGE','InMGE'),"bulk"] <- 'In'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$cell_type %in% c('Cycling'),"bulk"] <- 'Cycling'
table(macaque_multiome_ArchR$bulk)

#generate cell type list
cell_type_list <- c('NPC','Ex','In','OPC','others')

#for loop
for (i in cell_type_list) {
  
  #set cell type
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #subset peak matrix cells
  human_cell_list <- rownames(Greenleaf_ATAC_ArchR@cellColData)[c(Greenleaf_ATAC_ArchR$bulk == cell_type)]
  macaque_cell_list <- rownames(macaque_multiome_ArchR@cellColData)[c(macaque_multiome_ArchR$bulk == cell_type)]
  
  subset_human_peak_matrix <- human_peak_matrix[,(human_cell_list)]
  subset_macaque_peak_matrix <- macaque_peak_matrix[,c(macaque_cell_list)]
  
  #generate normalize factor
  norm_factor <- c(Greenleaf_ATAC_ArchR@cellColData[c(human_cell_list),'ReadsInPeaks'],macaque_multiome_ArchR@cellColData[c(macaque_cell_list),'ReadsInPeaks'])
  norm_factor <- median(norm_factor) * 10^4
  
  #normalize
  subset_human_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = human_cell_list,FUN = function(x){
    temp <- subset_human_peak_matrix[,x]/Greenleaf_ATAC_ArchR@cellColData[x,'ReadsInPeaks']*norm_factor
    return(temp)
  }))
  colnames(subset_human_peak_matrix) <- human_cell_list
  rownames(subset_human_peak_matrix) <- rownames(human_peak_matrix)
  
  subset_macaque_peak_matrix <- base::do.call(what = cbind,args = base::lapply(X = macaque_cell_list,FUN = function(x){
    temp <- subset_macaque_peak_matrix[,x]/macaque_multiome_ArchR@cellColData[x,'ReadsInPeaks']*norm_factor
    return(temp)
  }))
  colnames(subset_macaque_peak_matrix) <- macaque_cell_list
  rownames(subset_macaque_peak_matrix) <- rownames(macaque_peak_matrix)
  
  #get peak list
  
  #human
  peak_file <- list.files(path = './res/step_88_fig_221124/human_peak/PeakCalls')
  peak_file <- peak_file[c(grep(pattern = cell_type_dot,x = peak_file,fixed = TRUE))]
  peak_file <- paste('./res/step_88_fig_221124/human_peak/PeakCalls',peak_file,sep = '/')
  file.exists(peak_file) %>% print()
  peak_file <- readRDS(file = peak_file)
  
  names(human_SNC) <- paste(human_SNC@seqnames,as.character(human_SNC@ranges),sep = '-')
  temp <- countOverlaps(query = human_SNC,subject = peak_file)
  temp <- names(temp)[temp > 0]
  temp <- human_SNC[c(temp)]$human_SNC
  
  peak_list <- temp
  
  #macaque
  peak_file <- list.files(path = './res/step_88_fig_221124/macaque_peak/PeakCalls')
  peak_file <- peak_file[c(grep(pattern = cell_type_dot,x = peak_file,fixed = TRUE))]
  peak_file <- paste('./res/step_88_fig_221124/macaque_peak/PeakCalls',peak_file,sep = '/')
  file.exists(peak_file) %>% print()
  peak_file <- readRDS(file = peak_file)
  
  names(macaque_SNC) <- paste(macaque_SNC@seqnames,as.character(macaque_SNC@ranges),sep = '-')
  temp <- countOverlaps(query = macaque_SNC,subject = peak_file)
  temp <- names(temp)[temp > 0]
  temp <- macaque_SNC[c(temp)]$human_SNC
  
  peak_list <- c(peak_list,temp)
  
  #unify
  peak_list <- unique(peak_list)
  
  #subset peak matrix peaks
  
  #human
  names(human_SNC) <- paste(human_SNC@seqnames,as.character(human_SNC@ranges),sep = '-')
  (rownames(subset_human_peak_matrix) %in% names(human_SNC)) %>% table()
  rownames(subset_human_peak_matrix) <- human_SNC[rownames(subset_human_peak_matrix)]$human_SNC
  subset_human_peak_matrix <- subset_human_peak_matrix[c(peak_list),]
  
  #macaque
  names(macaque_SNC) <- paste(macaque_SNC@seqnames,as.character(macaque_SNC@ranges),sep = '-')
  (rownames(subset_macaque_peak_matrix) %in% names(macaque_SNC)) %>% table()
  rownames(subset_macaque_peak_matrix) <- macaque_SNC[rownames(subset_macaque_peak_matrix)]$human_SNC
  subset_macaque_peak_matrix <- subset_macaque_peak_matrix[c(peak_list),]
  
  #run wilcox
  DF_list <- my_DF_wilcox_test(mat1 = subset_human_peak_matrix,mat2 = subset_macaque_peak_matrix,
                               alternative = 'two.sided',paired = FALSE,workers = 1,
                               future.globals.maxSize = 200*(1024^3))
  DF_list$human_SNC <- rownames(DF_list)
  
  #human specific
  human_specific_DF_list <- DF_list[c(DF_list$log2FC > 1 & DF_list$fdr < 0.01),]
  
  #intersect with human peaks
  peak_file <- list.files(path = './res/step_88_fig_221124/human_peak/PeakCalls')
  peak_file <- peak_file[c(grep(pattern = cell_type_dot,x = peak_file,fixed = TRUE))]
  peak_file <- paste('./res/step_88_fig_221124/human_peak/PeakCalls',peak_file,sep = '/')
  file.exists(peak_file) %>% print()
  peak_file <- readRDS(file = peak_file)
  
  names(human_SNC) <- human_SNC$human_SNC
  temp <- countOverlaps(query = human_SNC,subject = peak_file)
  temp <- names(temp)[temp > 0]
  
  human_specific_DF_list <- human_specific_DF_list[c(human_specific_DF_list$human_SNC %in% temp),]
  human_specific_DF_list$group <- 'human_specific'
  
  #macaque specific
  macaque_specific_DF_list <- DF_list[c(DF_list$log2FC < -1 & DF_list$fdr < 0.01),]
  
  #intersect with macaque peaks
  peak_file <- list.files(path = './res/step_88_fig_221124/macaque_peak/PeakCalls')
  peak_file <- peak_file[c(grep(pattern = cell_type_dot,x = peak_file,fixed = TRUE))]
  peak_file <- paste('./res/step_88_fig_221124/macaque_peak/PeakCalls',peak_file,sep = '/')
  file.exists(peak_file) %>% print()
  peak_file <- readRDS(file = peak_file)
  
  names(macaque_SNC) <- macaque_SNC$human_SNC
  temp <- countOverlaps(query = macaque_SNC,subject = peak_file)
  temp <- names(temp)[temp > 0]
  
  macaque_specific_DF_list <- macaque_specific_DF_list[c(macaque_specific_DF_list$human_SNC %in% temp),]
  macaque_specific_DF_list$group <- 'macaque_specific'
  
  #species conserved
  species_conserved_DF_list <- DF_list[c(abs(DF_list$log2FC) <= 1),]
  
  #intersect with human peaks
  peak_file <- list.files(path = './res/step_88_fig_221124/human_peak/PeakCalls')
  peak_file <- peak_file[c(grep(pattern = cell_type_dot,x = peak_file,fixed = TRUE))]
  peak_file <- paste('./res/step_88_fig_221124/human_peak/PeakCalls',peak_file,sep = '/')
  file.exists(peak_file) %>% print()
  peak_file <- readRDS(file = peak_file)
  
  names(human_SNC) <- human_SNC$human_SNC
  temp <- countOverlaps(query = human_SNC,subject = peak_file)
  temp <- names(temp)[temp > 0]
  
  species_conserved_DF_list <- species_conserved_DF_list[c(species_conserved_DF_list$human_SNC %in% temp),]
  
  #intersect with macaque peaks
  peak_file <- list.files(path = './res/step_88_fig_221124/macaque_peak/PeakCalls')
  peak_file <- peak_file[c(grep(pattern = cell_type_dot,x = peak_file,fixed = TRUE))]
  peak_file <- paste('./res/step_88_fig_221124/macaque_peak/PeakCalls',peak_file,sep = '/')
  file.exists(peak_file) %>% print()
  peak_file <- readRDS(file = peak_file)
  
  names(macaque_SNC) <- macaque_SNC$human_SNC
  temp <- countOverlaps(query = macaque_SNC,subject = peak_file)
  temp <- names(temp)[temp > 0]
  
  species_conserved_DF_list <- species_conserved_DF_list[c(species_conserved_DF_list$human_SNC %in% temp),]
  species_conserved_DF_list$group <- 'species_conserved'
  
  #save data
  DF_list <- rbind(human_specific_DF_list,macaque_specific_DF_list,species_conserved_DF_list)
  char <- paste0('./res/step_88_fig_221124/DF_list/',cell_type_dot,'_DF_list.rds')
  saveRDS(object = DF_list,file = char)
  
  #report
  char <- paste(cell_type,'done!',sep = ' ')
  print(char)
}

# gene regulated by human specific SNC ------------------------------------

#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

cell_type_list <- c('NPC','Ex','In','OPC','others')

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#for loop
for (i in cell_type_list) {
  cell_type <- i
  
  #load DF list
  DF_list <- list.files(path = './res/step_88_fig_221124/DF_list')
  DF_list <- DF_list[grep(pattern = cell_type,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_88_fig_221124/DF_list',DF_list,sep = '/')
  DF_list <- readRDS(DF_list)
  
  #filter DF list
  DF_list <- DF_list[DF_list$group == 'human_specific',"human_SNC"]
  names(human_SNC) <- human_SNC$human_SNC
  DF_list <- human_SNC[DF_list]
  
  #get peak list overlap with DF list
  peak_list <- p2g@metadata$peakSet
  names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')
  
  DF_list <- countOverlaps(query = peak_list,subject = DF_list)
  DF_list <- names(DF_list)[DF_list > 0]
  
  #get ATAC idx
  idx_ATAC <- which(names(peak_list) %in% DF_list)
  temp <- which(p2g$idxATAC %in% idx_ATAC)
  idx_RNA <- p2g$idxRNA[temp]
  
  #get gene regulated by human specific SNC
  gene_list <- p2g@metadata$geneSet$name[idx_RNA]
  gene_list <- unique(gene_list)
  
  #save gene list
  char <- paste0('./res/step_88_fig_221124/gene_regulated_by_human_specific_SNC/',cell_type,'_gene_regulated_by_human_specific_SNC.rds')
  saveRDS(object = gene_list,file = char)
  char <- paste(cell_type,'done!',sep = ' ')
  print(char)
}

# generate affinity difference --------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get fimo out
human_fimo_out <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo_out <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)

human_fimo_out$idx <- paste(human_fimo_out$motif_id,human_fimo_out$sequence_name,human_fimo_out$strand,sep = '#')
macaque_fimo_out$idx <- paste(macaque_fimo_out$motif_id,macaque_fimo_out$sequence_name,macaque_fimo_out$strand,sep = '#')

#back up fimo out
human_fimo_out_backup <- human_fimo_out
macaque_fimo_out_backup <- macaque_fimo_out

#filter fimo out
human_fimo_out <- human_fimo_out[human_fimo_out$q.value < 0.05,]
macaque_fimo_out <- macaque_fimo_out[macaque_fimo_out$q.value < 0.05,]

idx <- which(human_fimo_out$start <= 31 & human_fimo_out$stop >= 31)
human_fimo_out <- human_fimo_out[idx,]

idx <- which(macaque_fimo_out$start <= 31 & macaque_fimo_out$stop >= 31)
macaque_fimo_out <- macaque_fimo_out[idx,]

#get highest motif-sequence-strand combination
human_fimo_out <- human_fimo_out[order(human_fimo_out$score,decreasing = TRUE),]
temp <- human_fimo_out$idx
temp <- which(!duplicated(temp))
human_fimo_out <- human_fimo_out[temp,]

macaque_fimo_out <- macaque_fimo_out[order(macaque_fimo_out$score,decreasing = TRUE),]
temp <- macaque_fimo_out$idx
temp <- which(!duplicated(temp))
macaque_fimo_out <- macaque_fimo_out[temp,]

#get full set
full_set <- unique(c(human_fimo_out$idx,macaque_fimo_out$idx))
affinity_difference <- data.frame(set = full_set,human_score = NA,macaque_score = NA)

#get human score
affinity_difference$human_score <- unlist(base::lapply(X = full_set,FUN = function(x){
  if(x %in% human_fimo_out_backup$idx){
    temp <- human_fimo_out_backup[human_fimo_out_backup$idx == x,"score"]
    return(max(temp))
  }else{
    temp <- strsplit(x = x,split = '#')
    temp <- temp[[1]][1]
    if(temp %in% human_fimo_out_backup$motif_id){
      temp <- human_fimo_out_backup[human_fimo_out_backup$motif_id == temp,"score"]
      return(min(temp))
    }else{
      return(NA)
    }
  }
}))

affinity_difference$macaque_score <- unlist(base::lapply(X = full_set,FUN = function(x){
  if(x %in% macaque_fimo_out_backup$idx){
    temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$idx == x,"score"]
    return(max(temp))
  }else{
    temp <- strsplit(x = x,split = '#')
    temp <- temp[[1]][1]
    if(temp %in% macaque_fimo_out_backup$motif_id){
      temp <- macaque_fimo_out_backup[macaque_fimo_out_backup$motif_id == temp,"score"]
      return(min(temp))
    }else{
      return(NA)
    }
  }
}))

affinity_difference <- na.omit(object = affinity_difference)
affinity_difference$difference <- affinity_difference$human_score - affinity_difference$macaque_score

#save data
saveRDS(object = affinity_difference,file = './res/step_88_fig_221124/affinity_difference.rds')
