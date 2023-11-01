#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: review genes affected by human SNC                              ##
## Data: 2022.11.26                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#gene list collected:
gene_list <- c('ANXA2','ADK','SFXN3','BIRC5','CSTB','ELP6','CSPG5','SST','HES1',
               'RTN4RL2','SNCB','UNC5A','STK33','DENND2B','TMEM9B','ATP1A3','BEGAIN',
               'FGF14','RBFOX2','MAPKAPK3','UCHL3','LMO7','IGSF9B','GRIP1','CAND1',
               'NXPH1','GFRA2','DMTN','LGI3')

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

# generate all candidate gene ---------------------------------------------

#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')
affinity_difference <- readRDS(file = './res/step_88_fig_221124/affinity_difference.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

#get p2g
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
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','InCGE','InMGE','OPC')

#candidate gene
candidate_gene <- c()

#for loop
for (i in cell_type_list) {
  #cell type
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load genes affected by human SNC
  gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC')
  gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
  gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
  gene_list <- readRDS(file = gene_list)
  
  #load DF list
  DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
  DF_list <- readRDS(file = DF_list)
  
  #loop each gene
  for (j in gene_list) {
    gene_name <- j
    
    #peaks related to this gene
    idx_gene <- which(p2g@metadata$geneSet$name == gene_name)
    idx <- which(p2g$idxRNA %in% idx_gene)
    idx_ATAC <- p2g$idxATAC[idx]
    
    peak_list <- p2g@metadata$peakSet[idx_ATAC]
    names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')
    
    #human specific SNC overlap with peak list
    SNC_list <- DF_list[DF_list$group == 'human_specific',"human_SNC"]
    temp <- human_SNC
    names(temp) <- temp$human_SNC
    SNC_list <- temp[SNC_list]
    
    SNC_list <- countOverlaps(query = SNC_list,subject = peak_list)
    SNC_list <- names(SNC_list)[SNC_list > 0]
    
    #check SNC_list has motif binding
    if(sum(SNC_list %in% affinity_difference$seq_name) >= 1){
      candidate_gene <- c(candidate_gene,gene_name)
    }
  }
}

candidate_gene <- unique(candidate_gene)

#save data
saveRDS(object = candidate_gene,file = './res/step_89_fig_221126/candidate_gene.rds')

# check gene --------------------------------------------------------------
#set work space
setwd(dir = '/home/sunym/temp')

#cell type list
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','InCGE','InMGE','OPC')

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = '/home/sunym/temp/merged_Greenleaf_ATAC_ArchR_221018')
macaque_multiome_ArchR <- loadArchRProject(path = '/home/sunym/temp/macaque_multiome_ArchR_221011')

candidate_gene <- readRDS(file = '/content/data/sunym/project/Brain/res/step_89_fig_221126/candidate_gene.rds')

human_SNC <- readRDS(file = '/content/data/sunym/project/Brain/res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = '/content/data/sunym/project/Brain/res/step_78_fig_221111/macaque_SNC.rds')

Greenleaf_RNA_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = '/content/data/sunym/project/Brain/processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

affinity_difference <- readRDS(file = '/content/data/sunym/project/Brain/res/step_88_fig_221124/affinity_difference.rds')
human_fimo_out <- read.table(file = '/content/data/sunym/project/Brain/res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo_out <- read.table(file = '/content/data/sunym/project/Brain/res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)

color_param <- readRDS(file = '/content/data/sunym/project/Brain/data/parameter/shared_param/MetaValue_color_param_220923.rds')
sequence_file <- readRDS(file = '/content/data/sunym/project/Brain/res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')

#add imputation
Greenleaf_ATAC_ArchR <- addImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR,reducedDims = 'IterativeLSI')
macaque_multiome_ArchR <- addImputeWeights(ArchRProj = macaque_multiome_ArchR,reducedDims = 'IterativeLSI')

#modify affinity difference
affinity_difference$seq_name <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))

## ANXA2 -------------------------------------------------------------------
#set param
gene_name <- 'ANXA2'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08377_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## TEDDM1 -------------------------------------------------------------------
#set param
gene_name <- 'TEDDM1'
extention <- 300000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08310_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GCTGTCCCTCCAGCTGACCTA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GCTGTCCCTCCAGCTGACCTA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ADK -------------------------------------------------------------------
#set param
gene_name <- 'ADK'
extention <- 200000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08915_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = human_seq,
                  subject = as('ACTTTCCCTCTCAGAGCTGG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = macaque_seq,
                  subject = as('ACTTTCCCTCTCAGAGCTGG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## SFXN3 -------------------------------------------------------------------
#set param
gene_name <- 'SFXN3'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08342_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = human_seq,
                  subject = as('GCTCGCTCTCCCCGG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = macaque_seq,
                  subject = as('GCTCGCTCTCCCCGG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## PLAAT4 -------------------------------------------------------------------
#set param
gene_name <- 'PLAAT4'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08365_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GCCAGGAGCCTGGAC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GCCAGGAGCCTGGAC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## PLAAT3 -------------------------------------------------------------------
#set param
gene_name <- 'PLAAT3'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08365_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GCCAGGAGCCTGGAC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GCCAGGAGCCTGGAC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## POLA2 -------------------------------------------------------------------
#set param
gene_name <- 'POLA2'
extention <- 200000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07671_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = human_seq,
                  subject = as('CAGCCACCTCAGGCCCCTTGGCAGGCTGCCACCCTAGGC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = macaque_seq,
                  subject = as('CAGCCACCTCAGGCCCCTTGGCAGGCTGCCACCCTAGGC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## SOCS1 -------------------------------------------------------------------
#set param
gene_name <- 'SOCS1'
extention <- 200000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08911_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GTGGAGCCTCAGCCCTAGTCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GTGGAGCCTCAGCCCTAGTCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## LITAF -------------------------------------------------------------------
#set param
gene_name <- 'LITAF'
extention <- 300000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08911_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GTGGAGCCTCAGCCCTAGTCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GTGGAGCCTCAGCCCTAGTCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## BIRC5 -------------------------------------------------------------------
#set param
gene_name <- 'BIRC5'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08377_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## METRNL -------------------------------------------------------------------
#set param
gene_name <- 'METRNL'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08911_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('ACCAGGGCTCCGCTCACAGCGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('ACCAGGGCTCCGCTCACAGCGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## CSTB -------------------------------------------------------------------
#set param
gene_name <- 'CSTB'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07691_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CCTGCGGGAGCTGCGTCCCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CCTGCGGGAGCTGCGTCCCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## GTSE1 -------------------------------------------------------------------
#set param
gene_name <- 'GTSE1'
extention <- 200000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08937_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GAGCCACCCCTCCTAGGCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GAGCCACCCCTCCTAGGCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## CELSR1 -------------------------------------------------------------------
#set param
gene_name <- 'CELSR1'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08937_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GAGCCACCCCTCCTAGGCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GAGCCACCCCTCCTAGGCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ELP6 -------------------------------------------------------------------
#set param
gene_name <- 'ELP6'
extention <- 200000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07587_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GGCCAGCCACACTCTCCTCCTGCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GGCCAGCCACACTCTCCTCCTGCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## CSPG5 -------------------------------------------------------------------
#set param
gene_name <- 'CSPG5'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08377_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## SST -------------------------------------------------------------------
#set param
gene_name <- 'SST'
extention <- 200000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07590_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('ATTCTGGTCTCCCTTTTATGCTCTTTCTTCCTG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('ATTCTGGTCTCCCTTTTATGCTCTTTCTTCCTG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## HES1 -------------------------------------------------------------------
#set param
gene_name <- 'HES1'
extention <- 200000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07610_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('TCTATATGGCCCCTCTTGGCCCCA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('TCTATATGGCCCCTCTTGGCCCCA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## CXCL2 -------------------------------------------------------------------
#set param
gene_name <- 'CXCL2'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08377_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('ATTTGTCTATGTATCTATGTAGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## LINC01503 -------------------------------------------------------------------
#set param
gene_name <- 'LINC01503'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

## PHKA1 -------------------------------------------------------------------
#set param
gene_name <- 'PHKA1'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08955_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CCTCCTGCTGCCCGGGGCTCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CCTCCTGCTGCCCGGGGCTCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ENO1 -------------------------------------------------------------------
#set param
gene_name <- 'ENO1'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07671_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('TGCTGGTCCTGTGCAGCCTGGGCTAATGGTGCCTCCCTA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('TGCTGGTCCTGTGCAGCCTGGGCTAATGGTGCCTCCCTA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## NEXN -------------------------------------------------------------------
#set param
gene_name <- 'NEXN'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08911_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GCCGAGGGACCGGCACCCACTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GCCGAGGGACCGGCACCCACTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## WFDC1 -------------------------------------------------------------------
#set param
gene_name <- 'WFDC1'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M09234_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GAAAGAGAAAGTGGCACCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GAAAGAGAAAGTGGCACCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ATP2C2 -------------------------------------------------------------------
#set param
gene_name <- 'ATP2C2'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M09234_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GAAAGAGAAAGTGGCACCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GAAAGAGAAAGTGGCACCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## KANK1 -------------------------------------------------------------------
#set param
gene_name <- 'KANK1'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07716_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CTCCCTCTTTCCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CTCCCTCTTTTCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ZDHHC21 -------------------------------------------------------------------
#set param
gene_name <- 'ZDHHC21'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07671_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('TTCTGGCTCATGGCCACTTGCACTGCCTCCATGCCCACC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('TTCTGGCTCATGGCCACTTGCACTGCCTCCATGCCCACC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## RTN4RL2 -------------------------------------------------------------------
#set param
gene_name <- 'RTN4RL2'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07773_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('TCCCTCCCTCCCACCCAGCCATACCGTGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('TCCCTCCCTCCCACCCAGCCATACCGTGGT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## IL17C -------------------------------------------------------------------
#set param
gene_name <- 'IL17C'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07691_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('CCTGCAGCCACTGCGGTCCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('CCTGCAGCCACTGCGGTCCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ZNF837 -------------------------------------------------------------------
#set param
gene_name <- 'ZNF837'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07773_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('TGTCTTAACCAGCTCGCTCCCAGGGGCTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('TGTCTTAACCAGCTCGCTCCCAGGGGCTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## CHMP2A -------------------------------------------------------------------
#set param
gene_name <- 'CHMP2A'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07773_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('TGTCTTAACCAGCTCGCTCCCAGGGGCTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('TGTCTTAACCAGCTCGCTCCCAGGGGCTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## MAP2 -------------------------------------------------------------------
#set param
gene_name <- 'MAP2'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08870_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('CTCACCTCGCCCTTTCTCCCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('CTCACCTCGCCCTTTCTCCCCT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## SNCB -------------------------------------------------------------------
#set param
gene_name <- 'SNCB'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07843_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CTCCTCCGAGAGGTCCTCCTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CTCCTCCGAGAGGTCCTCCTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## UNC5A -------------------------------------------------------------------
#set param
gene_name <- 'UNC5A'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07843_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CTCCTCCGAGAGGTCCTCCTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CTCCTCCGAGAGGTCCTCCTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## STK33 -------------------------------------------------------------------
#set param
gene_name <- 'STK33'
extention <- 300000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08310_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GCTTCGCCAGGCCCCCTTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GCTTCGCCAGGCCCCCTTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## DENND2B -------------------------------------------------------------------
#set param
gene_name <- 'DENND2B'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08310_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GCTTCGCCAGGCCCCCTTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GCTTCGCCAGGCCCCCTTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## TMEM9B -------------------------------------------------------------------
#set param
gene_name <- 'TMEM9B'
extention <- 300000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08310_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GCTTCGCCAGGCCCCCTTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GCTTCGCCAGGCCCCCTTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ATP1A3 -------------------------------------------------------------------
#set param
gene_name <- 'ATP1A3'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07605_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CAGTTGCTAACATTTATGTTATTTTATTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CAGTTGCTAACATTTATGTTATTTTATTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## BEGAIN -------------------------------------------------------------------
#set param
gene_name <- 'BEGAIN'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07610_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CCTGATAGCACCCCCCTTCCCCAA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CCTGATAGCACCCCCCTTCCCCAA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## FGF14 -------------------------------------------------------------------
#set param
gene_name <- 'FGF14'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07619_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GAACCCAGCTACGTCCCCCGAACC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GAACCCAGCTACGTCCCCCGAACC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## RBFOX2 -------------------------------------------------------------------
#set param
gene_name <- 'RBFOX2'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08892_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('GCTCTACCTGCCTCTGTCACAC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('GCTCTACCTGCCTCTGTCACAC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## MAPKAPK3 -------------------------------------------------------------------
#set param
gene_name <- 'MAPKAPK3'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07671_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('TGGCTCTGGGAGTCTCCAGCACTGGCGACTGCCCCCTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('TGGCTCTGGGAGTCTCCAGCACTGGCGACTGCCCCCTTT','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## UCHL3 -------------------------------------------------------------------
#set param
gene_name <- 'UCHL3'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08892_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CCTCTTCCCAGCTTTGCACGGC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CCTCTTCCCAGCTTTGCACGGC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## LMO7 -------------------------------------------------------------------
#set param
gene_name <- 'LMO7'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08892_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CCTCTTCCCAGCTTTGCACGGC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CCTCTTCCCAGCTTTGCACGGC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## IGSF9B -------------------------------------------------------------------
#set param
gene_name <- 'IGSF9B'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07625_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('GCTGCCACCGTCCCTCCCTCTGTCTCTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('GCTGCCACCGTCCCTCCCTCTGTCTCTCTC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## GRIP1 -------------------------------------------------------------------
#set param
gene_name <- 'GRIP1'
extention <- 400000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07590_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CATTTTTCCTCACTTTTTAGTGCTTACTCTTTG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CATTTTTCCTCACTTTTTAGTGCTTACTCTTTG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## CAND1 -------------------------------------------------------------------
#set param
gene_name <- 'CAND1'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07590_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CATTTTTCCTCACTTTTTAGTGCTTACTCTTTG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CATTTTTCCTCACTTTTTAGTGCTTACTCTTTG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## ARX -------------------------------------------------------------------
#set param
gene_name <- 'ARX'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07610_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('TTTATTTAACCCCTTGTTCCCCCA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('TTTATTTAACCCCTTGTTCCCCCA','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## NXPH1 -------------------------------------------------------------------
#set param
gene_name <- 'NXPH1'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M07605_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (human_seq),
                  subject = as('AGTTGTTATACATCTTTGTTATTTTCTTCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = (macaque_seq),
                  subject = as('AGTTGTTATACATCTTTGTTATTTTCTTCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## GFRA2 -------------------------------------------------------------------
#set param
gene_name <- 'GFRA2'
extention <- 300000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08362_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('AACGCACGGCCGCCTGGATCCCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('AACGCACGGCCGCCTGGATCCCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## DMTN -------------------------------------------------------------------
#set param
gene_name <- 'DMTN'
extention <- 100000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08362_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('AACGCACGGCCGCCTGGATCCCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('AACGCACGGCCGCCTGGATCCCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

## LGI3 -------------------------------------------------------------------
#set param
gene_name <- 'LGI3'
extention <- 250000
y_lim <- 0.999
loop_resolution = 1000

#gene expression
p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('human',gene_name,sep = ' ')) + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = paste('macaque',gene_name,sep = ' ')) + NoLegend()

p1+p2+plot_layout(ncol = 1)

#gene score
p1 <- plotEmbedding(ArchRProj = Greenleaf_ATAC_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = Greenleaf_ATAC_ArchR),
                    plotAs = 'points',size = 0.1)

p1 <- p1 + NoLegend() + NoAxes() + 
  labs(title = paste('human gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR,
                    colorBy = 'GeneScoreMatrix',
                    name = gene_name,
                    embedding = 'UMAP',
                    imputeWeights = getImputeWeights(ArchRProj = macaque_multiome_ArchR),
                    plotAs = 'points')

p2 <- p2 + NoLegend() + NoAxes() + 
  labs(title = paste('macaque gene score',gene_name,sep = ' ')) + 
  theme(plot.title = element_text(hjust = 0.5))

p1+p2+plot_layout(ncol = 1)

#SNC_list
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

idx_RNA <- which(p2g@metadata$geneSet$name == gene_name)
idx <- which(p2g$idxRNA %in% idx_RNA)
idx_ATAC <- p2g$idxATAC[idx]
peak_list <- p2g@metadata$peakSet[idx_ATAC]
names(peak_list) <- paste(peak_list@seqnames,as.character(peak_list@ranges),sep = '-')

temp <- human_SNC
names(temp) <- temp$human_SNC
SNC_list <- countOverlaps(query = temp,subject = peak_list)
SNC_list <- names(SNC_list)[SNC_list > 0]
temp <- affinity_difference
temp <- temp[temp$seq_name %in% SNC_list,]
SNC_list <- unique(temp$seq_name)
print(paste('SNC_list length:',length(SNC_list),sep = ' '))

#motif_file
motif_name <- 'M08362_2.00'
motif_file <- list.files(path = '/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_file <- motif_file[grep(pattern = motif_name,x = motif_file,fixed = TRUE)]
motif_file <- paste('/content/data/sunym/project/Brain/data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_file,sep = '/')
motif_file <- read_cisbp(file = motif_file)
motif_file@name <- motif_name
motif_file

#sequence alignment
NSM <- matrix(data = 0,nrow = 5,ncol = 5)
rownames(NSM) <- c('A','T','C','G','N')
colnames(NSM) <- c('A','T','C','G','N')
for (i in rownames(NSM)) {
  for (j in colnames(NSM)) {
    if(i == j){
      NSM[i,j] <- 1
    }
  }
}
NSM[,'N'] <- 1
NSM['N',] <- 1

human_seq <- sequence_file$human_seq[[SNC_list[1]]]
macaque_seq <- sequence_file$macaque_seq[[SNC_list[1]]]

human_fimo_out[human_fimo_out$sequence_name == SNC_list[1] & human_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('AACGCACGGCCGCCTGGATCCCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo_out[macaque_fimo_out$sequence_name == SNC_list[1] & macaque_fimo_out$motif_id == motif_name,]
pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('AACGCACGGCCGCCTGGATCCCCC','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#track plot
temp <- human_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])

temp <- macaque_SNC
names(temp) <- temp$human_SNC
p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,groupBy = 'cell_type',useGroups = cell_type_list,
                      plotSummary = c("bulkTrack","featureTrack","loopTrack","geneTrack"),
                      sizes = c(10,0.5,1,1),features = temp[SNC_list],
                      loops = getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = loop_resolution,returnLoops = TRUE),
                      geneSymbol = gene_name,upstream = extention,downstream = extention,
                      ylim = c(0,y_lim),pal = color_param$celltype)
grid::grid.newpage()
grid::grid.draw(p[[1]])