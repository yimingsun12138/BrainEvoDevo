#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: modern human fixed SNC on macaque genome                        ##
## Data: 2022.10.22                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#sleep
ii <- 1
while(1){
  i <- matrix(data = rnorm(n = 2000*2000),nrow = 2000,ncol = 2000)
  i <- i %*% i
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# create macaque GroupCoverage --------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

#create temp
temp <- getArrowFiles(ArchRProj = macaque_multiome_ArchR)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = macaque_multiome_ArchR),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = macaque_multiome_ArchR)
)
temp <- temp[rownames(macaque_multiome_ArchR@cellColData)]
temp$cell_type <- macaque_multiome_ArchR$cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = getPeakSet(ArchRProj = macaque_multiome_ArchR))
temp <- addPeakMatrix(ArchRProj = temp)

#export coverage
temp$ReadsInPromoter <- temp$ReadsInPeaks
coverage_file <- getGroupBW(ArchRProj = temp,groupBy = 'cell_type',normMethod = 'ReadsInPromoter',tileSize = 25,maxCells = NULL,verbose = TRUE)
ii <- gsub(pattern = 'ReadsInPromoter',replacement = 'ReadsInPeaks',fixed = TRUE,x = coverage_file)
for (i in 1:length(coverage_file)) {
  char <- paste('mv',coverage_file[i],ii[i],sep = ' ')
  system(char)
}

# SNC liftover to macaque -------------------------------------------------
#load data
SNC_GRanges_extent <- readRDS(file = './res/step_67_fig_221021/extent_plot/SNC_GRanges.rds')
SNC_GRanges <- readRDS(file = './res/step_67_fig_221021/SNC_GRanges.rds')

#liftover
SNC_GRanges_macaque <- my_rtracklayer_liftOver(ori_GRanges = SNC_GRanges,
                                               chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac10.over.chain',
                                               merge = FALSE)
table(duplicated(SNC_GRanges_macaque$ori_peak))

#add meta data
names(SNC_GRanges_macaque) <- paste(SNC_GRanges_macaque@seqnames,as.character(SNC_GRanges_macaque@ranges),sep = '-')
names(SNC_GRanges) <- paste(SNC_GRanges@seqnames,as.character(SNC_GRanges@ranges),sep = '-')
names(SNC_GRanges_extent) <- names(SNC_GRanges)
SNC_GRanges_macaque@elementMetadata <- cbind(SNC_GRanges_macaque@elementMetadata,
                                             SNC_GRanges_extent[SNC_GRanges_macaque$ori_peak]@elementMetadata)

#get macaque sequence
i <- BSgenome::getSeq(x = BSgenome.Mmulatta.UCSC.rheMac10,
                      name = SNC_GRanges_macaque,
                      as.character = TRUE)
SNC_GRanges_macaque$macaque_nucleotide <- i

#find sequence difference
table(SNC_GRanges_macaque$modern_nucleotide == SNC_GRanges_macaque$macaque_nucleotide)
table(SNC_GRanges_macaque$archaic_nucleotide == SNC_GRanges_macaque$macaque_nucleotide)


## display -----------------------------------------------------------------
temp_modern <- paste(names(SNC_GRanges_macaque),SNC_GRanges_macaque$modern_nucleotide)
temp_archaic <- paste(names(SNC_GRanges_macaque),SNC_GRanges_macaque$archaic_nucleotide)
temp_macaque <- paste(names(SNC_GRanges_macaque),SNC_GRanges_macaque$macaque_nucleotide)
temp <- list(modern_human = temp_modern,
             archaic_human = temp_archaic,
             macaque = temp_macaque)

#venn plot
pdf(file = './res/step_69_fig_221022/human_macaque_SNC_sequence_overlap_vennplot.pdf',width = 6,height = 6)
ggvenn(data = temp,columns = c('modern_human','archaic_human','macaque'))
dev.off()

#upset plot
pdf(file = './res/step_69_fig_221022/human_macaque_SNC_sequence_overlap_upset_plot.pdf',width = 6,height = 3)
upset(data = fromList(input = temp),
      sets = c('modern_human','archaic_human','macaque'),
      order.by = "freq",keep.order = TRUE)
dev.off()

#subset SNC GRanges
SNC_GRanges_macaque <- SNC_GRanges_macaque[!(SNC_GRanges_macaque$modern_nucleotide == SNC_GRanges_macaque$macaque_nucleotide)]


## modify SNC macaque ------------------------------------------------------
#load data
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#reduce SNC
cell_type_list <- names(color_param$celltype)
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)
SNC_GRanges_macaque@elementMetadata <- SNC_GRanges_macaque@elementMetadata[,!(colnames(SNC_GRanges_macaque@elementMetadata) %in% cell_type_list)]

#extent
temp <- rtracklayer::as.data.frame(x = SNC_GRanges_macaque)
temp <- temp[,c("seqnames","start","end")]
temp$start <- temp$start-250
temp$end <- temp$end+250
colnames(temp) <- c('chrom','start','end')
temp <- as(temp,'GRanges')
names(temp) <- names(SNC_GRanges_macaque)
temp@elementMetadata <- SNC_GRanges_macaque@elementMetadata
SNC_GRanges_macaque <- temp

#which cell type enrich these SNC
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

for (i in cell_type_list) {
  file_list <- paste0('./processed_data/221008_summary/macaque_multiome_ArchR_221011/PeakCalls/',i,'-reproduciblePeaks.gr.rds')
  print(file.exists(file_list))
  temp <- paste(i,'peak',sep = '.')
  assign(x = temp,value = readRDS(file_list))
  temp <- countOverlaps(query = SNC_GRanges_macaque,subject = get(temp))
  SNC_GRanges_macaque@elementMetadata[,i] <- 'NO'
  SNC_GRanges_macaque@elementMetadata[temp > 0,i] <- 'YES'
}

temp <- data.frame(cell_type = cell_type_list,
                   SNC = unlist(base::lapply(X = cell_type_list,FUN = function(x){
                     sum(SNC_GRanges_macaque@elementMetadata[,x] == 'YES')
                   })))
temp$cell_type <- gsub(pattern = '.',replacement = '-',x = temp$cell_type,fixed = TRUE)
temp$cell_type <- factor(temp$cell_type,levels = rev(temp$cell_type))

pdf(file = './res/step_69_fig_221022/human_fixed_SNC_macaque_cell_type_enrichment.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x = cell_type,y = SNC,fill = cell_type)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8) + 
  scale_fill_manual(values = color_param$celltype[as.character(temp$cell_type)]) + 
  theme_bw() + 
  ylim(c(0,800)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),
        aspect.ratio = 2.5) + 
  xlab('Cell type') + ylab('SNC number') + 
  coord_flip()
dev.off()

#save data
SNC_GRanges_macaque$macaque_ori_peak <- names(SNC_GRanges_macaque)
SNC_GRanges_extent$ori_peak <- names(SNC_GRanges_extent)
SNC_GRanges_extent <- SNC_GRanges_extent[SNC_GRanges_macaque$ori_peak]

saveRDS(object = SNC_GRanges_macaque,file = './res/step_69_fig_221022/SNC_macaque.rds')
saveRDS(object = SNC_GRanges_extent,file = './res/step_69_fig_221022/SNC_human.rds')

# human SNC distribution on macaque genome --------------------------------
#load data
SNC_macaque <- readRDS(file = './res/step_69_fig_221022/SNC_macaque.rds')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

## macaque SNC cell type enrichment intersect ----------------------------------------
#UPset plot
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- base::lapply(X = cell_type_list,FUN = function(x){
  temp <- names(SNC_macaque)[SNC_macaque@elementMetadata[,x] == 'YES']
  return(temp)
})

cell_type_list <- gsub(pattern = '.',replacement = '-',x = cell_type_list,fixed = TRUE)
names(temp) <- cell_type_list

#plot
pdf(file = './res/step_69_fig_221022/human_fixed_SNC_macaque_cell_type_enrichment_overlap.pdf',width = 12,height = 7)
upset(data = fromList(temp),sets = cell_type_list,order.by = "freq",keep.order = TRUE)
dev.off()

#overlap number
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- rtracklayer::as.data.frame(x = SNC_macaque)
temp <- temp[,cell_type_list]
temp <- base::lapply(X = rownames(temp),FUN = function(x){
  return(sum(temp[x,] == 'YES'))
})
temp <- unlist(temp)

temp <- as.data.frame(table(temp))
colnames(temp) <- c('Intersect','percent')
temp$percent <- temp$percent/sum(temp$percent)*100
temp$Intersect <- factor(temp$Intersect,levels = temp$Intersect)

pdf(file = './res/step_69_fig_221022/human_fixed_SNC_macaque_cell_type_active_number_barplot.pdf',width = 5,height = 3)
ggplot(data = temp,aes(x = Intersect,y = percent)) + 
  geom_bar(position = "stack", stat = "identity",width = 0.8,color = 'grey') + 
  theme_bw() + 
  scale_y_break(breaks = c(10,70)) + 
  theme(panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        panel.border = element_rect(size = 0.8),
        axis.ticks.y = element_blank(),
        text = element_text(size = 10),
        line = element_line(size = 0.8),
        legend.position = 'none') + 
  xlab('cell type size') + ylab('Percent %')
dev.off()

## human fixed SNC macaque genome feature ----------------------------------

#make TxDb
gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
gene_anno <- rtracklayer::as.data.frame(gene_anno)
gene_anno$seqnames <- as.character(gene_anno$seqnames)
gene_anno$seqnames <- unlist(base::lapply(X = gene_anno$seqnames,FUN = function(x){
  if(x %in% c(as.character(1:20),'X','Y','MT')){
    return(paste0('chr',x))
  }else{
    return(x)
  }
}))
gene_anno <- makeGRangesFromDataFrame(df = gene_anno,keep.extra.columns = TRUE)
gene_anno <- GenomicFeatures::makeTxDbFromGRanges(gene_anno)

#global distribution
peakAnno <- annotatePeak(SNC_macaque,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Mmu.eg.db")

pdf(file = './res/step_69_fig_221022/human_fixed_SNC_macaque_genome_distribution.pdf',width = 6,height = 4)
plotAnnoPie(peakAnno)
dev.off()

#active SNC in each cell type
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

SNC_list <- base::lapply(X = cell_type_list,FUN = function(x){
  return(SNC_macaque[SNC_macaque@elementMetadata[,x] == 'YES'])
})
cell_type_list <- gsub(pattern = '.',replacement = '-',x = cell_type_list,fixed = TRUE)
names(SNC_list) <- cell_type_list

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Mmu.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_69_fig_221022/macaque_cell_type_active_SNC_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

#SNC with different active intersect size
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type]
cell_type_list <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)

temp <- as.data.frame(SNC_macaque@elementMetadata[,cell_type_list])
temp <- rowSums(temp == 'YES')
SNC_macaque$intersect <- temp

SNC_list <- base::lapply(X = 0:15,FUN = function(x){
  return(SNC_macaque[SNC_macaque$intersect == x])
})
names(SNC_list) <- 0:15

peakAnno <- base::lapply(X = SNC_list,FUN = annotatePeak,tssRegion=c(-2000,2000),TxDb=gene_anno,annoDb="org.Mmu.eg.db")
p <- plotAnnoBar(peakAnno)

pdf(file = './res/step_69_fig_221022/macaque_SNC_with_intersect_size_genome_distribution.pdf',width = 6,height = 4.5)
p + theme(aspect.ratio = 1.5)
dev.off()

# accessibility difference between human and macaque ----------------------

## create peak matrix ------------------------------------------------------
#load data
SNC_macaque <- readRDS(file = './res/step_69_fig_221022/SNC_macaque.rds')
SNC_human <- readRDS(file = './res/step_69_fig_221022/SNC_human.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

#get macaque peak file
temp <- getArrowFiles(ArchRProj = macaque_multiome_ArchR)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/macaque_multiome_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = macaque_multiome_ArchR),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = macaque_multiome_ArchR)
)
temp <- temp[rownames(macaque_multiome_ArchR@cellColData)]
temp$cell_type <- macaque_multiome_ArchR$cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = SNC_macaque)
temp <- addPeakMatrix(ArchRProj = temp)
getAvailableMatrices(ArchRProj = temp)

macaque_peak_matrix <- getMatrixFromProject(ArchRProj = temp,useMatrix = 'PeakMatrix',verbose = TRUE)
gene_list <- paste(macaque_peak_matrix@rowRanges@seqnames,as.character(macaque_peak_matrix@rowRanges@ranges),sep = '-')
macaque_peak_matrix <- macaque_peak_matrix@assays@data$PeakMatrix
rownames(macaque_peak_matrix) <- gene_list

saveRDS(object = macaque_peak_matrix,file = './res/step_69_fig_221022/macaque_peak_matrix.rds')

#get human peak file
temp <- getArrowFiles(ArchRProj = Greenleaf_ATAC_ArchR)
file.exists(temp)
temp <- ArchRProject(
  ArrowFiles = temp,
  outputDirectory = '/home/sunym/temp/Greenleaf_ATAC_ArchR',
  copyArrows = TRUE,
  geneAnnotation = getGeneAnnotation(ArchRProj = Greenleaf_ATAC_ArchR),
  genomeAnnotation = getGenomeAnnotation(ArchRProj = Greenleaf_ATAC_ArchR)
)
temp <- temp[rownames(Greenleaf_ATAC_ArchR@cellColData)]
temp$cell_type <- Greenleaf_ATAC_ArchR$cell_type
table(temp$cell_type)
temp <- addPeakSet(ArchRProj = temp,peakSet = SNC_human)
temp <- addPeakMatrix(ArchRProj = temp)
getAvailableMatrices(ArchRProj = temp)

human_peak_matrix <- getMatrixFromProject(ArchRProj = temp,useMatrix = 'PeakMatrix',verbose = TRUE)
gene_list <- paste(human_peak_matrix@rowRanges@seqnames,as.character(human_peak_matrix@rowRanges@ranges),sep = '-')
human_peak_matrix <- human_peak_matrix@assays@data$PeakMatrix
rownames(human_peak_matrix) <- gene_list

saveRDS(object = human_peak_matrix,file = './res/step_69_fig_221022/human_peak_matrix.rds')

## find differentially active between human and macaque --------------------
#load data
SNC_human <- readRDS(file = './res/step_69_fig_221022/SNC_human.rds')
SNC_macaque <- readRDS(file = './res/step_69_fig_221022/SNC_macaque.rds')

human_peak_matrix <- readRDS(file = './res/step_69_fig_221022/human_peak_matrix.rds')
macaque_peak_matrix <- readRDS(file = './res/step_69_fig_221022/macaque_peak_matrix.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#rename peak matrix
names(SNC_human) <- paste(SNC_human@seqnames,as.character(SNC_human@ranges),sep = '-')
names(SNC_macaque) <- paste(SNC_macaque@seqnames,as.character(SNC_macaque@ranges),sep = '-')

table(rownames(macaque_peak_matrix) %in% names(SNC_macaque))
table(rownames(human_peak_matrix) %in% names(SNC_human))

rownames(macaque_peak_matrix) <- SNC_macaque[rownames(macaque_peak_matrix)]$ori_peak
rownames(human_peak_matrix) <- SNC_human[rownames(human_peak_matrix)]$ori_peak
gene_list <- rownames(human_peak_matrix)
macaque_peak_matrix <- macaque_peak_matrix[gene_list,]


### IP ----------------------------------------------------------------------
cell_type <- 'IP'

#subset data
temp_human <- human_peak_matrix[,as.character(Greenleaf_ATAC_ArchR@cellColData[colnames(human_peak_matrix),'cell_type']) == cell_type]
temp_macaque <- macaque_peak_matrix[,as.character(macaque_multiome_ArchR@cellColData[colnames(macaque_peak_matrix),"cell_type"]) == cell_type]

print(paste('human cell num:',ncol(temp_human),sep = ' '))
print(paste('macaque cell num:',ncol(temp_macaque),sep = ' '))

#DF
temp_DF <- my_DF_wilcox_test(mat1 = temp_human,mat2 = temp_macaque,alternative = 'two.sided',
                             paired = FALSE,workers = 4,future.globals.maxSize = 200*(1024^3))

#class
human_active <- rownames(temp_DF)[temp_DF$log2FC > 1 & temp_DF$fdr < 0.01]
temp <- SNC_human
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_active <- human_active[temp[human_active]@elementMetadata[,cell_type_dot] == 'YES']

macaque_active <- rownames(temp_DF)[temp_DF$log2FC < -1 & temp_DF$fdr < 0.01]
temp <- SNC_macaque
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
macaque_active <- macaque_active[temp[macaque_active]@elementMetadata[,cell_type_dot] == 'YES']

all_active <- rownames(temp_DF)[abs(temp_DF$log2FC) <= 1]
temp <- SNC_human
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
all_active <- all_active[temp[all_active]@elementMetadata[,cell_type_dot] == 'YES']
temp <- SNC_macaque
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
all_active <- all_active[temp[all_active]@elementMetadata[,cell_type_dot] == 'YES']

#save data
temp <- list(human = human_active,macaque = macaque_active,all = all_active)
char <- paste0('./res/step_69_fig_221022/DF_SNC/',cell_type_dot,'.SNC_list.rds')
saveRDS(object = temp,file = char)

#enrich heatmap
#human
SNC_list <- c(human_active,macaque_active,all_active)
SNC_list <- data.frame(chrom = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][1])
})),
start = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})),
end = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})))
SNC_list <- as(SNC_list,'GRanges')
SNC_list$group <- c(rep('human',times = length(human_active)),
                    rep('macaque',times = length(macaque_active)),
                    rep('all',times = length(all_active)))

signal_coverage <- list.files(path = './res/step_67_fig_221021/GroupBigWigs/cell_type/')
signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
signal_coverage <- paste('./res/step_67_fig_221021/GroupBigWigs/cell_type',signal_coverage,sep = '/')
file.exists(signal_coverage)

signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_list,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                      use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(4,'inches'))

#macaque
SNC_list <- c(human_active,macaque_active,all_active)
temp <- SNC_macaque
names(temp) <- temp$ori_peak
SNC_list <- temp[SNC_list]$macaque_ori_peak
SNC_list <- data.frame(chrom = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][1])
})),
start = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})),
end = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})))
SNC_list <- as(SNC_list,'GRanges')
SNC_list$group <- c(rep('human',times = length(human_active)),
                    rep('macaque',times = length(macaque_active)),
                    rep('all',times = length(all_active)))

signal_coverage <- list.files(path = './res/step_69_fig_221022/GroupBigWigs/cell_type/')
signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
signal_coverage <- paste('./res/step_69_fig_221022/GroupBigWigs/cell_type',signal_coverage,sep = '/')
file.exists(signal_coverage)

signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_list,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                      use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(4,'inches'))

#plot
col_fun <- colorRamp2(breaks = c(0,0.01,0.02),colors = c('#FFC30F','#C70039','#581845'))
char <- paste0('./res/step_69_fig_221022/cell_type_plot/',cell_type_dot,'.enrichheatmap_between_species.pdf')
row_order_list <- list(human = row_order(p1)$human,
                       macaque = row_order(p2)$macaque,
                       all = row_order(p1)$all)
row_order_list <- unlist(row_order_list)

pdf(file = char,width = 4,height = 5.5)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE,
                width = unit(1,'inches'),height = unit(4,'inches'),row_order = row_order_list) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  width = unit(1,'inches'),height = unit(4,'inches'),row_order = row_order_list)
dev.off()

#dot line plot
my_insertion_plot <- function(group,p1 = p1,p2 = p2){
  temp <- p1@matrix[SNC_list$group == group,]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[SNC_list$group == group,]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = color_param$species[c('human','macaque')]) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human <- my_insertion_plot(group = 'human',p1 = p1,p2 = p2)
p_macaque <- my_insertion_plot(group = 'macaque',p1 = p1,p2 = p2)
p_all <- my_insertion_plot(group = 'all',p1 = p1,p2 = p2)

char <- paste0('./res/step_69_fig_221022/cell_type_plot/',cell_type_dot,'.dot_line_plot_between_species.pdf')
pdf(file = char,width = 6,height = 8)
p_human + p_macaque + p_all + plot_layout(ncol = 1)
dev.off()

### RG-1 ----------------------------------------------------------------------
cell_type <- 'RG-1'

#subset data
temp_human <- human_peak_matrix[,as.character(Greenleaf_ATAC_ArchR@cellColData[colnames(human_peak_matrix),'cell_type']) == cell_type]
temp_macaque <- macaque_peak_matrix[,as.character(macaque_multiome_ArchR@cellColData[colnames(macaque_peak_matrix),"cell_type"]) == cell_type]

print(paste('human cell num:',ncol(temp_human),sep = ' '))
print(paste('macaque cell num:',ncol(temp_macaque),sep = ' '))

#DF
temp_DF <- my_DF_wilcox_test(mat1 = temp_human,mat2 = temp_macaque,alternative = 'two.sided',
                             paired = FALSE,workers = 4,future.globals.maxSize = 200*(1024^3))

#class
human_active <- rownames(temp_DF)[temp_DF$log2FC > 1 & temp_DF$fdr < 0.01]
temp <- SNC_human
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
human_active <- human_active[temp[human_active]@elementMetadata[,cell_type_dot] == 'YES']

macaque_active <- rownames(temp_DF)[temp_DF$log2FC < -1 & temp_DF$fdr < 0.01]
temp <- SNC_macaque
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
macaque_active <- macaque_active[temp[macaque_active]@elementMetadata[,cell_type_dot] == 'YES']

all_active <- rownames(temp_DF)[abs(temp_DF$log2FC) <= 1]
temp <- SNC_human
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
all_active <- all_active[temp[all_active]@elementMetadata[,cell_type_dot] == 'YES']
temp <- SNC_macaque
names(temp) <- temp$ori_peak
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
all_active <- all_active[temp[all_active]@elementMetadata[,cell_type_dot] == 'YES']

#save data
temp <- list(human = human_active,macaque = macaque_active,all = all_active)
char <- paste0('./res/step_69_fig_221022/DF_SNC/',cell_type_dot,'.SNC_list.rds')
saveRDS(object = temp,file = char)

#enrich heatmap
#human
SNC_list <- c(human_active,macaque_active,all_active)
SNC_list <- data.frame(chrom = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][1])
})),
start = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})),
end = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})))
SNC_list <- as(SNC_list,'GRanges')
SNC_list$group <- c(rep('human',times = length(human_active)),
                    rep('macaque',times = length(macaque_active)),
                    rep('all',times = length(all_active)))

signal_coverage <- list.files(path = './res/step_67_fig_221021/GroupBigWigs/cell_type/')
signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
signal_coverage <- paste('./res/step_67_fig_221021/GroupBigWigs/cell_type',signal_coverage,sep = '/')
file.exists(signal_coverage)

signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_list,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p1 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                      use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(5,'inches'))

#macaque
SNC_list <- c(human_active,macaque_active,all_active)
temp <- SNC_macaque
names(temp) <- temp$ori_peak
SNC_list <- temp[SNC_list]$macaque_ori_peak
SNC_list <- data.frame(chrom = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][1])
})),
start = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})),
end = unlist(base::lapply(X = SNC_list,FUN = function(x){
  temp <- strsplit(x = x,split = '-')
  return(temp[[1]][2])
})))
SNC_list <- as(SNC_list,'GRanges')
SNC_list$group <- c(rep('human',times = length(human_active)),
                    rep('macaque',times = length(macaque_active)),
                    rep('all',times = length(all_active)))

signal_coverage <- list.files(path = './res/step_69_fig_221022/GroupBigWigs/cell_type/')
signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
signal_coverage <- paste('./res/step_69_fig_221022/GroupBigWigs/cell_type',signal_coverage,sep = '/')
file.exists(signal_coverage)

signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_list,extend = 2000,w = 50,limit = NA,
                         value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)

p2 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                      use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(5,'inches'))

#plot
col_fun <- colorRamp2(breaks = c(0,0.01,0.02),colors = c('#FFC30F','#C70039','#581845'))
char <- paste0('./res/step_69_fig_221022/cell_type_plot/',cell_type_dot,'.enrichheatmap_between_species.pdf')
row_order_list <- list(human = row_order(p1)$human,
                       macaque = row_order(p2)$macaque,
                       all = row_order(p1)$all)
row_order_list <- unlist(row_order_list)

pdf(file = char,width = 4,height = 5.5)
EnrichedHeatmap(mat = p1@matrix,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE,
                width = unit(1,'inches'),height = unit(5,'inches'),row_order = row_order_list) + 
  EnrichedHeatmap(mat = p2@matrix,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                  use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                  width = unit(1,'inches'),height = unit(5,'inches'),row_order = row_order_list)
dev.off()

#dot line plot
my_insertion_plot <- function(group,p1 = p1,p2 = p2){
  temp <- p1@matrix[SNC_list$group == group,]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'human'
  insertion_matrix <- temp
  
  temp <- p2@matrix[SNC_list$group == group,]
  temp <- colMeans(temp)
  temp <- data.frame(temp)
  temp$pos <- c(-40:-1,1:40)
  colnames(temp) <- c('insertion','pos')
  temp$species <- 'macaque'
  insertion_matrix <- rbind(insertion_matrix,temp)
  
  p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
    geom_point(alpha = 0) + 
    geom_line() + 
    theme_cowplot() + 
    scale_color_manual(values = color_param$species[c('human','macaque')]) + 
    theme(aspect.ratio = 1,
          panel.background = element_rect(fill = NA,colour = 'black'),
          axis.line = element_blank(),
          plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
    labs(title = group) + xlab('normalized insertion') + ylab('position')
  return(p)
}

p_human <- my_insertion_plot(group = 'human',p1 = p1,p2 = p2)
p_macaque <- my_insertion_plot(group = 'macaque',p1 = p1,p2 = p2)
p_all <- my_insertion_plot(group = 'all',p1 = p1,p2 = p2)

char <- paste0('./res/step_69_fig_221022/cell_type_plot/',cell_type_dot,'.dot_line_plot_between_species.pdf')
pdf(file = char,width = 6,height = 8)
p_human + p_macaque + p_all + plot_layout(ncol = 1)
dev.off()


### for loop ----------------------------------------------------------------
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% macaque_multiome_ArchR$cell_type & cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type]

for (i in cell_type_list) {
  
  #cell type
  cell_type <- i
  
  #subset data
  temp_human <- human_peak_matrix[,as.character(Greenleaf_ATAC_ArchR@cellColData[colnames(human_peak_matrix),'cell_type']) == cell_type]
  temp_macaque <- macaque_peak_matrix[,as.character(macaque_multiome_ArchR@cellColData[colnames(macaque_peak_matrix),"cell_type"]) == cell_type]
  
  print(paste('human cell num:',ncol(temp_human),sep = ' '))
  print(paste('macaque cell num:',ncol(temp_macaque),sep = ' '))
  
  #DF
  temp_DF <- my_DF_wilcox_test(mat1 = temp_human,mat2 = temp_macaque,alternative = 'two.sided',
                               paired = FALSE,workers = 4,future.globals.maxSize = 200*(1024^3))
  
  #class
  human_active <- rownames(temp_DF)[temp_DF$log2FC > 1 & temp_DF$fdr < 0.01]
  temp <- SNC_human
  names(temp) <- temp$ori_peak
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  human_active <- human_active[temp[human_active]@elementMetadata[,cell_type_dot] == 'YES']
  
  macaque_active <- rownames(temp_DF)[temp_DF$log2FC < -1 & temp_DF$fdr < 0.01]
  temp <- SNC_macaque
  names(temp) <- temp$ori_peak
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  macaque_active <- macaque_active[temp[macaque_active]@elementMetadata[,cell_type_dot] == 'YES']
  
  all_active <- rownames(temp_DF)[abs(temp_DF$log2FC) <= 1]
  temp <- SNC_human
  names(temp) <- temp$ori_peak
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  all_active <- all_active[temp[all_active]@elementMetadata[,cell_type_dot] == 'YES']
  temp <- SNC_macaque
  names(temp) <- temp$ori_peak
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  all_active <- all_active[temp[all_active]@elementMetadata[,cell_type_dot] == 'YES']
  
  #save data
  temp <- list(human = human_active,macaque = macaque_active,all = all_active)
  char <- paste0('./res/step_69_fig_221022/DF_SNC/',cell_type_dot,'.SNC_list.rds')
  saveRDS(object = temp,file = char)
  
  #enrich heatmap
  #human
  SNC_list <- c(human_active,macaque_active,all_active)
  SNC_list <- data.frame(chrom = unlist(base::lapply(X = SNC_list,FUN = function(x){
    temp <- strsplit(x = x,split = '-')
    return(temp[[1]][1])
  })),
  start = unlist(base::lapply(X = SNC_list,FUN = function(x){
    temp <- strsplit(x = x,split = '-')
    return(temp[[1]][2])
  })),
  end = unlist(base::lapply(X = SNC_list,FUN = function(x){
    temp <- strsplit(x = x,split = '-')
    return(temp[[1]][2])
  })))
  SNC_list <- as(SNC_list,'GRanges')
  SNC_list$group <- c(rep('human',times = length(human_active)),
                      rep('macaque',times = length(macaque_active)),
                      rep('all',times = length(all_active)))
  print(table(SNC_list$group))
  
  signal_coverage <- list.files(path = './res/step_67_fig_221021/GroupBigWigs/cell_type/')
  signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
  signal_coverage <- paste('./res/step_67_fig_221021/GroupBigWigs/cell_type',signal_coverage,sep = '/')
  file.exists(signal_coverage)
  
  signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
  mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  
  p1 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                        use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(5.2,'inches'))
  
  #macaque
  SNC_list <- c(human_active,macaque_active,all_active)
  temp <- SNC_macaque
  names(temp) <- temp$ori_peak
  SNC_list <- temp[SNC_list]$macaque_ori_peak
  SNC_list <- data.frame(chrom = unlist(base::lapply(X = SNC_list,FUN = function(x){
    temp <- strsplit(x = x,split = '-')
    return(temp[[1]][1])
  })),
  start = unlist(base::lapply(X = SNC_list,FUN = function(x){
    temp <- strsplit(x = x,split = '-')
    return(temp[[1]][2])
  })),
  end = unlist(base::lapply(X = SNC_list,FUN = function(x){
    temp <- strsplit(x = x,split = '-')
    return(temp[[1]][2])
  })))
  SNC_list <- as(SNC_list,'GRanges')
  SNC_list$group <- c(rep('human',times = length(human_active)),
                      rep('macaque',times = length(macaque_active)),
                      rep('all',times = length(all_active)))
  
  signal_coverage <- list.files(path = './res/step_69_fig_221022/GroupBigWigs/cell_type/')
  signal_coverage <- signal_coverage[grep(pattern = cell_type_dot,x = signal_coverage,fixed = TRUE)]
  signal_coverage <- paste('./res/step_69_fig_221022/GroupBigWigs/cell_type',signal_coverage,sep = '/')
  file.exists(signal_coverage)
  
  signal_coverage <- rtracklayer::import.bw(con = signal_coverage)
  mat <- normalizeToMatrix(signal = signal_coverage,target = SNC_list,extend = 2000,w = 50,limit = NA,
                           value_column = 'score',background = 0,mean_mode = 'w0',smooth = TRUE,verbose = TRUE)
  
  p2 <- EnrichedHeatmap(mat = mat,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                        use_raster = TRUE,raster_resize_mat = mean,width = unit(1,'inches'),height = unit(5.2,'inches'))
  
  #plot
  col_fun <- colorRamp2(breaks = c(0,0.01,0.02),colors = c('#FFC30F','#C70039','#581845'))
  char <- paste0('./res/step_69_fig_221022/cell_type_plot/',cell_type_dot,'.enrichheatmap_between_species.pdf')
  row_order_list <- list(human = row_order(p1)$human,
                         macaque = row_order(p2)$macaque,
                         all = row_order(p1)$all)
  row_order_list <- unlist(row_order_list)
  
  pdf(file = char,width = 4,height = 7.5)
  print(EnrichedHeatmap(mat = p1@matrix,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                        use_raster = TRUE,raster_resize_mat = mean,col = col_fun,name = 'insertion',pos_line = FALSE,
                        width = unit(1,'inches'),height = unit(5.2,'inches'),row_order = row_order_list) + 
          EnrichedHeatmap(mat = p2@matrix,row_split = factor(SNC_list$group,levels = c('human','macaque','all')),
                          use_raster = TRUE,raster_resize_mat = mean,col = col_fun,show_heatmap_legend = FALSE,pos_line = FALSE,
                          width = unit(1,'inches'),height = unit(5.2,'inches'),row_order = row_order_list))
  dev.off()
  
  #dot line plot
  my_insertion_plot <- function(group,p1 = p1,p2 = p2){
    temp <- p1@matrix[SNC_list$group == group,]
    temp <- colMeans(temp)
    temp <- data.frame(temp)
    temp$pos <- c(-40:-1,1:40)
    colnames(temp) <- c('insertion','pos')
    temp$species <- 'human'
    insertion_matrix <- temp
    
    temp <- p2@matrix[SNC_list$group == group,]
    temp <- colMeans(temp)
    temp <- data.frame(temp)
    temp$pos <- c(-40:-1,1:40)
    colnames(temp) <- c('insertion','pos')
    temp$species <- 'macaque'
    insertion_matrix <- rbind(insertion_matrix,temp)
    
    p <- ggplot(data = insertion_matrix,aes(x=pos,y=insertion,color=species)) + 
      geom_point(alpha = 0) + 
      geom_line() + 
      theme_cowplot() + 
      scale_color_manual(values = color_param$species[c('human','macaque')]) + 
      theme(aspect.ratio = 1,
            panel.background = element_rect(fill = NA,colour = 'black'),
            axis.line = element_blank(),
            plot.title = element_text(face = 'bold',size = 14,hjust = 0.5)) + 
      labs(title = group) + xlab('normalized insertion') + ylab('position')
    return(p)
  }
  
  p_human <- my_insertion_plot(group = 'human',p1 = p1,p2 = p2)
  p_macaque <- my_insertion_plot(group = 'macaque',p1 = p1,p2 = p2)
  p_all <- my_insertion_plot(group = 'all',p1 = p1,p2 = p2)
  
  char <- paste0('./res/step_69_fig_221022/cell_type_plot/',cell_type_dot,'.dot_line_plot_between_species.pdf')
  pdf(file = char,width = 6,height = 8)
  print(p_human + p_macaque + p_all + plot_layout(ncol = 1))
  dev.off()
  
  #end
  print(paste(cell_type,'done!',sep = ' '))
}
