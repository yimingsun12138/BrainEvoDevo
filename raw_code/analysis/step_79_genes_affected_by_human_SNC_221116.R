#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: genes affected by human SNC                                     ##
## Data: 2022.11.16                                                                ##
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

# gene regulated by human specific SNC intersect between cell type --------

#load data
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)

for (i in cell_type_list) {
  #cell type
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type)
  
  #get gene list
  gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
  gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
  gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
  gene_list <- readRDS(file = gene_list)
  
  #assign
  char <- paste(cell_type_dot,'_gene_list',sep = '')
  assign(x = char,value = gene_list)
}

#upset plot
gene_list <- list(RG.1_gene_list,RG.2_gene_list,IP_gene_list,Ex.1_gene_list,Ex.2_gene_list,Ex.3_gene_list,Ex.4_gene_list,OPC_gene_list)
names(gene_list) <- cell_type_list_dot

pdf(file = './res/step_79_fig_221116/gene_regulated_by_human_specific_SNC_cell_type_upset_plot.pdf',width = 10,height = 5)
upset(fromList(gene_list),order.by = "freq",sets = cell_type_list_dot,keep.order = TRUE)
dev.off()

# try transplotR ----------------------------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

human_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
human_anno <- rtracklayer::as.data.frame(x = human_anno)

#add gene quantile
temp <- readRDS(file = './res/step_79_fig_221112/Greenleaf_RNA_quantile_matrix.rds')
temp <- CreateAssayObject(data = temp,min.cells = 0,min.features = 0)
Greenleaf_RNA_Seurat[['quantile']] <- temp
rm(temp)
gc()

temp <- readRDS(file = './res/step_79_fig_221112/macque_RNA_quantile_matrix.rds')
temp <- CreateAssayObject(data = temp,min.cells = 0,min.features = 0)
macaque_multiome_Seurat[['quantile']] <- temp
rm(temp)
gc()

#RG-1 up regulated gene
cell_type <- 'RG-1'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
DF_list <- readRDS(file = DF_list)

#get gene list
gene_list <- list.files(path = './res/step_79_fig_221112/up_gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/up_gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(file = gene_list)

temp_human <- Greenleaf_RNA_Seurat@assays$quantile@data[gene_list,Greenleaf_RNA_Seurat$ReAnno_celltype == cell_type]
temp_macaque <- macaque_multiome_Seurat@assays$quantile@data[gene_list,macaque_multiome_Seurat$cell_type == cell_type]
temp <- cbind(temp_human,temp_macaque)
temp <- CreateSeuratObject(counts = temp,project = 'quantile',assay = 'RNA',min.cells = 0,min.features = 0)
temp$species <- NA
temp@meta.data[colnames(temp_human),"species"] <- 'human'
temp@meta.data[colnames(temp_macaque),"species"] <- 'macaque'

pdf(file = '/home/sunym/temp/temp.pdf',width = 16,height = (round(length(gene_list)/4)+1)*4)
VlnPlot(object = temp,features = gene_list,cols = color_param$species[c('human','macaque')],pt.size = 0,assay = 'RNA',group.by = 'species',ncol = 4,slot = 'counts')
dev.off()

rm(Greenleaf_RNA_Seurat)
rm(macaque_multiome_Seurat)
gc()

#try gene CSTB
gene_name <- 'CSTB'
extend_size <- 100000

#load bw
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)

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

#get target site
TSS_site <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]

chrom <- as.character(TSS_site@seqnames)
start_site <- TSS_site@ranges@start - extend_size
end_site <- TSS_site@ranges@start + extend_size

#plot gene annotation
p_gene_anno <- trancriptVis(gtfFile = human_anno,collapse = TRUE,
                            Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                            posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                            addNormalArrow = FALSE,text.pos = 'middle',
                            newStyleArrow = FALSE,xAxis.info = TRUE)

#plot peak
temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
temp <- human_SNC
names(temp) <- temp$human_SNC
temp <- temp[DF_list[DF_list$group == 'human_specific',"human_SNC"]]
rtracklayer::export(object = temp,con = '/home/sunym/temp/SNC.bed',format = 'bed')

p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed','/home/sunym/temp/SNC.bed'),
                      chr = chrom,region.min = start_site,region.max = end_site,
                      show.legend = FALSE,fill = ggsci::pal_d3()(2),track.width = 0.3)

#plot track
p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = human_anno,
                         chr = chrom,region.min = start_site,region.max = end_site,
                         sample.order = cell_type_list,space.y = 0,
                         y.max = 0.03,color = color_param$celltype[cell_type_list],
                         theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)

#plot link
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

temp <- which(p2g@metadata$geneSet$name == gene_name)
temp <- which(p2g$idxRNA == temp)
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
links_meta_data <- links_meta_data[links_meta_data$start >= start_site & links_meta_data$end <= end_site,]

p_p2g <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                 facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),xAixs.info = FALSE)

#plot
pdf(file = '/home/sunym/temp/temp.pdf',width = 16,height = 8)
p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.1) %>% insert_bottom(plot = p_p2g,height = 0.2) %>% insert_bottom(plot = p_gene_anno,height = 0.1)
dev.off()

#try ArchR
setwd('/home/sunym/temp')

gene_list <- c(gene_name)
gene_anno <- getGeneAnnotation(ArchRProj = Greenleaf_ATAC_ArchR)

p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g_loop <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = TRUE
)

names(human_SNC) <- human_SNC$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,geneSymbol = gene_list,groupBy = 'cell_type',
                      useGroups = c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC'),
                      features = human_SNC[DF_list[DF_list$group == 'human_specific',"human_SNC"]],
                      upstream = 100000,downstream = 100000,
                      geneAnnotation = gene_anno,
                      loops = p2g_loop,
                      sizes = c(10,1.5,1.5,1.5))

pdf(file = '/home/sunym/temp/test.pdf',width = 7,height = 4.5)
grid::grid.newpage()
grid::grid.draw(p$CSTB)
dev.off()

#yeah, transplotR worked!

# human vs macaque re write ---------------------------------------------------

## human -------------------------------------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

human_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
human_anno <- rtracklayer::as.data.frame(x = human_anno)

#load param
cell_type <- 'RG-1'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
gene_name <- 'CSTB'
extend_size <- 100000
resolution <- 1000
ylim_track <- 0.03

#load DF_list
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
DF_list <- readRDS(file = DF_list)

#load BW
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)

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

#get target site
TSS_site <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]

chrom <- as.character(TSS_site@seqnames)
start_site <- TSS_site@ranges@start - extend_size
end_site <- TSS_site@ranges@start + extend_size

#plot gene annotation
p_gene_anno <- trancriptVis(gtfFile = human_anno,collapse = TRUE,
                            Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                            posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                            addNormalArrow = FALSE,text.pos = 'middle',
                            newStyleArrow = FALSE,xAxis.info = TRUE)

#plot peak
temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
temp <- human_SNC
names(temp) <- temp$human_SNC
temp <- temp[DF_list[DF_list$group == 'human_specific',"human_SNC"]]
rtracklayer::export(object = temp,con = '/home/sunym/temp/SNC.bed',format = 'bed')

p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed','/home/sunym/temp/SNC.bed'),
                      chr = chrom,region.min = start_site,region.max = end_site,
                      show.legend = FALSE,fill = ggsci::pal_d3()(2),track.width = 0.3)

#plot track
p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = human_anno,
                         chr = chrom,region.min = start_site,region.max = end_site,
                         sample.order = cell_type_list,space.y = 0,
                         y.max = ylim_track,color = color_param$celltype[cell_type_list],
                         theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)

#plot link
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = resolution,
  returnLoops = FALSE
)

temp <- which(p2g@metadata$geneSet$name == gene_name)
temp <- which(p2g$idxRNA == temp)
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
links_meta_data <- links_meta_data[links_meta_data$start >= start_site & links_meta_data$end <= end_site,]

p_p2g <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                 facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),xAixs.info = FALSE)

#plot
pdf(file = './res/step_79_fig_221116/human_CSTB_trackplot.pdf',width = 8,height = 8)
p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.05) %>% insert_bottom(plot = p_p2g,height = 0.1) %>% insert_bottom(plot = p_gene_anno,height = 0.05)
dev.off()

## macaque -----------------------------------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

macaque_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_anno <- rtracklayer::as.data.frame(x = macaque_anno)

#load param
cell_type <- 'RG-1'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
gene_name <- 'CSTB'
extend_size <- 100000
resolution <- 1000
ylim_track <- 0.03

#load DF_list
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
DF_list <- readRDS(file = DF_list)

#load BW
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)

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

#get target site
TSS_site <- getPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]

chrom <- as.character(TSS_site@seqnames)
start_site <- TSS_site@ranges@start - extend_size
end_site <- TSS_site@ranges@start + extend_size

#plot gene annotation
p_gene_anno <- trancriptVis(gtfFile = macaque_anno,collapse = TRUE,
                            Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                            posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                            addNormalArrow = FALSE,text.pos = 'middle',
                            newStyleArrow = FALSE,xAxis.info = TRUE)

#plot peak
temp <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
temp <- macaque_SNC
names(temp) <- temp$human_SNC
temp <- temp[DF_list[DF_list$group == 'human_specific',"human_SNC"]]
rtracklayer::export(object = temp,con = '/home/sunym/temp/SNC.bed',format = 'bed')

p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed','/home/sunym/temp/SNC.bed'),
                      chr = chrom,region.min = start_site,region.max = end_site,
                      show.legend = FALSE,fill = ggsci::pal_d3()(2),track.width = 0.3)

#plot track
p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = macaque_anno,
                         chr = chrom,region.min = start_site,region.max = end_site,
                         sample.order = cell_type_list,space.y = 0,
                         y.max = ylim_track,color = color_param$celltype[cell_type_list],
                         theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)

#plot link
p2g <- getPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.45,
  resolution = resolution,
  returnLoops = FALSE
)

temp <- which(p2g@metadata$geneSet$name == gene_name)
temp <- which(p2g$idxRNA == temp)
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
links_meta_data <- links_meta_data[links_meta_data$start >= start_site & links_meta_data$end <= end_site,]

p_p2g <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                 facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),xAixs.info = FALSE)

#plot
pdf(file = './res/step_79_fig_221116/macaque_CSTB_trackplot.pdf',width = 8,height = 8)
p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.05) %>% insert_bottom(plot = p_p2g,height = 0.1) %>% insert_bottom(plot = p_gene_anno,height = 0.05)
dev.off()

# CSTB ---------------------------------------------------

## human -------------------------------------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

human_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
human_anno <- rtracklayer::as.data.frame(x = human_anno)

#load param
cell_type <- 'RG-1'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
gene_name <- 'CSTB'
extend_size <- 100000
resolution <- 1000
ylim_track <- 0.03

#load DF_list
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
DF_list <- readRDS(file = DF_list)

#load BW
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)

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

#get target site
TSS_site <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]

chrom <- as.character(TSS_site@seqnames)
start_site <- TSS_site@ranges@start - extend_size
end_site <- TSS_site@ranges@start + extend_size

#plot gene annotation
p_gene_anno <- trancriptVis(gtfFile = human_anno,collapse = TRUE,
                            Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                            posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                            addNormalArrow = FALSE,text.pos = 'middle',
                            newStyleArrow = FALSE,xAxis.info = TRUE)

#plot peak
temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
temp <- human_SNC
names(temp) <- temp$human_SNC
temp <- temp[DF_list[DF_list$group == 'human_specific',"human_SNC"]]
rtracklayer::export(object = temp,con = '/home/sunym/temp/SNC.bed',format = 'bed')

p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed','/home/sunym/temp/SNC.bed'),
                      chr = chrom,region.min = start_site,region.max = end_site,
                      show.legend = FALSE,fill = ggsci::pal_d3()(2),track.width = 0.3)

#plot track
p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = human_anno,
                         chr = chrom,region.min = start_site,region.max = end_site,
                         sample.order = cell_type_list,space.y = 0,
                         y.max = ylim_track,color = color_param$celltype[cell_type_list],
                         theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)

#plot link
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = resolution,
  returnLoops = FALSE
)

temp <- which(p2g@metadata$geneSet$name == gene_name)
temp <- which(p2g$idxRNA == temp)
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
links_meta_data <- links_meta_data[links_meta_data$start >= start_site & links_meta_data$end <= end_site,]

p_p2g <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                 facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),xAixs.info = FALSE)

#plot
pdf(file = './res/step_79_fig_221116/human_CSTB_trackplot.pdf',width = 8,height = 8)
p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.05) %>% insert_bottom(plot = p_p2g,height = 0.1) %>% insert_bottom(plot = p_gene_anno,height = 0.05)
dev.off()

## macaque -----------------------------------------------------------------
#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

macaque_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
macaque_anno <- rtracklayer::as.data.frame(x = macaque_anno)

#load param
cell_type <- 'RG-1'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')
gene_name <- 'CSTB'
extend_size <- 100000
resolution <- 1000
ylim_track <- 0.03

#load DF_list
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
DF_list <- readRDS(file = DF_list)

#load BW
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list)

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

#get target site
TSS_site <- getPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)
TSS_site <- TSS_site@metadata$geneSet[which(TSS_site@metadata$geneSet$name == gene_name)]

chrom <- as.character(TSS_site@seqnames)
start_site <- TSS_site@ranges@start - extend_size
end_site <- TSS_site@ranges@start + extend_size

#plot gene annotation
macaque_anno <- macaque_anno[which(macaque_anno$gene_name != 'RRP1'),]
p_gene_anno <- trancriptVis(gtfFile = macaque_anno,collapse = TRUE,
                            Chr = gsub(pattern = 'chr',replacement = '',x = chrom,fixed = TRUE),
                            posStart = start_site,posEnd = end_site,textLabel = 'gene_name',
                            addNormalArrow = FALSE,text.pos = 'middle',
                            newStyleArrow = FALSE,xAxis.info = TRUE)

#plot peak
temp <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
rtracklayer::export(object = temp,con = '/home/sunym/temp/peak.bed',format = 'bed')
temp <- macaque_SNC
names(temp) <- temp$human_SNC
temp <- temp[DF_list[DF_list$group == 'human_specific',"human_SNC"]]
rtracklayer::export(object = temp,con = '/home/sunym/temp/SNC.bed',format = 'bed')

p_peak_anno <- bedVis(bdFile = c('/home/sunym/temp/peak.bed','/home/sunym/temp/SNC.bed'),
                      chr = chrom,region.min = start_site,region.max = end_site,
                      show.legend = FALSE,fill = ggsci::pal_d3()(2),track.width = 0.3)

#plot track
p_track_plot <- trackVis(bWData = cell_type_bw,gtf.file = macaque_anno,
                         chr = chrom,region.min = start_site,region.max = end_site,
                         sample.order = cell_type_list,space.y = 0,
                         y.max = ylim_track,color = color_param$celltype[cell_type_list],
                         theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)

#plot link
links_meta_data <- data.frame(chr = chrom,start = start_site - 3000,end = start_site + 3000,cor = 1,group = 'group1')

p_p2g <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                 facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),xAixs.info = FALSE)

#plot
pdf(file = './res/step_79_fig_221116/macaque_CSTB_trackplot.pdf',width = 8,height = 8)
p_track_plot %>% insert_bottom(plot = p_peak_anno,height = 0.05) %>% insert_bottom(plot = p_p2g,height = 0.1) %>% insert_bottom(plot = p_gene_anno,height = 0.05)
dev.off()

## RNA expression ----------------------------------------------------------
#load data
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')
Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
Idents(macaque_multiome_Seurat) <- 'cell_type'
Idents(Greenleaf_RNA_Seurat) <- 'ReAnno_celltype'

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = 'human CSTB') + NoLegend()

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = gene_name,pt.size = 0.1,
                  slot = 'data',label = FALSE,cols = ArchRPalettes$whitePurple[c(-1,-2)],
                  order = TRUE) + 
  theme(aspect.ratio = 1) + NoAxes() + 
  labs(title = 'macaque CSTB') + NoLegend()

pdf(file = './res/step_79_fig_221116/CSTB_human_macaque_RNA_express_featureplot.pdf',width = 5,height = 10)
p1+p2+plot_layout(ncol = 1)
dev.off()

## see motif changes -------------------------------------------------------
#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
sequence_list <- readRDS(file = './res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')

#load affinity change
affinity_difference <- readRDS(file = './res/step_77_fig_221108/affinity_difference_data_frame.rds')
human_fimo <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
motif_info <- read.csv(file = './data/reference/motif_DB/CISBP/human_221108/TF_information.csv')

#SNC related to CSTB
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')

#find human SNC related to CSTB
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

temp <- which(p2g@metadata$geneSet$name == gene_name)
temp <- which(p2g$idxRNA == temp)
temp <- p2g$idxATAC[temp]
temp <- p2g@metadata$peakSet[temp]

names(human_SNC) <- human_SNC$human_SNC
SNC_list <- countOverlaps(query = human_SNC,subject = temp)
SNC_list <- names(SNC_list)[SNC_list > 0]

#see the affinity change
affinity_difference$human_SNC <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))
affinity_difference <- affinity_difference[affinity_difference$human_SNC %in% SNC_list & affinity_difference$cell_type == cell_type,]

#motif
motif_1 <- list.files(path = './data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs')
motif_1 <- motif_1[grep(pattern = 'M07691_2.00',x = motif_1,fixed = TRUE)]
motif_1 <- paste('./data/reference/motif_DB/CISBP/human_221108/pwms_all_motifs',motif_1,sep = '/')
motif_1 <- read_cisbp(file = motif_1)
motif_1@name <- 'M07691_2.00'

pdf(file = './res/step_79_fig_221116/M07691_2.00_motif_view.pdf',width = 8,height = 2.5)
view_motifs(motifs = motif_1,show.names = TRUE,names.pos = 'right')
dev.off()

#get sequence
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

human_seq <- sequence_list$human_seq[[SNC_list]]
macaque_seq <- sequence_list$macaque_seq[[SNC_list]]

human_fimo[human_fimo$motif_id == motif_1@name & human_fimo$sequence_name == SNC_list,]

pairwiseAlignment(pattern = reverseComplement(human_seq),
                  subject = as('CCTGCGGGAGCTGCGTCCCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

macaque_fimo[macaque_fimo$motif_id == motif_1@name & macaque_fimo$sequence_name == SNC_list,]

pairwiseAlignment(pattern = reverseComplement(macaque_seq),
                  subject = as('CCTGCGGGAGCTGCGTCCCAG','DNAString'),
                  type="global",substitutionMatrix = NSM,gapOpening=10,
                  gapExtension=4,scoreOnly = FALSE)

#check motif information
ii <- motif_info[motif_info$Motif_ID == motif_1@name,]

# filter RG-1 gene --------------------------------------------------------
cell_type <- 'RG-1'
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4','OPC')

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
sequence_list <- readRDS(file = './res/step_77_fig_221108/human_SNC_extented_aligned_sequence.rds')

affinity_difference <- readRDS(file = './res/step_77_fig_221108/affinity_difference_data_frame.rds')
human_fimo <- read.table(file = './res/step_77_fig_221108/human_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
macaque_fimo <- read.table(file = './res/step_77_fig_221108/macaque_SNC_seq_fimo_out/fimo.tsv',sep = '\t',header = TRUE)
motif_info <- read.csv(file = './data/reference/motif_DB/CISBP/human_221108/TF_information.csv')

human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')

#load gene list
cell_type_dot <- gsub(pattern = '-',replacement = '.',fixed = TRUE,x = cell_type)
gene_list <- list.files(path = './res/step_79_fig_221112/up_gene_regulated_by_human_specific_SNC')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/up_gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(file = gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#modify affinity difference
affinity_difference$human_SNC <- unlist(base::lapply(X = affinity_difference$set,FUN = function(x){
  temp <- strsplit(x = x,split = '#')
  return(temp[[1]][2])
}))

#filter gene list
gene_list <- base::lapply(X = gene_list,FUN = function(x){
  #filter human SNC overlap with p2g
  temp <- which(p2g@metadata$geneSet$name == x)
  temp <- which(p2g$idxRNA %in% temp)
  temp <- p2g$idxATAC[temp]
  temp <- p2g@metadata$peakSet[temp]
  
  names(human_SNC) <- human_SNC$human_SNC
  SNC_list <- countOverlaps(query = human_SNC,subject = temp)
  SNC_list <- names(SNC_list)[SNC_list > 0]
  
  #see the affinity change
  temp <- affinity_difference[affinity_difference$human_SNC %in% SNC_list,]
  
  #return
  if(nrow(temp) > 0){
    return(x)
  }else{
    return(NA)
  }
})

gene_list <- unlist(gene_list)
gene_list <- gene_list[!is.na(gene_list)]
