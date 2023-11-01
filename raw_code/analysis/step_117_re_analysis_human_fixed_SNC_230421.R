#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: re_analysis_human_fixed_SNC                                     ##
## Data: 2023.04.21                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.0/')
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
chat <- Init_chat_session(model = 'gpt-4-0314')

# pre-process modern human fixed SNC --------------------------------------
SNC_GRanges <- readRDS(file = './data/public/The_cis_regulatory_effects_of_modern_human_specific_variants/processed_data/SNC_GRanges.rds')

#check the base in table
table(BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,name = SNC_GRanges,as.character = TRUE) == SNC_GRanges$modern_nucleotide)
#basicly correct

#re create the SNC_GRanges
MH_SNC <- as.character(SNC_GRanges)
names(MH_SNC) <- NULL
MH_SNC <- as(MH_SNC,'GRanges')
names(MH_SNC) <- as.character(MH_SNC)

names(SNC_GRanges) <- as.character(SNC_GRanges)
table(names(MH_SNC) == names(SNC_GRanges))

#assign base character
MH_SNC@elementMetadata$modern_nucleotide <- SNC_GRanges@elementMetadata$modern_nucleotide
MH_SNC@elementMetadata$archaic_nucleotide <- SNC_GRanges@elementMetadata$archaic_nucleotide
MH_SNC@elementMetadata$hg38_nucleotide <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,names = MH_SNC,as.character = TRUE)

#liftover to Mmul_10
temp <- rtracklayer::import.chain(con = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac10.over.chain')
MH_SNC_macaque <- rtracklayer::liftOver(x = MH_SNC,chain = temp)
MH_SNC_macaque <- base::unlist(MH_SNC_macaque)
table(duplicated(names(MH_SNC_macaque)))

MH_SNC <- MH_SNC[names(MH_SNC) %in% names(MH_SNC_macaque)]
MH_SNC_macaque <- MH_SNC_macaque[names(MH_SNC)]
MH_SNC@elementMetadata$human_coord <- names(MH_SNC)
MH_SNC@elementMetadata$macaque_coord <- as.character(MH_SNC_macaque)

#get macaque nucleotide
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

extend_human <- paste0(as.character(MH_SNC@seqnames),":",MH_SNC@ranges@start - 250,"-",MH_SNC@ranges@start + 250)
extend_human <- as(extend_human,"GRanges")

extend_macaque <- paste0(as.character(MH_SNC_macaque@seqnames),":",MH_SNC_macaque@ranges@start - 250,"-",MH_SNC_macaque@ranges@start + 250)
extend_macaque <- as(extend_macaque,"GRanges")

macaque_base_list <- base::lapply(X = 1:length(extend_human),FUN = function(x){
  human_seq <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,names = extend_human[x])
  macaque_seq <- BSgenome::getSeq(x = BSgenome.Mmulatta.UCSC.rheMac10,names = extend_macaque[x])
  
  #alignment
  default_score <- pairwiseAlignment(pattern = human_seq,
                                     subject = macaque_seq,
                                     type="global",
                                     substitutionMatrix = NSM,
                                     gapOpening=10,gapExtension=4,
                                     scoreOnly = TRUE)
  reverse_score <- pairwiseAlignment(pattern = human_seq,
                                     subject = Biostrings::reverse(macaque_seq),
                                     type="global",
                                     substitutionMatrix = NSM,
                                     gapOpening=10,gapExtension=4,
                                     scoreOnly = TRUE)
  complement_score <- pairwiseAlignment(pattern = human_seq,
                                        subject = Biostrings::complement(macaque_seq),
                                        type="global",
                                        substitutionMatrix = NSM,
                                        gapOpening=10,gapExtension=4,
                                        scoreOnly = TRUE)
  reverseComplement_score <- pairwiseAlignment(pattern = human_seq,
                                               subject = Biostrings::reverseComplement(macaque_seq),
                                               type="global",
                                               substitutionMatrix = NSM,
                                               gapOpening=10,gapExtension=4,
                                               scoreOnly = TRUE)
  
  #get the correct macaque seq
  idx <- which.max(c(default_score,reverse_score,complement_score,reverseComplement_score))
  if(idx == 1){
    macaque_seq = macaque_seq
  }else if(idx == 2){
    macaque_seq = Biostrings::reverse(macaque_seq)
  }else if(idx == 3){
    macaque_seq = Biostrings::complement(macaque_seq)
  }else if(idx == 4){
    macaque_seq = Biostrings::reverseComplement(macaque_seq)
  }else{
    stop('something wrong!')
  }
  
  #return
  temp <- as.character(macaque_seq[[1]][251])
  return(temp)
})
macaque_base_list <- unlist(macaque_base_list)

MH_SNC@elementMetadata$Mmul_10_nucleotide <- macaque_base_list

#save MH_SNC
MH_SNC <- MH_SNC[MH_SNC$hg38_nucleotide != MH_SNC$Mmul_10_nucleotide]
MH_SNC@elementMetadata <- MH_SNC@elementMetadata[,c("human_coord","macaque_coord","modern_nucleotide","archaic_nucleotide","hg38_nucleotide","Mmul_10_nucleotide")]
names(MH_SNC) <- MH_SNC$human_coord
saveRDS(object = MH_SNC,file = './res/step_117_fig_230421/modern_human_fixed_SNC_GRanges.rds')

# get human seq and macaque seq -------------------------------------------
#SNC +- 30bp
MH_SNC <- readRDS(file = './res/step_117_fig_230421/modern_human_fixed_SNC_GRanges.rds')

extend_human <- base::lapply(X = MH_SNC$human_coord,FUN = function(x){
  SNC_site <- strsplit(x = x,split = ':')
  start_site <- as.numeric(SNC_site[[1]][2]) - 30
  end_site <- as.numeric(SNC_site[[1]][2]) + 30
  return(paste0(SNC_site[[1]][1],':',start_site,"-",end_site))
})
extend_human <- as(unlist(extend_human),"GRanges")
names(extend_human) <- MH_SNC$human_coord

extend_macaque <- base::lapply(X = MH_SNC$macaque_coord,FUN = function(x){
  SNC_site <- strsplit(x = x,split = ':')
  start_site <- as.numeric(SNC_site[[1]][2]) - 30
  end_site <- as.numeric(SNC_site[[1]][2]) + 30
  return(paste0(SNC_site[[1]][1],':',start_site,"-",end_site))
})
extend_macaque <- as(unlist(extend_macaque),"GRanges")
names(extend_macaque) <- MH_SNC$macaque_coord

#get seq
extend_human <- BSgenome::getSeq(x = BSgenome.Hsapiens.UCSC.hg38,names = extend_human)
extend_macaque <- BSgenome::getSeq(x = BSgenome.Mmulatta.UCSC.rheMac10,names = extend_macaque)

writeXStringSet(x = extend_human,filepath = './res/step_117_fig_230421/human_extended_30bp.fasta',format = 'fasta')
writeXStringSet(x = extend_macaque,filepath = './res/step_117_fig_230421/macaque_extended_30bp.fasta',format = 'fasta')

# get motif file and turn into meme format --------------------------------
motif_list <- list.files('./data/reference/motif_DB/JASPAR2022/pfms_meme')
motif_list <- paste('./data/reference/motif_DB/JASPAR2022/pfms_meme',motif_list,sep = '/')

motif_file <- base::lapply(X = motif_list,FUN = function(x){
  temp <- universalmotif::read_meme(file = x)
  return(temp)
})

universalmotif::write_meme(motifs = motif_file,file = './res/step_117_fig_230421/JASPAR2022_motifs.meme')


# modify fimo results -----------------------------------------------------
#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

MH_SNC <- readRDS(file = './res/step_117_fig_230421/modern_human_fixed_SNC_GRanges.rds')
consensus_peakset <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')

human_fimo <- read.table(file = './res/step_117_fig_230421/human_SNC_fimo_out/fimo.tsv',header = TRUE)
macaque_fimo <- read.table(file = './res/step_117_fig_230421/macaque_SNC_fimo_out/fimo.tsv',header = TRUE)

#filter fimo out
idx <- base::lapply(X = 1:nrow(human_fimo),FUN = function(x){
  if(human_fimo[x,"start"] <= 31 & human_fimo[x,"stop"] >= 31){
    return(x)
  }else{
    return(NULL)
  }
})
idx <- base::unlist(idx)
human_fimo <- human_fimo[idx,]

idx <- base::lapply(X = 1:nrow(macaque_fimo),FUN = function(x){
  if(macaque_fimo[x,"start"] <= 31 & macaque_fimo[x,"stop"] >= 31){
    return(x)
  }else{
    return(NULL)
  }
})
idx <- base::unlist(idx)
macaque_fimo <- macaque_fimo[idx,]

#overlaped consensus peaks
#human
human_fimo$consensus_peakset <- NULL
temp_query <- as(human_fimo$sequence_name,'GRanges')
temp_subject <- as(consensus_peakset$human_coord,'GRanges')
names(temp_subject) <- names(consensus_peakset)

idx <- findOverlaps(query = temp_query,subject = temp_subject)
human_fimo[queryHits(idx),'consensus_peakset'] <- names(temp_subject)[subjectHits(idx)]

human_fimo <- human_fimo[!is.na(human_fimo$consensus_peakset),]

#macaque
macaque_fimo$consensus_peakset <- NULL
temp_query <- as(macaque_fimo$sequence_name,'GRanges')
temp_subject <- as(consensus_peakset$macaque_coord,'GRanges')
names(temp_subject) <- names(consensus_peakset)

idx <- findOverlaps(query = temp_query,subject = temp_subject)
macaque_fimo[queryHits(idx),'consensus_peakset'] <- names(temp_subject)[subjectHits(idx)]

macaque_fimo <- macaque_fimo[!(is.na(macaque_fimo$consensus_peakset)),]

#linked gene
human_P2G <- getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = 1,returnLoops = FALSE)
human_gene_list <- base::lapply(X = human_fimo$consensus_peakset,FUN = function(x){
  #get consensus peak coord
  peak_coord <- consensus_peakset[x]$human_coord
  peak_coord <- as(peak_coord,'GRanges')
  
  #overlapped peak
  idx <- countOverlaps(query = human_P2G@metadata$peakSet,subject = peak_coord)
  idx <- which(idx > 0)
  
  if(length(idx) == 0){
    return('')
  }
  
  #overlapped linkage
  idx <- which(human_P2G$idxATAC %in% idx)
  
  if(length(idx) == 0){
    return('')
  }
  
  #regulated genes
  idx <- human_P2G$idxRNA[idx]
  idx <- human_P2G@metadata$geneSet$name[idx]
  return(paste(idx,collapse = ','))
})
human_fimo$affected_human_genes <- unlist(human_gene_list)

macaque_P2G <- getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = 1,returnLoops = FALSE)
macaque_gene_list <- base::lapply(X = macaque_fimo$consensus_peakset,FUN = function(x){
  #get consensus peak coord
  peak_coord <- consensus_peakset[x]$macaque_coord
  peak_coord <- as(peak_coord,'GRanges')
  
  #overlapped peak
  idx <- countOverlaps(query = macaque_P2G@metadata$peakSet,subject = peak_coord)
  idx <- which(idx > 0)
  
  if(length(idx) == 0){
    return('')
  }
  
  #overlapped linkage
  idx <- which(macaque_P2G$idxATAC %in% idx)
  
  if(length(idx) == 0){
    return('')
  }
  
  #regulated genes
  idx <- macaque_P2G$idxRNA[idx]
  idx <- macaque_P2G@metadata$geneSet$name[idx]
  return(paste(idx,collapse = ','))
})
macaque_fimo$affected_macaque_genes <- unlist(macaque_gene_list)

#save data
saveRDS(object = human_fimo,file = './res/step_117_fig_230421/modified_human_fimo_out.rds')
saveRDS(object = macaque_fimo,file = './res/step_117_fig_230421/modified_macaque_fimo_out.rds')

# RG HGARs caused by MH SNC -----------------------------------------------
# something wrong here, may be the plot function
#load data
cell_type <- 'RG-1'
resolution <- 500
ylim_track <- 0.03
cell_type_list <- c('RG-1','RG-2','IP','Ex-1','Ex-2','Ex-3','Ex-4')
human_fimo <- readRDS(file = './res/step_117_fig_230421/modified_human_fimo_out.rds')
macaque_fimo <- readRDS(file = './res/step_117_fig_230421/modified_macaque_fimo_out.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

consensus_peakset <- readRDS(file = './res/step_118_fig_230422/consensus_peakset_GRanges.rds')

human_gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Homo_sapiens.GRCh38.105.gtf',format = 'gtf')
macaque_gene_anno <- rtracklayer::import(con = './data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf',format = 'gtf')
human_gene_anno <- rtracklayer::as.data.frame(human_gene_anno)
macaque_gene_anno <- rtracklayer::as.data.frame(macaque_gene_anno)

MH_SNC <- readRDS(file = './res/step_117_fig_230421/modern_human_fixed_SNC_GRanges.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_230320.rds')

#get p2g
human_P2G <- getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = 1,returnLoops = FALSE)
macaque_P2G <- getPeak2GeneLinks(ArchRProj = macaque_multiome_ArchR,resolution = 1,returnLoops = FALSE)

#load DAP
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
DAP_list <- list.files(path = './res/step_118_fig_230422/Wilcox')
DAP_list <- DAP_list[grep(pattern = cell_type_dot,x = DAP_list,fixed = TRUE)]
DAP_list <- paste('./res/step_118_fig_230422/Wilcox',DAP_list,sep = '/')
DAP_list <- readRDS(DAP_list)

DAP_list <- DAP_list[DAP_list$group == 'HGARs',]

#add overlapped DAP to human fimo
human_fimo$DAP <- NULL
human_fimo[which(human_fimo$consensus_peakset %in% DAP_list$peak),'DAP'] <- 'HGARs'
human_fimo <- human_fimo[which(human_fimo$affected_human_genes != '' & human_fimo$DAP == 'HGARs'),]

temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
human_fimo <- human_fimo[which(countOverlaps(query = as(human_fimo$sequence_name,'GRanges'),subject = temp) > 0),]

#gene_list
gene_list <- base::lapply(X = human_fimo$affected_human_genes,FUN = function(x){
  temp <- strsplit(x = x,split = ',')
  return(temp[[1]])
})
gene_list <- unique(unlist(gene_list))

#load bw
cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)
cell_type_bw <- list.files(path = './res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type')
cell_type_bw <- unlist(base::lapply(X = cell_type_list_dot,FUN = function(x){
  temp <- cell_type_bw[grep(pattern = x,x = cell_type_bw,fixed = TRUE)]
  return(temp)
}))
cell_type_bw <- paste('./res/step_75_fig_221103/norm_by_ReadsInPeaks/human/GroupBigWigs/cell_type',cell_type_bw,sep = '/')
cell_type_bw <- loadBigWig(bwFile = cell_type_bw)
cell_type_bw$fileName <- gsub(pattern = '-TileSize-25-normMethod-ReadsInPeaks-ArchR',replacement = '',x = cell_type_bw$fileName,fixed = TRUE)
cell_type_bw$fileName <- gsub(pattern = '.',replacement = '-',x = cell_type_bw$fileName,fixed = TRUE)
gc()
human_bw <- cell_type_bw

cell_type_list_dot <- gsub(pattern = '-',replacement = '.',x = cell_type_list,fixed = TRUE)
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
macaque_bw <- cell_type_bw

rm(cell_type_bw)
gc()

#trackplot
for (gene_symbol in gene_list) {
  
  #get affected consensus peakset
  human_affected_peak <- human_fimo[grep(pattern = gene_symbol,x = human_fimo$affected_human_genes,fixed = TRUE),"consensus_peakset"]
  human_affected_peak <- unique(human_affected_peak)
  
  macaque_affected_peak <- consensus_peakset[which(consensus_peakset$human_coord %in% human_affected_peak)]$macaque_coord
  
  #get TSS
  temp_chrom <- as(human_affected_peak,'GRanges')
  temp_chrom <- as.character(unique(temp_chrom@seqnames))
  if(length(temp_chrom) > 1){
    stop('multi-chrom detected!')
  }
  human_TSS <- human_P2G@metadata$geneSet
  human_TSS <- human_TSS[which(human_TSS@seqnames == temp_chrom & human_TSS$name == gene_symbol)]
  
  temp_chrom <- as(macaque_affected_peak,'GRanges')
  temp_chrom <- as.character(unique(temp_chrom@seqnames))
  if(length(temp_chrom) > 1){
    stop('multi-chrom detected!')
  }
  macaque_TSS <- macaque_P2G@metadata$geneSet
  macaque_TSS <- macaque_TSS[which(macaque_TSS@seqnames == temp_chrom & macaque_TSS$name == gene_symbol)]
  
  #modify ranges
  temp <- as(human_affected_peak,'GRanges')
  temp <- append(temp,human_TSS)
  human_start_site <- min(temp@ranges@start) - 10000
  human_end_site <- max(temp@ranges@start + temp@ranges@width - 1) + 10000
  
  temp <- as(macaque_affected_peak,'GRanges')
  temp <- append(temp,macaque_TSS)
  macaque_start_site <- min(temp@ranges@start) - 10000
  macaque_end_site <- max(temp@ranges@start + temp@ranges@width - 1) + 10000
  
  #gene plot
  human_gene_plot <- trancriptVis(gtfFile = human_gene_anno,collapse = TRUE,gene = gene_symbol,
                                  Chr = gsub(pattern = 'chr',replacement = '',x = as.character(human_TSS@seqnames),fixed = TRUE),
                                  posStart = human_start_site,posEnd = human_end_site,textLabel = 'gene_name',
                                  addNormalArrow = TRUE,text.pos = 'middle',newStyleArrow = FALSE,xAxis.info = TRUE)
  macaque_gene_plot <- trancriptVis(gtfFile = macaque_gene_anno,collapse = TRUE,gene = gene_symbol,
                                    Chr = gsub(pattern = 'chr',replacement = '',x = as.character(macaque_TSS@seqnames),fixed = TRUE),
                                    posStart = macaque_start_site,posEnd = macaque_end_site,textLabel = 'gene_name',
                                    addNormalArrow = TRUE,text.pos = 'middle',newStyleArrow = FALSE,xAxis.info = TRUE)
  
  #peak plot
  temp <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
  non_overlapped_peak <- temp[which(countOverlaps(query = temp,subject = as(MH_SNC$human_coord,'GRanges')) == 0)]
  overlapped_peak <- temp[which(countOverlaps(query = temp,subject = as(MH_SNC$human_coord,'GRanges')) > 0)]
  rtracklayer::export(object = non_overlapped_peak,con = '/home/sunym/temp/non_overlapped_peak.bed',format = 'bed')
  rtracklayer::export(object = overlapped_peak,con = '/home/sunym/temp/overlapped_peak.bed',format = 'bed')
  
  human_peak_plot <- bedVis(bdFile = c('/home/sunym/temp/non_overlapped_peak.bed','/home/sunym/temp/overlapped_peak.bed'),
                            chr = as.character(human_TSS@seqnames),region.min = human_start_site,region.max = human_end_site,
                            show.legend = FALSE,fill = ggsci::pal_d3()(2),track.width = 0.15,collapse = TRUE)
  
  temp <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
  non_overlapped_peak <- temp[which(countOverlaps(query = temp,subject = as(MH_SNC$macaque_coord,'GRanges')) == 0)]
  overlapped_peak <- temp[which(countOverlaps(query = temp,subject = as(MH_SNC$macaque_coord,'GRanges')) > 0)]
  rtracklayer::export(object = non_overlapped_peak,con = '/home/sunym/temp/non_overlapped_peak.bed',format = 'bed')
  rtracklayer::export(object = overlapped_peak,con = '/home/sunym/temp/overlapped_peak.bed',format = 'bed')
  
  macaque_peak_plot <- bedVis(bdFile = c('/home/sunym/temp/non_overlapped_peak.bed','/home/sunym/temp/overlapped_peak.bed'),
                              chr = as.character(macaque_TSS@seqnames),region.min = macaque_start_site,region.max = macaque_end_site,
                              show.legend = FALSE,fill = ggsci::pal_d3()(2),track.width = 0.15,collapse = TRUE)
  
  #track plot
  human_track_plot <- trackVis(bWData = human_bw,gtf.file = human_gene_anno,
                               chr = as.character(human_TSS@seqnames),region.min = human_start_site,region.max = human_end_site,
                               sample.order = cell_type_list,space.y = 0,
                               y.max = ylim_track,color = color_param$celltype[cell_type_list],
                               theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)
  
  macaque_track_plot <- trackVis(bWData = macaque_bw,gtf.file = macaque_gene_anno,
                                 chr = as.character(macaque_TSS@seqnames),region.min = macaque_start_site,region.max = macaque_end_site,
                                 sample.order = cell_type_list,space.y = 0,
                                 y.max = ylim_track,color = color_param$celltype[cell_type_list],
                                 theme = 'bw',xAxis.info = FALSE,yAxis.info = FALSE,new.yaxis = TRUE)
  
  #link plot
  p2g <- getPeak2GeneLinks(
    ArchRProj = Greenleaf_ATAC_ArchR,
    corCutOff = 0.45,
    resolution = resolution,
    returnLoops = FALSE
  )
  
  temp <- which(p2g@metadata$geneSet$name == gene_symbol)
  temp <- which(p2g$idxRNA == temp)
  idx_ATAC <- p2g$idxATAC[temp]
  idx_ATAC <- p2g@metadata$peakSet[idx_ATAC]
  idx_ATAC <- rtracklayer::as.data.frame(idx_ATAC)
  idx_ATAC <- round((idx_ATAC$start + idx_ATAC$end) / 2)
  links_meta_data <- data.frame(chr = as.character(human_TSS@seqnames),start = NA,end = NA,cor = c(p2g$Correlation[temp],0,1),group = 'group1')
  for (i in 1:length(idx_ATAC)) {
    if(human_TSS@ranges@start <= idx_ATAC[i]){
      links_meta_data[i,'start'] <- human_TSS@ranges@start
      links_meta_data[i,"end"] <- idx_ATAC[i]
    }else{
      links_meta_data[i,"start"] <- idx_ATAC[i]
      links_meta_data[i,'end'] <- human_TSS@ranges@start
    }
  }
  human_link <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                        facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),
                        xAixs.info = FALSE,curvature = 0.1)
  
  p2g <- getPeak2GeneLinks(
    ArchRProj = macaque_multiome_ArchR,
    corCutOff = 0.45,
    resolution = resolution,
    returnLoops = FALSE
  )
  
  temp <- which(p2g@metadata$geneSet$name == gene_symbol)
  temp <- which(p2g$idxRNA == temp)
  idx_ATAC <- p2g$idxATAC[temp]
  idx_ATAC <- p2g@metadata$peakSet[idx_ATAC]
  idx_ATAC <- rtracklayer::as.data.frame(idx_ATAC)
  idx_ATAC <- round((idx_ATAC$start + idx_ATAC$end) / 2)
  links_meta_data <- data.frame(chr = as.character(macaque_TSS@seqnames),start = NA,end = NA,cor = c(p2g$Correlation[temp],0,1),group = 'group1')
  for (i in 1:length(idx_ATAC)) {
    if(macaque_TSS@ranges@start <= idx_ATAC[i]){
      links_meta_data[i,'start'] <- macaque_TSS@ranges@start
      links_meta_data[i,"end"] <- idx_ATAC[i]
    }else{
      links_meta_data[i,"start"] <- idx_ATAC[i]
      links_meta_data[i,'end'] <- macaque_TSS@ranges@start
    }
  }
  macaque_link <- linkVis(linkData = links_meta_data,start = 'start',end = 'end',
                          facet = FALSE,link.aescolor = 'cor',link.color = c('#FFD9CBFF','#FF0000FF'),
                          xAixs.info = FALSE,curvature = 0.1)
  
  #plot
  human_track_plot %>% insert_bottom(plot = human_peak_plot,height = 0.05) %>% insert_bottom(plot = human_link,height = 0.1) %>% insert_bottom(plot = human_gene_plot,height = 0.05)
  macaque_track_plot %>% insert_bottom(plot = macaque_peak_plot,height = 0.05) %>% insert_bottom(plot = macaque_link,height = 0.1) %>% insert_bottom(plot = macaque_gene_plot,height = 0.05)
}
