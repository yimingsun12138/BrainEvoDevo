#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: RG_1 HGAR overlap with public human gain enhancer               ##
## Data: 2023.02.16                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

#public data usage: Evolutionary changes in promoter and enhancer activity during human corticogenesis

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#general setting
setwd('/content/data/sunym/project/Brain/')
.libPaths('/content/data/sunym/software/R_lib/R_4.2.2/')
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
library(biomaRt)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

# export species homologous peaks -------------------------------------------
peak_list <- readRDS(file = './res/step_73_fig_221102/Brain_ATAC_peak.rds')

peak_human <- peak_list$human
peak_macaque <- peak_list$macaque
peak_mouse <- peak_list$mouse

rtracklayer::export.bed(object = peak_human,con = '/home/sunym/trash/human.bed')
rtracklayer::export.bed(object = peak_macaque,con = '/home/sunym/trash/macaque.bed')
rtracklayer::export.bed(object = peak_mouse,con = '/home/sunym/trash/mouse.bed')

#human
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/human.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/hg38/hg38.chrom.sizes',
              './res/step_104_fig_230216/export_bigBed/human_homologous_peak.bigBed')
system(command = char)

#macaque
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/macaque.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/rheMac10/rheMac10.chrom.sizes',
              './res/step_104_fig_230216/export_bigBed/macaque_homologous_peak.bigBed')
system(command = char)

#mouse
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/mouse.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/mm10/mm10.chrom.sizes',
              './res/step_104_fig_230216/export_bigBed/mouse_homologous_peak.bigBed')
system(command = char)

# export public H3K27ac peak ----------------------------------------------

#merge peakset
human_rep_1 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1554672_Hu_12Fpcw_H3K27ac_rep1_regions.bed')
human_rep_2 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1554673_Hu_12Fpcw_H3K27ac_rep2_regions.bed')

peak_human <- my_bedtools_merge(peakset_x = human_rep_1,peakset_y = human_rep_2,
                                bedtools_path = '/content/data/sunym/software/bedtools2/bedtools',
                                d = 0,bedtools_param = NULL,tmp_path = '/home/sunym/trash')
export.bed(object = peak_human,con = '/home/sunym/temp/human.bed')

macaque_rep_1 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1489728_Rh_12Fpcw_H3K27ac_rep1_regions.bed')
macaque_rep_2 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1489729_Rh_12Fpcw_H3K27ac_rep2_regions.bed')

peak_macaque <- my_bedtools_merge(peakset_x = macaque_rep_1,peakset_y = macaque_rep_2,
                                  bedtools_path = '/content/data/sunym/software/bedtools2/bedtools',
                                  d = 0,bedtools_param = NULL,tmp_path = '/home/sunym/trash')
export.bed(object = peak_macaque,con = '/home/sunym/temp/macaque.bed')

mouse_rep_1 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1489707_Mm_e14_H3K27ac_rep1_regions.bed')
mouse_rep_2 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1489708_Mm_e14_H3K27ac_rep2_regions.bed')

peak_mouse <- my_bedtools_merge(peakset_x = mouse_rep_1,peakset_y = mouse_rep_2,
                                bedtools_path = '/content/data/sunym/software/bedtools2/bedtools',
                                d = 0,bedtools_param = NULL,tmp_path = '/home/sunym/trash')
export.bed(object = peak_mouse,con = '/home/sunym/temp/mouse_e14.bed')

mouse_rep_1 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1489713_Mm_e17F_H3K27ac_rep1_regions.bed')
mouse_rep_2 <- rtracklayer::import.bed(con = './res/step_104_fig_230216/public_data/GSM1489714_Mm_e17F_H3K27ac_rep2_regions.bed')

peak_mouse <- my_bedtools_merge(peakset_x = mouse_rep_1,peakset_y = mouse_rep_2,
                                bedtools_path = '/content/data/sunym/software/bedtools2/bedtools',
                                d = 0,bedtools_param = NULL,tmp_path = '/home/sunym/trash')
export.bed(object = peak_mouse,con = '/home/sunym/temp/mouse_e17.bed')

#convert to bigbed
#human
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/temp/human.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/hg19/hg19.chrom.sizes',
              './res/step_104_fig_230216/export_bigBed/public_human_H3K27ac.bigBed')
system(command = char)

#macaque
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/temp/macaque.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/rheMac2/rheMac2.chrom.sizes',
              './res/step_104_fig_230216/export_bigBed/public_macaque_H3K27ac.bigBed')
system(command = char)

#mouse e14
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/temp/mouse_e14.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/mm9/mm9.chrom.sizes',
              './res/step_104_fig_230216/export_bigBed/public_mouse_e14_H3K27ac.bigBed')
system(command = char)

#mouse e17
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/temp/mouse_e17.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/mm9/mm9.chrom.sizes',
              './res/step_104_fig_230216/export_bigBed/public_mouse_e17_H3K27ac.bigBed')
system(command = char)

# public results module 10 ECM related gene -------------------------------
gene_list <- c('ENSG00000139567','ENSG00000116962','ENSG00000163520','ENSG00000047648','ENSG00000184371','ENSG00000108821',
               'ENSG00000134013','ENSG00000107562','ENSG00000050555','ENSG00000140092','ENSG00000125810','ENSG00000106991',
               'ENSG00000146648','ENSG00000204291','ENSG00000087303','ENSG00000162733','ENSG00000197467','ENSG00000179776',
               'ENSG00000137809','ENSG00000148180','ENSG00000126561','ENSG00000140416','ENSG00000198959','ENSG00000100767',
               'ENSG00000134871','ENSG00000143382','ENSG00000187498','ENSG00000125378','ENSG00000182492','ENSG00000132561',
               'ENSG00000109625','ENSG00000164692','ENSG00000185585','ENSG00000104324','ENSG00000069702','ENSG00000142910',
               'ENSG00000158270','ENSG00000168079','ENSG00000176692','ENSG00000137273','ENSG00000028137','ENSG00000149090',
               'ENSG00000184163','ENSG00000164303','ENSG00000124491','ENSG00000123243')

#use biomart release 105
listEnsembl()
listEnsemblArchives()

ensembl_human <- useEnsembl(biomart = 'genes',host = 'https://dec2021.archive.ensembl.org',version = 105,verbose = TRUE)
listDatasets(mart = ensembl_human)
ensembl_human <- useDataset(dataset = 'hsapiens_gene_ensembl',mart = ensembl_human,verbose = TRUE)

searchAttributes(mart = ensembl_human,pattern = 'symbol')
searchAttributes(mart = ensembl_human,pattern = 'ensembl_gene_id')
searchFilters(mart = ensembl_human,pattern = 'id')
gene_table <- getBM(attributes = c('ensembl_gene_id','hgnc_symbol'),
                    filters = c('ensembl_gene_id'),values = gene_list,
                    mart = ensembl_human,verbose = TRUE)
gene_list <- gene_table$hgnc_symbol


#overlap with our own marker gene
RG_1_marker <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')
intersect(x = gene_list,y = RG_1_marker)

RG_1_marker <- readRDS(file = './res/step_102_fig_230211/RG_1_primate_specific_gene_list.rds')
intersect(x = gene_list,y = RG_1_marker)

RG_1_marker <- readRDS(file = './res/step_102_fig_230211/RG_1_species_conserved_gene_list.rds')
intersect(x = gene_list,y = RG_1_marker)

#the only overlapped gene is CPQ, which is primate conserved.
#GO:0005615 extracellular space
#seems not what we want most

# human specific RG-1 marker ----------------------------------------------
#load ArchR object
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
RG_1_marker <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')

#RG-1 marker related peaks
P2G <- getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = 1,returnLoops = FALSE)
names(P2G@metadata$peakSet) <- paste(P2G@metadata$peakSet@seqnames,as.character(P2G@metadata$peakSet@ranges),sep = '-')
idx <- which(P2G@metadata$geneSet$name %in% RG_1_marker)
idx <- which(P2G$idxRNA %in% idx)
idx <- P2G$idxATAC[idx]
idx <- unique(idx)
peak_set <- P2G@metadata$peakSet[idx]

#overlap with human specific peaks
human_specific_peaks <- readRDS(file = './res/step_73_fig_221102/Brain_ATAC_peak.rds')
human_specific_peaks <- human_specific_peaks$human
temp <- readRDS(file = './res/step_74_fig_221104/RG.1_DAP_list_by_wilcox.rds')
human_specific_peaks <- human_specific_peaks[names(temp)[temp == 'human_specific']]
peak_set <- peak_set[countOverlaps(query = peak_set,subject = human_specific_peaks) > 0]

#overlap with human H3K27ac peaks
human_H3K27ac <- rtracklayer::import.bed(con = '/home/sunym/temp/human.bed')
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToHg19.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp_peak_set <- temp_peak_set[countOverlaps(query = temp_peak_set,subject = human_H3K27ac) > 0]
final_list <- temp_peak_set$ori_peak

#overlap with macaque H3K27ac peaks
macaque_H3K27ac <- rtracklayer::import.bed(con = '/home/sunym/temp/macaque.bed')
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac2.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp_peak_set <- temp_peak_set[countOverlaps(query = temp_peak_set,subject = macaque_H3K27ac) == 0]
final_list <- intersect(x = final_list,y = temp_peak_set$ori_peak)

#overlap with mouse H3K27ac peaks
mouse_H3K27ac <- rtracklayer::import.bed(con = '/home/sunym/temp/mouse_e14.bed')
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToMm10.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp <- data.frame(mouse_id = paste(temp_peak_set@seqnames,as.character(temp_peak_set@ranges),sep = '-'),human_id = temp_peak_set$ori_peak)
rownames(temp) <- temp$mouse_id
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = temp_peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/mm10ToMm9.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp_peak_set$ori_peak <- temp[temp_peak_set$ori_peak,"human_id"]
temp_peak_set <- temp_peak_set[countOverlaps(query = temp_peak_set,subject = mouse_H3K27ac) == 0]
final_list <- intersect(x = final_list,y = temp_peak_set$ori_peak)

#export
#hg38
peak_human <- P2G@metadata$peakSet[final_list]
export.bed(object = peak_human,con = '/home/sunym/trash/hg38_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/hg38_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/hg38/hg38.chrom.sizes',
              '/home/sunym/trash/hg38_peak_list.bb')
system(command = char)

#hg19
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToHg19.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/hg19_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/hg19_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/hg19/hg19.chrom.sizes',
              '/home/sunym/trash/hg19_peak_list.bb')
system(command = char)

#rheMac10
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac10.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/rheMac10_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/rheMac10_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/rheMac10/rheMac10.chrom.sizes',
              '/home/sunym/trash/rheMac10_peak_list.bb')
system(command = char)

#rheMac2
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac2.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/rheMac2_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/rheMac2_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/rheMac2/rheMac2.chrom.sizes',
              '/home/sunym/trash/rheMac2_peak_list.bb')
system(command = char)

#mm10
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToMm10.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/mm10_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/mm10_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/mm10/mm10.chrom.sizes',
              '/home/sunym/trash/mm10_peak_list.bb')
system(command = char)

#mm9
temp <- my_rtracklayer_liftOver(ori_GRanges = temp,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/mm10ToMm9.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/mm9_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/mm9_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/mm9/mm9.chrom.sizes',
              '/home/sunym/trash/mm9_peak_list.bb')
system(command = char)

# RG_1 HGAR that can not be found using bulk data -------------------------

#load ArchR object
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
RG_1_marker <- readRDS(file = './res/step_102_fig_230211/RG_1_human_specific_gene_list.rds')

#RG-1 marker related peaks
P2G <- getPeak2GeneLinks(ArchRProj = Greenleaf_ATAC_ArchR,resolution = 1,returnLoops = FALSE)
names(P2G@metadata$peakSet) <- paste(P2G@metadata$peakSet@seqnames,as.character(P2G@metadata$peakSet@ranges),sep = '-')
idx <- which(P2G@metadata$geneSet$name %in% RG_1_marker)
idx <- which(P2G$idxRNA %in% idx)
idx <- P2G$idxATAC[idx]
idx <- unique(idx)
peak_set <- P2G@metadata$peakSet[idx]

#overlap with human specific peaks
human_specific_peaks <- readRDS(file = './res/step_73_fig_221102/Brain_ATAC_peak.rds')
human_specific_peaks <- human_specific_peaks$human
temp <- readRDS(file = './res/step_74_fig_221104/RG.1_DAP_list_by_wilcox.rds')
human_specific_peaks <- human_specific_peaks[names(temp)[temp == 'human_specific']]
peak_set <- peak_set[countOverlaps(query = peak_set,subject = human_specific_peaks) > 0]

#overlap with human H3K27ac peaks
human_H3K27ac <- rtracklayer::import.bed(con = '/home/sunym/temp/human.bed')
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToHg19.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp_peak_set <- temp_peak_set[countOverlaps(query = temp_peak_set,subject = human_H3K27ac) == 0]
final_list <- temp_peak_set$ori_peak

#overlap with macaque H3K27ac peaks
macaque_H3K27ac <- rtracklayer::import.bed(con = '/home/sunym/temp/macaque.bed')
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac2.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp_peak_set <- temp_peak_set[countOverlaps(query = temp_peak_set,subject = macaque_H3K27ac) == 0]
final_list <- intersect(x = final_list,y = temp_peak_set$ori_peak)

#overlap with mouse H3K27ac peaks
mouse_H3K27ac <- rtracklayer::import.bed(con = '/home/sunym/temp/mouse_e14.bed')
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToMm10.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp <- data.frame(mouse_id = paste(temp_peak_set@seqnames,as.character(temp_peak_set@ranges),sep = '-'),human_id = temp_peak_set$ori_peak)
rownames(temp) <- temp$mouse_id
temp_peak_set <- my_rtracklayer_liftOver(ori_GRanges = temp_peak_set,chain_file = './data/reference/UCSC_chain_file_for_liftOver/mm10ToMm9.over.chain',
                                         merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
temp_peak_set$ori_peak <- temp[temp_peak_set$ori_peak,"human_id"]
temp_peak_set <- temp_peak_set[countOverlaps(query = temp_peak_set,subject = mouse_H3K27ac) == 0]
final_list <- intersect(x = final_list,y = temp_peak_set$ori_peak)

#export
#hg38
peak_human <- P2G@metadata$peakSet[final_list]
export.bed(object = peak_human,con = '/home/sunym/trash/hg38_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/hg38_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/hg38/hg38.chrom.sizes',
              '/home/sunym/trash/hg38_peak_list.bb')
system(command = char)

#hg19
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToHg19.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/hg19_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/hg19_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/hg19/hg19.chrom.sizes',
              '/home/sunym/trash/hg19_peak_list.bb')
system(command = char)

#rheMac10
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac10.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/rheMac10_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/rheMac10_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/rheMac10/rheMac10.chrom.sizes',
              '/home/sunym/trash/rheMac10_peak_list.bb')
system(command = char)

#rheMac2
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToRheMac2.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/rheMac2_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/rheMac2_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/rheMac2/rheMac2.chrom.sizes',
              '/home/sunym/trash/rheMac2_peak_list.bb')
system(command = char)

#mm10
temp <- my_rtracklayer_liftOver(ori_GRanges = peak_human,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/hg38ToMm10.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/mm10_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/mm10_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/mm10/mm10.chrom.sizes',
              '/home/sunym/trash/mm10_peak_list.bb')
system(command = char)

#mm9
temp <- my_rtracklayer_liftOver(ori_GRanges = temp,
                                chain_file = './data/reference/UCSC_chain_file_for_liftOver/mm10ToMm9.over.chain',
                                merge = TRUE,workers = 6,future.globals.maxSize = 200*(1024^3))
export.bed(object = temp,con = '/home/sunym/trash/mm9_peak_list.bed')
char <- paste('/content/data/sunym/script/bed2bigBed','/home/sunym/trash/mm9_peak_list.bed',
              '/content/data/sunym/project/Brain/data/reference/UCSC_chrom_size/mm9/mm9.chrom.sizes',
              '/home/sunym/trash/mm9_peak_list.bb')
system(command = char)
