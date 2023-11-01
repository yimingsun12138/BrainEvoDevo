#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: validate species DAP method                                     ##
## Data: 2022.05.27                                                                ##
## Author: Yiming Sun                                                              ##
#####################################################################################

# description:
# continue the step 30 script to validate the parameter while doing liftover
# in last script, I validate the 1.macaque to human 2.mouse to human 3.human to mouse liftover parameter
# in this script, I am going to validate the human to macaque liftover parameter

#sleep
ii <- 1
while(1){
  cat(paste("round",ii),sep = "\n")
  ii <- ii+1
  Sys.sleep(30)
}

#general setting
setwd('/data/User/sunym/project/Brain/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(parallel)
library(Seurat)
library(ArchR)
library(rtracklayer)
library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(ComplexHeatmap)
library(parallel)
library(patchwork)
library(ggpubr)
library(ggvenn)
library(circlize)
library(ChIPpeakAnno)
library(GenomicFeatures)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(BSgenome.Mmusculus.UCSC.mm10)
library(EnrichedHeatmap)
library(circlize)
library(scales)
library(topGO)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggrepel)
library(ggpointdensity)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# validate liftover mismatch parameter ------------------------------------
#load data
Greenleaf_ATAC_ArchR <- readRDS(file = './data/public/Chromatin_and_gene_regulatory_dynamics_of_the_developing_human_cerebral_cortex_at_single_cell_resolution/scATAC_seq/ArchR/processed_data/Greenleaf_ATAC_ArchR_220412/Save-ArchR-Project.rds')
macaque_multiome_ArchR <- readRDS(file = './ArchR/processed_data/macaque_multiome_ArchR_220411/Save-ArchR-Project.rds')
mouse_ATAC_ArchR <- readRDS(file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/ArchR/processed_data/mouse_ATAC_ArchR_220414/Save-ArchR-Project.rds')

human_peakset <- getPeakSet(ArchRProj = Greenleaf_ATAC_ArchR)
macaque_peakset <- getPeakSet(ArchRProj = macaque_multiome_ArchR)
mouse_peakset <- getPeakSet(ArchRProj = mouse_ATAC_ArchR)

## create mouse_human_merged_peakset ---------------------------------------
#mouse liftover to human, mismatch = 0.4, delta_length_ratio = 0.4
lifted_mouse_peakset <- my_unique_peakset_liftover(ori_GRanges = mouse_peakset,
                                                   UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                   chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/mm10ToHg38.over.chain',
                                                   liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                   chr_filter = TRUE,mapped_chr = unique(as.character(human_peakset@seqnames)),
                                                   overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220527')

mouse_human_merged_peakset <- my_bedtools_merge(peakset_x = lifted_mouse_peakset$mapped,peakset_y = human_peakset,
                                                bedtools_path = '/data/User/sunym/software/bedtools/bedtools',
                                                d = 0,bedtools_param = NULL,tmp_path = '/data/User/sunym/temp/tmp_220527')

mouse_human_merged_peakset <- my_unique_peakset_liftover(ori_GRanges = mouse_human_merged_peakset,
                                                         UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                         chain_file = './data/public/Multimodal_profiling_of_the_transcriptional_regulatory_landscape_of_the_developing_mouse_cortex_identifies_Neurog2_as_a_key_epigenome_remodeler/scATAC_seq/reference/hg38ToMm10.over.chain',
                                                         liftOver_mismatch = 0.4,length_filter = TRUE,length_mismatch = 0.4,
                                                         chr_filter = TRUE,mapped_chr = unique(as.character(mouse_peakset@seqnames)),
                                                         overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220527')

# saveRDS(mouse_human_merged_peakset,file = './res/step_30_fig_220527/mouse_human_merged_peakset.rds')

## mouse_human_merged_peakset liftover to macaque with different mismatch--------------------------
lifted_human_peakset <- list()
for(i in c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9)){
  char <- paste('mismatch',as.character(i),sep = '_')
  lifted_human_peakset[[char]] <- my_unique_peakset_liftover(ori_GRanges = mouse_human_merged_peakset$ori,
                                                             UCSC_liftOver_path = '/data/User/sunym/software/UCSC_LiftOver/liftOver',
                                                             chain_file = './data/reference/hg38ToRheMac10.over.chain',
                                                             liftOver_mismatch = i,length_filter = FALSE,length_mismatch = 0.1,
                                                             chr_filter = TRUE,mapped_chr = unique(as.character(macaque_peakset@seqnames)),
                                                             overlap_filter = TRUE,tmp_path = '/data/User/sunym/temp/tmp_220527')
}

my_send_sms('liftover done!')

# #save data
# saveRDS(lifted_human_peakset,file = './res/step_30_fig_220527/lifted_mouse_human_merged_peakset_to_macaque_with_different_mismatch.rds')

#load data
mouse_human_merged_peakset <- readRDS(file = './res/step_30_fig_220527/mouse_human_merged_peakset.rds')
lifted_human_peakset <- readRDS(file = './res/step_30_fig_220527/lifted_mouse_human_merged_peakset_to_macaque_with_different_mismatch.rds')

temp <- base::lapply(X = c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9),FUN = function(x){
  char <- paste('mismatch',as.character(x),sep = '_')
  ratio <- length(lifted_human_peakset[[char]]$mapped)/length(mouse_human_merged_peakset$ori)
  return(ratio)
})
temp <- unlist(temp)

temp <- data.frame(mismatch = c(0.01,0.05,0.1,0.2,0.3,0.4,0.5,0.7,0.9),ratio = temp)

pdf(file = './res/step_30_fig_220527/lifted_mouse_human_merged_peakset_to_macaque_mismatch_preserved_ratio_dotplot.pdf',width = 4,height = 3)
ggplot(data = temp,aes(x=mismatch,y=ratio)) + 
  geom_point(size = 2,color = 'black') + 
  geom_line(color = 'grey',linetype = 'dashed') + 
  ylim(c(0.7,1)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.6,
        axis.title = element_text(face = 'bold',size = 14)) + 
  xlab('mismatch') + ylab('preserved ratio')
dev.off()


## sequence pairwise alignment ---------------------------------------------
my_pairwise_alighment <- function(peak_query,peak_subject,query_genome,subject_genome){
  #create NSM
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
  subject_length <- peak_subject@ranges@width
  
  #get sequence
  peak_query <- getSeq(x = query_genome,names = peak_query)
  peak_subject <- getSeq(x = subject_genome,names = peak_subject)
  
  #alighment score
  alighment_score <- c()
  #normal alighment
  temp <- pairwiseAlignment(pattern = peak_query,subject = peak_subject,type="global",
                            substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  alighment_score <- append(alighment_score,temp)
  #reverse alighment
  temp <- pairwiseAlignment(pattern = reverse(peak_query),subject = peak_subject,type="global",
                            substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  alighment_score <- append(alighment_score,temp)
  #reverse complement alighment
  temp <- pairwiseAlignment(pattern = reverseComplement(peak_query),subject = peak_subject,type="global",
                            substitutionMatrix = NSM,gapOpening=10, gapExtension=4,scoreOnly = TRUE)
  alighment_score <- append(alighment_score,temp)
  
  #return
  alighment_score <- max(alighment_score)/subject_length
  return(alighment_score)
}

#mismatch_0.01
alignment_0.01 <- base::lapply(X = 1:100,FUN = function(x){
  temp <- my_pairwise_alighment(peak_query = lifted_human_peakset$mismatch_0.01$mapped[x],
                                peak_subject = lifted_human_peakset$mismatch_0.01$ori[x],
                                query_genome = BSgenome.Mmulatta.UCSC.rheMac10,
                                subject_genome = BSgenome.Hsapiens.UCSC.hg38)
  return(temp)
})

#the function runs too slow, try to do it on cluster!

#load data
lifted_human_peakset <- readRDS(file = './res/step_30_fig_220527/lifted_mouse_human_merged_peakset_to_macaque_with_different_mismatch.rds')
mouse_human_merged_peakset <- readRDS(file = './res/step_30_fig_220527/mouse_human_merged_peakset.rds')
alignment_score <- readRDS('./res/step_30_fig_220527/lifted_mouse_human_merged_peak_alignment_score_matrix.rds')

#alignment score distribution
pdf(file = './res/step_30_fig_220527/lifted_mouse_human_merged_peak_alignment_score_distribution.pdf',width = 6,height = 4)
ggplot(data = alignment_score,aes(x=exp_align,color=group)) + 
  geom_density() + 
  geom_vline(xintercept = 0.25,color = 'grey',linetype = 'dashed') + 
  theme_cowplot() + 
  scale_color_manual(values = as.character(ArchRPalettes$circus)) + 
  theme(aspect.ratio = 0.65,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'alignment score distribution') + 
  xlab('exp(alignment score)') + ylab('density')
dev.off()

#let the liftover mismatch = 0.1

#delta length and alignment score
alignment_score <- alignment_score[alignment_score$group == 'mismatch_0.1',]
temp <- alignment_score$peak
temp_ori <- lifted_human_peakset$mismatch_0.1$ori
temp_mapped <- lifted_human_peakset$mismatch_0.1$mapped
names(temp_mapped) <- temp_mapped$name
delta_length_ratio <- (temp_mapped[temp]@ranges@width - temp_ori[temp]@ranges@width)/temp_ori[temp]@ranges@width
temp <- data.frame(exp_align = alignment_score$exp_align,length_ratio = delta_length_ratio)

pdf(file = './res/step_30_fig_220527/lifted_mouse_human_merged_peak_alignment_score_vs_delta_length_ratio_dotplot.pdf',width = 6,height = 4)
ggplot(data = temp,aes(x=exp_align,y=length_ratio)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  ylim(c(-1,1)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.65,
        axis.title = element_text(size = 12,face = 'bold')) + 
  geom_hline(yintercept = 0.1,color = 'red',linetype = 'dashed') + 
  geom_hline(yintercept = -0.1,color = 'red',linetype = 'dashed') + 
  xlab('exp(alignment score)') + ylab('delta length ratio')
dev.off()


#let the liftover mismatch = 0.4
alignment_score <- readRDS('./res/step_30_fig_220527/lifted_mouse_human_merged_peak_alignment_score_matrix.rds')

#delta length and alignment score
alignment_score <- alignment_score[alignment_score$group == 'mismatch_0.4',]
temp <- alignment_score$peak
temp_ori <- lifted_human_peakset$mismatch_0.4$ori
temp_mapped <- lifted_human_peakset$mismatch_0.4$mapped
names(temp_mapped) <- temp_mapped$name
delta_length_ratio <- (temp_mapped[temp]@ranges@width - temp_ori[temp]@ranges@width)/temp_ori[temp]@ranges@width
temp <- data.frame(exp_align = alignment_score$exp_align,length_ratio = delta_length_ratio)

pdf(file = './res/step_30_fig_220527/lifted_mouse_human_merged_peak_alignment_score_vs_delta_length_ratio_dotplot_with_mismatch_0.4.pdf',width = 6,height = 4)
ggplot(data = temp,aes(x=exp_align,y=length_ratio)) + 
  geom_pointdensity(size = 0.1) + 
  scale_color_viridis() + 
  ylim(c(-1,1)) + 
  theme_cowplot() + 
  theme(aspect.ratio = 0.65,
        axis.title = element_text(size = 12,face = 'bold')) + 
  geom_hline(yintercept = 0.4,color = 'red',linetype = 'dashed') + 
  geom_hline(yintercept = -0.4,color = 'red',linetype = 'dashed') + 
  geom_vline(xintercept = 0.25,color = 'grey',linetype = 'dashed') + 
  xlab('exp(alignment score)') + ylab('delta length ratio')
dev.off()
