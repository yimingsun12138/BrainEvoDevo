#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: validate my modified wilcoxon test using macaque multiome data  ##
## Data: 2022.06.19                                                                ##
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
library(networkD3)
library(htmlwidgets)
library(DESeq2)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# load data ---------------------------------------------------------------
macaque_multiome_ArchR <- loadArchRProject(path = './ArchR/processed_data/macaque_multiome_ArchR_220411/')
macaque_peak_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix',verbose = TRUE)

meta_data <- macaque_peak_matrix@colData
ii <- macaque_peak_matrix@rowRanges
ii <- paste(ii@seqnames,as.character(ii@ranges),sep = '-')
macaque_peak_matrix <- macaque_peak_matrix@assays@data$PeakMatrix
rownames(macaque_peak_matrix) <- ii
macaque_peak_matrix <- CreateSeuratObject(counts = macaque_peak_matrix,project = 'ATAC',assay = 'RNA',meta.data = as.data.frame(meta_data),min.cells = 0,min.features = 0)
gc()

# wilcoxon test compare my method and ArchR default------------------------------------------

## RG vs Ex-1 --------------------------------------------------------------
cell_type_1 <- 'RG'
cell_type_2 <- 'Ex-1'


### default NF same cell ----------------------------------------------------
#ArchR method
marker_ArchR <- getMarkerFeatures(ArchRProj = macaque_multiome_ArchR,
                                  groupBy = 'cell_type',
                                  useGroups = cell_type_1,
                                  bgdGroups = cell_type_2,
                                  useMatrix = 'PeakMatrix',
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  normBy = NULL,
                                  testMethod = 'wilcoxon',
                                  maxCells = 1000,
                                  scaleTo = 10^4,
                                  k = 100,
                                  bufferRatio = 0.8,
                                  binarize = FALSE,
                                  verbose = TRUE)

ii <- marker_ArchR@elementMetadata
ii <- paste(ii$seqnames,ii$start,ii$end,sep = '-')
marker_ArchR <- data.frame(log2FC = marker_ArchR@assays@data$Log2FC$x,
                           fdr = marker_ArchR@assays@data$FDR$x,
                           pval = marker_ArchR@assays@data$Pval$x)
rownames(marker_ArchR) <- ii

#my method
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8)

mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

#normalize
temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
normfactor <- (normfactor * 10^4)/median(normfactor[colnames(macaque_peak_matrix@assays$RNA@counts)] * colSums(macaque_peak_matrix@assays$RNA@counts))
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
normfactor <- (normfactor * 10^4)/median(normfactor[colnames(macaque_peak_matrix@assays$RNA@counts)] * colSums(macaque_peak_matrix@assays$RNA@counts))
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_2@assays$RNA@counts)
colnames(temp) <- colnames(mat_2@assays$RNA@counts)
mat_2 <- temp

table(rownames(mat_1) == rownames(mat_2))
temp <- rowSums(mat_1) + rowSums(mat_2)
gc()

#wilcoxon
marker_my <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare my and ArchR wilcoxon test
table(rownames(marker_ArchR) %in% rownames(marker_my))
ii <- dplyr::intersect(rownames(marker_ArchR),rownames(marker_my))
marker_ArchR <- marker_ArchR[ii,]
marker_my <- marker_my[ii,]

pdf(file = './res/step_38_fig_220619/Log2FC_default_NF_same_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$log2FC,my = marker_my$log2FC),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'Log2FC default NF same cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()

pdf(file = './res/step_38_fig_220619/fdr_default_NF_same_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$fdr,my = marker_my$fdr),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'fdr default NF same cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()


### default DF different cells ----------------------------------------------
marker_ArchR <- getMarkerFeatures(ArchRProj = macaque_multiome_ArchR,
                                  groupBy = 'cell_type',
                                  useGroups = cell_type_1,
                                  bgdGroups = cell_type_2,
                                  useMatrix = 'PeakMatrix',
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  normBy = NULL,
                                  testMethod = 'wilcoxon',
                                  maxCells = 1000,
                                  scaleTo = 10^4,
                                  k = 100,
                                  bufferRatio = 0.8,
                                  binarize = FALSE,
                                  verbose = TRUE)

ii <- marker_ArchR@elementMetadata
ii <- paste(ii$seqnames,ii$start,ii$end,sep = '-')
marker_ArchR <- data.frame(log2FC = marker_ArchR@assays@data$Log2FC$x,
                           fdr = marker_ArchR@assays@data$FDR$x,
                           pval = marker_ArchR@assays@data$Pval$x)
rownames(marker_ArchR) <- ii

match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 34)

mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

#normalize
temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
normfactor <- (normfactor * 10^4)/median(normfactor[colnames(macaque_peak_matrix@assays$RNA@counts)] * colSums(macaque_peak_matrix@assays$RNA@counts))
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
normfactor <- (normfactor * 10^4)/median(normfactor[colnames(macaque_peak_matrix@assays$RNA@counts)] * colSums(macaque_peak_matrix@assays$RNA@counts))
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_2@assays$RNA@counts)
colnames(temp) <- colnames(mat_2@assays$RNA@counts)
mat_2 <- temp

table(rownames(mat_1) == rownames(mat_2))
temp <- rowSums(mat_1) + rowSums(mat_2)
gc()

#wilcoxon
marker_my <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare my and ArchR wilcoxon test
table(rownames(marker_ArchR) %in% rownames(marker_my))
ii <- dplyr::intersect(rownames(marker_ArchR),rownames(marker_my))
marker_ArchR <- marker_ArchR[ii,]
marker_my <- marker_my[ii,]

pdf(file = './res/step_38_fig_220619/Log2FC_default_NF_different_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$log2FC,my = marker_my$log2FC),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'Log2FC default NF different cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()

pdf(file = './res/step_38_fig_220619/pval_default_NF_different_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$pval,my = marker_my$pval),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'pval default NF different cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()

#define marker
#log2fc > 1,fdr < 0.01
marker_ArchR <- rownames(marker_ArchR)[marker_ArchR$log2FC > 1 & marker_ArchR$fdr < 0.01]
marker_my <- rownames(marker_my)[marker_my$log2FC > 1 & marker_my$fdr < 0.01]
temp <- list(marker_ArchR = marker_ArchR,marker_my = marker_my)

pdf(file = './res/step_38_fig_220619/RG_vs_Ex_1_marker_different_cells_vennplot.pdf',width = 5,height = 5)
ggvenn(data = temp,c('marker_ArchR','marker_my'))
dev.off()


### different NF same cells -------------------------------------------------
marker_ArchR <- getMarkerFeatures(ArchRProj = macaque_multiome_ArchR,
                                  groupBy = 'cell_type',
                                  useGroups = cell_type_1,
                                  bgdGroups = cell_type_2,
                                  useMatrix = 'PeakMatrix',
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  normBy = NULL,
                                  testMethod = 'wilcoxon',
                                  maxCells = 1000,
                                  scaleTo = 10^4,
                                  k = 100,
                                  bufferRatio = 0.8,
                                  binarize = FALSE,
                                  verbose = TRUE)

ii <- marker_ArchR@elementMetadata
ii <- paste(ii$seqnames,ii$start,ii$end,sep = '-')
marker_ArchR <- data.frame(log2FC = marker_ArchR@assays@data$Log2FC$x,
                           fdr = marker_ArchR@assays@data$FDR$x,
                           pval = marker_ArchR@assays@data$Pval$x)
rownames(marker_ArchR) <- ii

match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8)

mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

#normalize
temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_2@assays$RNA@counts)
colnames(temp) <- colnames(mat_2@assays$RNA@counts)
mat_2 <- temp

table(rownames(mat_1) == rownames(mat_2))
temp <- rowSums(mat_1) + rowSums(mat_2)
gc()

#wilcoxon
marker_my <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare my and ArchR wilcoxon test
table(rownames(marker_ArchR) %in% rownames(marker_my))
ii <- dplyr::intersect(rownames(marker_ArchR),rownames(marker_my))
marker_ArchR <- marker_ArchR[ii,]
marker_my <- marker_my[ii,]

pdf(file = './res/step_38_fig_220619/Log2FC_different_NF_same_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$log2FC,my = marker_my$log2FC),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'Log2FC different NF same cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()

pdf(file = './res/step_38_fig_220619/pval_different_NF_same_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$pval,my = marker_my$pval),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'pval different NF same cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()

#define marker
#log2fc > 1,fdr < 0.01
marker_ArchR <- rownames(marker_ArchR)[marker_ArchR$log2FC > 1 & marker_ArchR$fdr < 0.01]
marker_my <- rownames(marker_my)[marker_my$log2FC > 1 & marker_my$fdr < 0.01]
temp <- list(marker_ArchR = marker_ArchR,marker_my = marker_my)

pdf(file = './res/step_38_fig_220619/RG_vs_Ex_1_marker_different_NF_vennplot.pdf',width = 5,height = 5)
ggvenn(data = temp,c('marker_ArchR','marker_my'))
dev.off()

### different NF all cells -------------------------------------------------
marker_ArchR <- getMarkerFeatures(ArchRProj = macaque_multiome_ArchR,
                                  groupBy = 'cell_type',
                                  useGroups = cell_type_1,
                                  bgdGroups = cell_type_2,
                                  useMatrix = 'PeakMatrix',
                                  bias = c("TSSEnrichment", "log10(nFrags)"),
                                  normBy = NULL,
                                  testMethod = 'wilcoxon',
                                  maxCells = 1000,
                                  scaleTo = 10^4,
                                  k = 100,
                                  bufferRatio = 0.8,
                                  binarize = FALSE,
                                  verbose = TRUE)

ii <- marker_ArchR@elementMetadata
ii <- paste(ii$seqnames,ii$start,ii$end,sep = '-')
marker_ArchR <- data.frame(log2FC = marker_ArchR@assays@data$Log2FC$x,
                           fdr = marker_ArchR@assays@data$FDR$x,
                           pval = marker_ArchR@assays@data$Pval$x)
rownames(marker_ArchR) <- ii

mat_1 <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type_1]
mat_2 <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type_2]

#normalize
temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_2@assays$RNA@counts)
colnames(temp) <- colnames(mat_2@assays$RNA@counts)
mat_2 <- temp

table(rownames(mat_1) == rownames(mat_2))
temp <- rowSums(mat_1) + rowSums(mat_2)
gc()

#wilcoxon
marker_my <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare my and ArchR wilcoxon test
table(rownames(marker_ArchR) %in% rownames(marker_my))
ii <- dplyr::intersect(rownames(marker_ArchR),rownames(marker_my))
marker_ArchR <- marker_ArchR[ii,]
marker_my <- marker_my[ii,]

pdf(file = './res/step_38_fig_220619/Log2FC_different_NF_all_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$log2FC,my = marker_my$log2FC),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'Log2FC different NF all cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()

pdf(file = './res/step_38_fig_220619/pval_different_NF_all_cells_RG_vs_Ex_1.pdf',width = 5,height = 5)
ggplot(data = data.frame(ArchR = marker_ArchR$pval,my = marker_my$pval),aes(x = ArchR,y = my)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'pval different NF all cells') + 
  xlab('ArchR method') + ylab('my method')
dev.off()

#define marker
#log2fc > 1,fdr < 0.01
marker_ArchR <- rownames(marker_ArchR)[marker_ArchR$log2FC > 1 & marker_ArchR$fdr < 0.01]
marker_my <- rownames(marker_my)[marker_my$log2FC > 1 & marker_my$fdr < 0.01]
temp <- list(marker_ArchR = marker_ArchR,marker_my = marker_my)

pdf(file = './res/step_38_fig_220619/RG_vs_Ex_1_marker_different_NF_all_cells_vennplot.pdf',width = 5,height = 5)
ggvenn(data = temp,c('marker_ArchR','marker_my'))
dev.off()