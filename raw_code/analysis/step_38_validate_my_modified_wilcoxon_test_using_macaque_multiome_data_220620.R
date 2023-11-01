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

temp <- readRDS(file = './res/step_35_fig_220606/macaque_peak_matrix_Seurat.rds')
macaque_peak_matrix$Reads_In_Cell_Type_Peaks <- temp@meta.data[colnames(macaque_peak_matrix@assays$RNA@counts),"Reads_In_Cell_Type_Peaks"]
macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks <- macaque_peak_matrix@meta.data[rownames(macaque_multiome_ArchR@cellColData),"Reads_In_Cell_Type_Peaks"]
rm(temp)
gc()

# validate parameter while doing wilcoxon ---------------------------------

## RG vs Ex-1 --------------------------------------------------------------
cell_type_1 <- 'RG'
cell_type_2 <- 'Ex-1'

### validate NF (median/i) --------------------------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#ArchR default normalize
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

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

marker_1 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

# median/i normalize
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

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

marker_2 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare
table(rownames(marker_1) %in% rownames(marker_2))
ii <- dplyr::intersect(rownames(marker_1),rownames(marker_2))
marker_1 <- marker_1[ii,]
marker_2 <- marker_2[ii,]

char <- paste(cell_type_1,'vs',cell_type_2,'Log2FC my NF same cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$log2FC,m2 = marker_2$log2FC),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('default NF') + ylab('median/i NF')
dev.off()

char <- paste(cell_type_1,'vs',cell_type_2,'pval my NF same cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$pval,m2 = marker_2$pval),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('default NF') + ylab('median/i NF')
dev.off()

#log2fc > 1,fdr < 0.01
marker_1 <- rownames(marker_1)[marker_1$log2FC > 1 & marker_1$fdr < 0.01]
marker_2 <- rownames(marker_2)[marker_2$log2FC > 1 & marker_2$fdr < 0.01]
temp <- list(marker_1 = marker_1,marker_2 = marker_2)

char <- paste(cell_type_1,'vs',cell_type_2,'marker overlap my NF same cells',sep = ' ')
char <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
char <- gsub(pattern = '-',replacement = '_',x = char,fixed = TRUE)
char <- paste0('./res/step_38_fig_220620/',char,'.pdf')

pdf(file = char,width = 5,height = 5)
ggvenn(data = temp,c('marker_1','marker_2'))
dev.off()

### validate NF (very small) --------------------------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#ArchR default normalize
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

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

marker_1 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

# very small normalize
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- (median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS)/(10^4)
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- (median(macaque_multiome_ArchR$ReadsInTSS)/macaque_multiome_ArchR$ReadsInTSS)/(10^4)
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

marker_2 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare
table(rownames(marker_1) %in% rownames(marker_2))
ii <- dplyr::intersect(rownames(marker_1),rownames(marker_2))
marker_1 <- marker_1[ii,]
marker_2 <- marker_2[ii,]

char <- paste(cell_type_1,'vs',cell_type_2,'Log2FC small NF same cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$log2FC,m2 = marker_2$log2FC),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('default NF') + ylab('small NF')
dev.off()

char <- paste(cell_type_1,'vs',cell_type_2,'pval small NF same cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$pval,m2 = marker_2$pval),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('default NF') + ylab('small NF')
dev.off()

#log2fc > 1,fdr < 0.01
marker_1 <- rownames(marker_1)[marker_1$log2FC > 1 & marker_1$fdr < 0.01]
marker_2 <- rownames(marker_2)[marker_2$log2FC > 1 & marker_2$fdr < 0.01]
temp <- list(marker_1 = marker_1,marker_2 = marker_2)

char <- paste(cell_type_1,'vs',cell_type_2,'marker overlap small NF same cells',sep = ' ')
char <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
char <- gsub(pattern = '-',replacement = '_',x = char,fixed = TRUE)
char <- paste0('./res/step_38_fig_220620/',char,'.pdf')

pdf(file = char,width = 5,height = 5)
ggvenn(data = temp,c('marker_1','marker_2'))
dev.off()

### validate cell sampling --------------------------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#my NF
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

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

marker_1 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#redo sampling
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 34)

#my NF
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

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

marker_2 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare
table(rownames(marker_1) %in% rownames(marker_2))
ii <- dplyr::intersect(rownames(marker_1),rownames(marker_2))
marker_1 <- marker_1[ii,]
marker_2 <- marker_2[ii,]

char <- paste(cell_type_1,'vs',cell_type_2,'Log2FC fifferent cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$log2FC,m2 = marker_2$log2FC),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('seed 1') + ylab('seed 34')
dev.off()

char <- paste(cell_type_1,'vs',cell_type_2,'pval different cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$pval,m2 = marker_2$pval),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('seed 1') + ylab('seed 34')
dev.off()

#log2fc > 1,fdr < 0.01
marker_1 <- rownames(marker_1)[marker_1$log2FC > 1 & marker_1$fdr < 0.01]
marker_2 <- rownames(marker_2)[marker_2$log2FC > 1 & marker_2$fdr < 0.01]
temp <- list(marker_1 = marker_1,marker_2 = marker_2)

char <- paste(cell_type_1,'vs',cell_type_2,'marker overlap different cells',sep = ' ')
char <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
char <- gsub(pattern = '-',replacement = '_',x = char,fixed = TRUE)
char <- paste0('./res/step_38_fig_220620/',char,'.pdf')

pdf(file = char,width = 5,height = 5)
ggvenn(data = temp,c('marker_1','marker_2'))
dev.off()

### validate all cell --------------------------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#my NF
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

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

marker_1 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#my NF
mat_1 <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type_1]
mat_2 <- macaque_peak_matrix[,macaque_peak_matrix$cell_type == cell_type_2]

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

marker_2 <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#compare
table(rownames(marker_1) %in% rownames(marker_2))
ii <- dplyr::intersect(rownames(marker_1),rownames(marker_2))
marker_1 <- marker_1[ii,]
marker_2 <- marker_2[ii,]

char <- paste(cell_type_1,'vs',cell_type_2,'Log2FC all cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$log2FC,m2 = marker_2$log2FC),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('sampled cells') + ylab('all cells')
dev.off()

char <- paste(cell_type_1,'vs',cell_type_2,'pval all cells',sep = ' ')
temp <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
temp <- gsub(pattern = '-',replacement = '_',x = temp,fixed = TRUE)
temp <- paste0('./res/step_38_fig_220620/',temp,'.pdf')

pdf(file = temp,width = 5,height = 5)
ggplot(data = data.frame(m1 = marker_1$pval,m2 = marker_2$pval),aes(x = m1,y = m2)) + 
  geom_pointdensity(size = 0.1) + 
  theme_cowplot() + 
  scale_color_viridis() + 
  geom_smooth(method = 'lm',size = 0.5,color = 'blue') + 
  geom_abline(slope = 1,intercept = 0,size = 0.5,color = 'red',linetype = 'dashed') + 
  stat_cor(method = 'pearson') + 
  theme(aspect.ratio = 1,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = char) + 
  xlab('sampled cells') + ylab('all cells')
dev.off()

#log2fc > 1,fdr < 0.01
marker_1_backup <- marker_1
marker_2_backup <- marker_2
marker_1 <- rownames(marker_1)[marker_1$log2FC > 1 & marker_1$fdr < 0.01]
marker_2 <- rownames(marker_2)[marker_2$log2FC > 1 & marker_2$fdr < 0.01]
temp <- list(marker_1 = marker_1,marker_2 = marker_2)

char <- paste(cell_type_1,'vs',cell_type_2,'marker overlap all cells',sep = ' ')
char <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
char <- gsub(pattern = '-',replacement = '_',x = char,fixed = TRUE)
char <- paste0('./res/step_38_fig_220620/',char,'.pdf')

pdf(file = char,width = 5,height = 5)
ggvenn(data = temp,c('marker_1','marker_2'))
dev.off()

#why much more marker enriched using all cells?
temp <- marker_2[!(marker_2 %in% marker_1)]
temp <- marker_1_backup[temp,]

pdf(file = './res/step_38_fig_220620/RG_vs_Ex_1_all_cell_specific_marker_fail_to_detect_while_sampling_log2fc_densityplot.pdf',width = 6,height = 4)
ggplot(data = temp) + 
  geom_density(aes(x = log2FC),color = 'black',fill = 'grey') + 
  theme_cowplot() + 
  geom_vline(xintercept = 1,color = 'red',size = 0.5,linetype = 'dashed') + 
  theme(aspect.ratio = 0.6,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'all cell specific marker fail to detect while sampling')
dev.off()

temp$logfdr <- -log10(temp$fdr)
pdf(file = './res/step_38_fig_220620/RG_vs_Ex_1_all_cell_specific_marker_fail_to_detect_while_sampling_fdr_densityplot.pdf',width = 6,height = 4)
ggplot(data = temp) + 
  geom_density(aes(x = logfdr),color = 'black',fill = 'grey') + 
  theme_cowplot() + 
  geom_vline(xintercept = 2,color = 'red',size = 0.5,linetype = 'dashed') + 
  theme(aspect.ratio = 0.6,
        plot.title = element_text(size = 16,face = 'bold',hjust = 0.5),
        axis.title = element_text(size = 12,face = 'bold')) + 
  labs(title = 'all cell specific marker fail to detect while sampling') + 
  xlab('-log10(fdr)')
dev.off()

### compare NF factor using RIT RIP nFrags ----------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#ReadsInTSS
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

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

marker_ReadsInTSS <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#ReadsInPeaks
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
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

marker_ReadsInPeaks <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#nFrags
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd$RG$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
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

marker_nFrags <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#venn plot
marker_ReadsInTSS <- rownames(marker_ReadsInTSS)[marker_ReadsInTSS$log2FC > 1 & marker_ReadsInTSS$fdr < 0.01]
marker_ReadsInPeaks <- rownames(marker_ReadsInPeaks)[marker_ReadsInPeaks$log2FC > 1 & marker_ReadsInPeaks$fdr < 0.01]
marker_nFrags <- rownames(marker_nFrags)[marker_nFrags$log2FC > 1 & marker_nFrags$fdr < 0.01]

temp <- list(marker_ReadsInTSS = marker_ReadsInTSS,
             marker_ReadsInPeaks = marker_ReadsInPeaks,
             marker_nFrags = marker_nFrags)

pdf(file = './res/step_38_fig_220620/RG_vs_Ex_1_marker_vennplot_using_RIT_RIP_nFrags_as_NF.pdf',width = 7,height = 7)
ggvenn(data = temp,c('marker_ReadsInTSS','marker_ReadsInPeaks','marker_nFrags'))
dev.off()

## Ex-1 vs Ex-2 --------------------------------------------------------------
cell_type_1 <- 'Ex-1'
cell_type_2 <- 'Ex-2'

### compare NF factor using RIT RIP nFrags ----------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#ReadsInTSS
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

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

marker_ReadsInTSS <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#ReadsInPeaks
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
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

marker_ReadsInPeaks <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#nFrags
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
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

marker_nFrags <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#Reads_In_Cell_Type_Peaks
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks)/macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks)/macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
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

marker_Reads_In_Cell_Type_Peaks <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#venn plot
marker_ReadsInTSS <- rownames(marker_ReadsInTSS)[marker_ReadsInTSS$log2FC > 1 & marker_ReadsInTSS$fdr < 0.01]
marker_ReadsInPeaks <- rownames(marker_ReadsInPeaks)[marker_ReadsInPeaks$log2FC > 1 & marker_ReadsInPeaks$fdr < 0.01]
marker_nFrags <- rownames(marker_nFrags)[marker_nFrags$log2FC > 1 & marker_nFrags$fdr < 0.01]
marker_Reads_In_Cell_Type_Peaks <- rownames(marker_Reads_In_Cell_Type_Peaks)[marker_Reads_In_Cell_Type_Peaks$log2FC > 1 & marker_Reads_In_Cell_Type_Peaks$fdr < 0.01]

temp <- list(marker_ReadsInTSS = marker_ReadsInTSS,
             marker_ReadsInPeaks = marker_ReadsInPeaks,
             marker_nFrags = marker_nFrags,
             marker_Reads_In_Cell_Type_Peaks = marker_Reads_In_Cell_Type_Peaks)

if(is.null(cell_type_2)){
  cell_type_2 <- 'others'
}
char <- paste(cell_type_1,'vs',cell_type_2,'marker_vennplot_using_RIT_RIP_nFrags_as_NF',sep = ' ')
char <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
char <- gsub(pattern = '-',replacement = '_',x = char,fixed = TRUE)
char <- paste0('./res/step_38_fig_220620/',char,'.pdf')

pdf(file = char,width = 10,height = 10)
ggvenn(data = temp,c('marker_ReadsInTSS','marker_ReadsInPeaks','marker_nFrags','marker_Reads_In_Cell_Type_Peaks'))
dev.off()

## InCGE vs InMGE --------------------------------------------------------------
cell_type_1 <- 'InCGE'
cell_type_2 <- 'InMGE'

### compare NF factor using RIT RIP nFrags ----------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#ReadsInTSS
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

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

marker_ReadsInTSS <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#ReadsInPeaks
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
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

marker_ReadsInPeaks <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#nFrags
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
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

marker_nFrags <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#Reads_In_Cell_Type_Peaks
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks)/macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks)/macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
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

marker_Reads_In_Cell_Type_Peaks <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#venn plot
marker_ReadsInTSS <- rownames(marker_ReadsInTSS)[marker_ReadsInTSS$log2FC > 1 & marker_ReadsInTSS$fdr < 0.01]
marker_ReadsInPeaks <- rownames(marker_ReadsInPeaks)[marker_ReadsInPeaks$log2FC > 1 & marker_ReadsInPeaks$fdr < 0.01]
marker_nFrags <- rownames(marker_nFrags)[marker_nFrags$log2FC > 1 & marker_nFrags$fdr < 0.01]
marker_Reads_In_Cell_Type_Peaks <- rownames(marker_Reads_In_Cell_Type_Peaks)[marker_Reads_In_Cell_Type_Peaks$log2FC > 1 & marker_Reads_In_Cell_Type_Peaks$fdr < 0.01]

temp <- list(marker_ReadsInTSS = marker_ReadsInTSS,
             marker_ReadsInPeaks = marker_ReadsInPeaks,
             marker_nFrags = marker_nFrags,
             marker_Reads_In_Cell_Type_Peaks = marker_Reads_In_Cell_Type_Peaks)

if(is.null(cell_type_2)){
  cell_type_2 <- 'others'
}
char <- paste(cell_type_1,'vs',cell_type_2,'marker_vennplot_using_RIT_RIP_nFrags_as_NF',sep = ' ')
char <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
char <- gsub(pattern = '-',replacement = '_',x = char,fixed = TRUE)
char <- paste0('./res/step_38_fig_220620/',char,'.pdf')

pdf(file = char,width = 10,height = 10)
ggvenn(data = temp,c('marker_ReadsInTSS','marker_ReadsInPeaks','marker_nFrags','marker_Reads_In_Cell_Type_Peaks'))
dev.off()

## RG vs others --------------------------------------------------------------
cell_type_1 <- 'RG'
cell_type_2 <- NULL

### compare NF factor using RIT RIP nFrags ----------------------------------
match_bias <- ArchR:::.matchBiasCellGroups(input = macaque_multiome_ArchR@cellColData,
                                           groups = macaque_multiome_ArchR$cell_type,
                                           useGroups = cell_type_1,
                                           bgdGroups = cell_type_2,
                                           bias = c("TSSEnrichment", "log10(nFrags)"),
                                           k = 100,n = 1000,bufferRatio = 0.8,seed = 1)

#ReadsInTSS
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

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

marker_ReadsInTSS <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#ReadsInPeaks
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$ReadsInPeaks)/macaque_multiome_ArchR$ReadsInPeaks
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

marker_ReadsInPeaks <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#nFrags
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$nFrags)/macaque_multiome_ArchR$nFrags
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

marker_nFrags <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#Reads_In_Cell_Type_Peaks
mat_1 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$cells]]
mat_2 <- macaque_peak_matrix[,rownames(macaque_multiome_ArchR@cellColData)[match_bias$matchbgd[[1]]$bgd]]

temp <- mat_1@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks)/macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
names(normfactor) <- rownames(macaque_multiome_ArchR@cellColData)
temp <- base::do.call(what = cbind,args = base::lapply(X = colnames(temp),FUN = function(i){
  ii <- temp[,i] * normfactor[i]
  return(ii)
}))
rownames(temp) <- rownames(mat_1@assays$RNA@counts)
colnames(temp) <- colnames(mat_1@assays$RNA@counts)
mat_1 <- temp

temp <- mat_2@assays$RNA@counts
normfactor <- median(macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks)/macaque_multiome_ArchR$Reads_In_Cell_Type_Peaks
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

marker_Reads_In_Cell_Type_Peaks <- my_DF_wilcox_test(mat1 = mat_1,mat2 = mat_2,alternative = 'two.sided',paired = FALSE,workers = 6,future.globals.maxSize = 200*(1024^3))

#venn plot
marker_ReadsInTSS <- rownames(marker_ReadsInTSS)[marker_ReadsInTSS$log2FC > 1 & marker_ReadsInTSS$fdr < 0.01]
marker_ReadsInPeaks <- rownames(marker_ReadsInPeaks)[marker_ReadsInPeaks$log2FC > 1 & marker_ReadsInPeaks$fdr < 0.01]
marker_nFrags <- rownames(marker_nFrags)[marker_nFrags$log2FC > 1 & marker_nFrags$fdr < 0.01]
marker_Reads_In_Cell_Type_Peaks <- rownames(marker_Reads_In_Cell_Type_Peaks)[marker_Reads_In_Cell_Type_Peaks$log2FC > 1 & marker_Reads_In_Cell_Type_Peaks$fdr < 0.01]

temp <- list(marker_ReadsInTSS = marker_ReadsInTSS,
             marker_ReadsInPeaks = marker_ReadsInPeaks,
             marker_nFrags = marker_nFrags,
             marker_Reads_In_Cell_Type_Peaks = marker_Reads_In_Cell_Type_Peaks)

if(is.null(cell_type_2)){
  cell_type_2 <- 'others'
}
char <- paste(cell_type_1,'vs',cell_type_2,'marker_vennplot_using_RIT_RIP_nFrags_as_NF',sep = ' ')
char <- gsub(pattern = ' ',replacement = '_',x = char,fixed = TRUE)
char <- gsub(pattern = '-',replacement = '_',x = char,fixed = TRUE)
char <- paste0('./res/step_38_fig_220620/',char,'.pdf')

pdf(file = char,width = 10,height = 10)
ggvenn(data = temp,c('marker_ReadsInTSS','marker_ReadsInPeaks','marker_nFrags','marker_Reads_In_Cell_Type_Peaks'))
dev.off()