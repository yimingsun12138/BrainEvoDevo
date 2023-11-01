#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: Using bulk RNAseq to define brain layer marker                  ##
## Data: 2021.06.21                                                                ##
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

#library
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

#load data
macaque_express_matrix <- readRDS(file = './processed_data/macaque_bulk_RNA_seq_from_liuyt_210621.rds')
macaque_coldata <- data.frame(species='macaque',sample=colnames(macaque_express_matrix))
rownames(macaque_coldata) <- colnames(macaque_express_matrix)
donor_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][1]
  return(temp)
})
donor_list <- unlist(donor_list)
macaque_coldata$donor <- donor_list
layer_list <- base::lapply(macaque_coldata$sample,function(x){
  temp <- strsplit(x,split = '_')
  temp <- temp[[1]][2]
  return(temp)
})
layer_list <- unlist(layer_list)
macaque_coldata$layer <- layer_list

#coldata for DE analysis
layer_list <- unique(layer_list)
temp <- do.call(cbind,base::lapply(layer_list,function(x){
  temp_list <- as.character(macaque_coldata$layer == x)
  temp_list[which(temp_list == 'TRUE')] <- 'YES'
  temp_list[which(temp_list == 'FALSE')] <- 'NO'
  return(temp_list)
}))
colnames(temp) <- layer_list
rownames(temp) <- rownames(macaque_coldata)
macaque_coldata <- cbind(macaque_coldata,temp)

#DESeq2 for all layer
dds <- DESeqDataSetFromMatrix(countData = macaque_express_matrix,colData = macaque_coldata,design = ~donor+layer)
rld <- rlog(dds, blind = FALSE)
#Sample distances
sampleDists <- dist(t(assay(rld)))
sampleDists
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- rld$sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9,"Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colors,
         main = 'sample euclidean distance',
         clustering_method = 'ward.D')
#PCAplot
pcaData <- plotPCA(rld, intgroup = c('donor','layer'), returnData = TRUE)
pcaData
percentVar <- round(100 * attr(pcaData, "percentVar"))
percentVar
ggplot(pcaData, aes(x = PC1, y = PC2, color = layer, shape = donor)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("PCA with VST data")

#find DE
layer_list <- unique(macaque_coldata$layer)
DE_list <- base::lapply(layer_list,function(x){
  formula_char <- paste0('~','donor','+',as.character(x))
  formula_char <- as.formula(formula_char)
  temp <- DESeqDataSetFromMatrix(countData = macaque_express_matrix,colData = macaque_coldata,design = formula_char)
  temp <- DESeq(temp)
  temp <- results(temp,contrast = c(as.character(x),'YES','NO'),lfcThreshold = 0,alpha = 0.05)
  temp <- temp[!is.na(temp$padj),]
  temp <- temp[!is.na(temp$log2FoldChange),]
  temp <- temp[(temp$log2FoldChange > 0) & (temp$padj < 0.05),]
  return(rownames(temp))
})
names(DE_list) <- layer_list
saveRDS(DE_list,file = './res/step_2_fig_210622/macaque_layer_DE_list_210622.rds')
