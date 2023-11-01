#####################################################################################
## Project: macaque fetal Brain multiome ATAC_seq                                  ##
## Script Purpose: preprocess of macaque multiome ATAC_seq                         ##
## Data: 2022.03.07                                                                ##
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
setwd('/data/User/sunym/project/Brain/ArchR/')
.libPaths('/data/User/sunym/software/R_lib/R_4.1.3/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(parallel)
library(ArchR)
library(Seurat)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(SingleCellExperiment)
library(dplyr)
library(patchwork)
library(harmony)
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(rtracklayer)
library(org.Mmu.eg.db)
library(clusterProfiler)
library(scibet)
library(viridis)
library(networkD3)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')
#initialize ArchR
addArchRThreads(threads = 5)

# genomeAnnotation <- createGenomeAnnotation(genome = BSgenome.Mmulatta.UCSC.rheMac10)
# 
# geneAnnotation <- rtracklayer::import('/data/User/sunym/project/Brain/data/reference/ensembl_gtf_for_mapping/Macaca_mulatta.Mmul_10.103.gtf')
# geneAnnotation <- rtracklayer::as.data.frame(geneAnnotation)
# geneAnnotation$seqnames <- as.character(geneAnnotation$seqnames)
# geneAnnotation <- geneAnnotation[geneAnnotation$seqnames %in% c(as.character(1:20),'X','Y','MT'),]
# geneAnnotation$seqnames <- paste0('chr',geneAnnotation$seqnames)
# geneAnnotation <- GenomicRanges::makeGRangesFromDataFrame(geneAnnotation,keep.extra.columns = TRUE)
# 
# TxDb_object <- GenomicFeatures::makeTxDbFromGRanges(geneAnnotation)
# geneAnnotation <- createGeneAnnotation(TxDb = TxDb_object, OrgDb = org.Mmu.eg.db)

# create arrow file -------------------------------------------------------
inputFiles <- list.files(path = './fragment_file/macaque_multiome_data/')
temp <- base::lapply(X = inputFiles,FUN = function(x){
  return(strsplit(x = x,split = '_')[[1]][1])
})
temp <- unlist(temp)

inputFiles <- paste('./fragment_file/macaque_multiome_data',inputFiles,sep = '/')
names(inputFiles) <- temp
file.exists(inputFiles)

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4,
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation,
  TileMatParams = list(tileSize = 5000),
  subThreading = FALSE
)

my_send_sms('create arrow file done!')

# doublet inference -------------------------------------------------------
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# create archrproject -----------------------------------------------------
macaque_multiome_ArchR <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = './processed_data/macaque_multiome_ArchR_220307',
  copyArrows = TRUE,
  geneAnnotation = geneAnnotation,
  genomeAnnotation = genomeAnnotation
)

#add metadata
macaque_multiome_ArchR@cellColData$species <- 'macaque'
macaque_multiome_ArchR@cellColData$donor <- as.character(macaque_multiome_ArchR$Sample)
table(macaque_multiome_ArchR$donor)

#filter cells
macaque_integration_Seurat <- readRDS(file = '../processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
table(rownames(macaque_multiome_ArchR@cellColData) %in% colnames(macaque_integration_Seurat))
table(macaque_integration_Seurat$tech)
macaque_multiome_ArchR <- macaque_multiome_ArchR[rownames(macaque_multiome_ArchR@cellColData) %in% colnames(macaque_integration_Seurat)]

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220307/", load = FALSE, overwrite = TRUE)

#QC plot
df <- macaque_multiome_ArchR@cellColData[,c('nFrags','TSSEnrichment','Sample')]
df$nFrags <- log10(df$nFrags)
colnames(df) <- replace(x = colnames(df),list = c(colnames(df) == 'nFrags'),values = 'log10(nFrags)')
df

for (i in unique(df$Sample)) {
  p <- df[df$Sample == i,]
  assign(x = paste('QC',as.character(i),sep = '_'),
         value = ggPoint(
           x = p[,1], 
           y = p[,2], 
           colorDensity = TRUE,
           continuousSet = "sambaNight",
           xlabel = "Log10 Unique Fragments",
           ylabel = "TSS Enrichment",
           title = as.character(i),
           xlim = c(log10(500), quantile(df[,1], probs = 0.99)),
           ylim = c(0, quantile(df[,2], probs = 0.99))) + 
           theme(plot.title = element_text(hjust = 0.5)) + 
           geom_hline(yintercept = 4, lty = "dashed") + 
           geom_vline(xintercept = 3, lty = "dashed"))
}

pdf(file = './res/step_7_fig_220307/TSS_enrichment_vs_nFrag_dotplot_after_filter.pdf',width = 16,height = 5)
QC_A50A+QC_A50B+QC_A82A+QC_A82B+plot_layout(ncol = 4)
dev.off()

#tss enrichment
p <- plotGroups(ArchRProj = macaque_multiome_ArchR, 
                groupBy = "Sample", 
                colorBy = "cellColData", 
                name = "TSSEnrichment",
                plotAs = "violin", 
                alpha = 0.4, 
                addBoxPlot = TRUE) + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12,face = 'bold')) + xlab('') + 
  scale_y_continuous(breaks = seq(0,30,5))

pdf(file = './res/step_7_fig_220307/TSS_enrichment_vlnplot_after_filter.pdf',width = 5,height = 5)
p
dev.off()

#nFrag distribution
p <- plotGroups(ArchRProj = macaque_multiome_ArchR, 
                groupBy = "Sample", 
                colorBy = "cellColData", 
                name = "log10(nFrags)",
                plotAs = "violin", 
                alpha = 0.4, 
                addBoxPlot = TRUE) + 
  theme(aspect.ratio = 1, 
        axis.text = element_text(size = 10), 
        axis.title = element_text(size = 12,face = 'bold')) + xlab('')

pdf(file = './res/step_7_fig_220307/nFrag_vlnplot_after_filter.pdf',width = 5,height = 5)
p
dev.off()

#fragment length distribution
p <- plotFragmentSizes(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample', 
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_7_fig_220307/frag_length_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

#tss accessibility
p <- plotTSSEnrichment(ArchRProj = macaque_multiome_ArchR,
                       groupBy = 'Sample',
                       returnDF = FALSE) + 
  theme(aspect.ratio = 1)

pdf(file = './res/step_7_fig_220307/TSS_accessibility_by_sample_after_filter.pdf',width = 7,height = 5)
p + theme(legend.text = element_text(size = 10), legend.position = 'right')
dev.off()

# dim reduction -----------------------------------------------------------
#LSI
macaque_multiome_ArchR <- addIterativeLSI(ArchRProj = macaque_multiome_ArchR,
                                          useMatrix = "TileMatrix", 
                                          name = "IterativeLSI", 
                                          iterations = 3, 
                                          clusterParams = list(
                                            resolution = c(0.6), 
                                            sampleCells = 10000, 
                                            n.start = 10,
                                            maxClusters = 30
                                          ), 
                                          varFeatures = 20000, 
                                          totalFeatures = 500000, 
                                          dimsToUse = 1:30, 
                                          corCutOff = 0.5, 
                                          force = TRUE)
#cluster
macaque_multiome_ArchR <- addClusters(
  input = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:12, 
  force = TRUE
)

cM <- confusionMatrix(i = macaque_multiome_ArchR$Clusters,j = macaque_multiome_ArchR$Sample)
cM
cM <- cM / Matrix::rowSums(cM)
p <- pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

#umap
macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:12, 
  force = TRUE
)

macaque_multiome_ArchR$cell_type <- as.character(macaque_integration_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"cell_type"])

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220307/", load = FALSE, overwrite = TRUE)

# label transfer ----------------------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_integration_Seurat <- readRDS(file = '../processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_RNA_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'scRNA']
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

#recreate macaque RNA seurat
temp <- macaque_RNA_seurat
macaque_RNA_seurat <- temp@assays$RNA@counts
meta_data <- temp@meta.data
meta_data <- meta_data[,!(colnames(meta_data) %in% c('orig.ident','nCount_RNA','nFeature_RNA'))]
macaque_RNA_seurat <- CreateSeuratObject(counts = macaque_RNA_seurat,project = 'macaque',assay = 'RNA',
                                         meta.data = meta_data,min.cells = 0,min.features = 0)

rm(list = c('temp','meta_data'))
gc()

#constrained integration
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p1+p2+p3+plot_layout(ncol = 3)
#cluster C15 OPC
#cluster C16 RG-2 and RG-3
#cluster C17 RG-1
#cluster C19 Cycling
#cluster C18 IP
#cluster C1 C2 C3 C4 C5 C6 Ex
#cluster C13 End/Per
#cluster C14 Mic
group_list <- SimpleList(
  non_neuron = SimpleList(
    RNA = colnames(macaque_RNA_seurat)[macaque_RNA_seurat$cell_type %in% c('End','Mic','OPC','Per')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C15','C13','C14')]
  ),
  In = SimpleList(
    RNA = colnames(macaque_RNA_seurat)[macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C7','C8','C9','C10','C11','C12')]
  ),
  Ex = SimpleList(
    RNA = colnames(macaque_RNA_seurat)[macaque_RNA_seurat$cell_type %in% c('Ex-1','Ex-2','Ex-3','RG-1','RG-2','RG-3','IP','Cyc-G2M','Cyc-S')],
    ATAC = rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C16','C17','C19','C18','C1','C2','C3','C4','C5','C6')]
  )
)

macaque_multiome_ArchR <- addImputeWeights(macaque_multiome_ArchR)
macaque_multiome_ArchR <- addGeneIntegrationMatrix(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = macaque_RNA_seurat,
  addToArrow = TRUE, 
  groupList = group_list,
  groupRNA = "cell_type",
  nameCell = "predictedCell",
  nameGroup = "predictedGroup",
  nameScore = "predictedScore",
  force = TRUE
)

my_send_sms('ArchR integration done!')

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220307/", load = FALSE, overwrite = TRUE)

#compare with multiome annotation
scibet::Confusion_heatmap(ori = macaque_multiome_ArchR$cell_type,prd = macaque_multiome_ArchR$predictedGroup)

#sankey plot
confusion_matrix <- my_confusion_matrix(ori = paste(macaque_multiome_ArchR$cell_type,'RNA',sep = '.'),prd = paste(macaque_multiome_ArchR$predictedGroup,'ATAC',sep = '.'))
confusion_matrix <- confusion_matrix[confusion_matrix$value > 20,]
nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

p <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                   Source = "IDori", Target = "IDprd",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, nodeWidth= 20, fontSize=13, nodePadding=20)
p

# add gene expression matrix ----------------------------------------------
seRNA <- import10xFeatureMatrix(input = c('/data/User/sunym/project/Brain/data/Multiome/A50A/outs/filtered_feature_bc_matrix.h5',
                                          '/data/User/sunym/project/Brain/data/Multiome/A50B/outs/filtered_feature_bc_matrix.h5',
                                          '/data/User/sunym/project/Brain/data/Multiome/A82A/outs/filtered_feature_bc_matrix.h5',
                                          '/data/User/sunym/project/Brain/data/Multiome/A82B/outs/filtered_feature_bc_matrix.h5'),
                                names = c('A50A','A50B','A82A','A82B'))

temp <- rowRanges(seRNA[[1]])
name_list <- seqlevels(temp)
name_list <- paste0('chr',name_list)
names(name_list) <- seqlevels(temp)
temp <- renameSeqlevels(x = temp,value = name_list)

seRNA <- cbind(assay(seRNA[[1]]),assay(seRNA[[2]]),assay(seRNA[[3]]),assay(seRNA[[4]]))
seRNA <- seRNA[,rownames(macaque_multiome_ArchR@cellColData)]
seRNA <- SummarizedExperiment(assays = list(counts=seRNA),rowRanges = temp)

macaque_multiome_ArchR <- addGeneExpressionMatrix(input = macaque_multiome_ArchR,seRNA = seRNA,verbose = TRUE,force = TRUE)
getAvailableMatrices(macaque_multiome_ArchR)

temp <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneExpressionMatrix',verbose = TRUE)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220307/", load = FALSE, overwrite = TRUE)

# compare multiome and label transfer -------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_integration_Seurat <- readRDS(file = '../processed_data/220305_summary/macaque_integration_Seurat_220307.rds')
macaque_RNA_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')

macaque_multiome_ArchR$Gex_cell_type <- as.character(macaque_RNA_seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"cell_type"])
table(macaque_multiome_ArchR$Gex_cell_type == macaque_multiome_ArchR$cell_type)

macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Gex_cell_type %in% c('Per','End'),"cell_type"] <- 'End/Per'
macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Gex_cell_type %in% c('Cyc-S','Cyc-G2M'),"cell_type"] <- 'Cycling'

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")

pdf(file = './res/step_7_fig_220307/macaque_multiome_RNA_cell_type_label_transfer_dimplot.pdf',width = 16,height = 7)
p1+p2+p3+plot_layout(ncol = 3)
dev.off()

temp <- function(ori,prd){
  cross.validation.filt <- tibble(ori = ori, prd = prd) %>% dplyr::count(ori, prd) %>% tidyr::spread(key = prd, value = n)
  cross.validation.filt[is.na(cross.validation.filt)] = 0
  cross.validation.filt[, -1] <- round(cross.validation.filt[, -1]/rowSums(cross.validation.filt[, -1]), 2)
  cross.validation.filt <- cross.validation.filt %>% tidyr::gather(key = "prd", value = "Prob", -ori)
  return(cross.validation.filt)
}

temp <- temp(ori = macaque_multiome_ArchR$cell_type,prd = macaque_multiome_ArchR$predictedGroup)
temp$ori <- factor(temp$ori,levels = rev(c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-3','RG-2','RG-1','OPC','Ependymal','Mic','End/Per')))
temp$prd <- factor(temp$prd,levels = rev(c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','Cyc-S','Cyc-G2M','RG-3','RG-2','RG-1','OPC','Mic','End','Per')))

pdf(file = './res/step_7_fig_220307/macaque_multiome_RNA_cell_type_label_transfer_confusion_heatmap.pdf',width = 4.5,height = 4.5)
ggplot(data = temp,aes(x=ori,y=prd,fill=Prob)) + 
  geom_tile() + theme(axis.title = element_text(size = 0)) + 
  theme(axis.text = element_text(size = 10)) + 
  theme(legend.title = element_text(size = 0)) + 
  theme(legend.text = element_text(size = 10)) + 
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title = element_blank()) + 
  theme(axis.text.y = element_text(color = "black"), 
        axis.text.x = element_text(color = "black", angle = 45, hjust = 1)) + 
  scale_fill_viridis() + 
  theme(aspect.ratio = 15/14,axis.title = element_text(size = 12,face = 'bold')) + 
  xlab('Gex_cell_type') + ylab('predict label')
dev.off()

# cell type correlation ---------------------------------------------------
#gene express matrix
gene_express_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneExpressionMatrix',verbose = TRUE)
temp <- gene_express_matrix@colData
gene_list <- gene_express_matrix@elementMetadata$name
gene_express_matrix <- gene_express_matrix@assays@data$GeneExpressionMatrix
rownames(gene_express_matrix) <- gene_list
gene_express_matrix <- CreateSeuratObject(counts = gene_express_matrix,project = 'temp',assay = 'RNA',meta.data = as.data.frame(temp),min.cells = 0,min.features = 0)

#gene integratin matrix
gene_integration_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneIntegrationMatrix',verbose = TRUE)
temp <- as.data.frame(gene_integration_matrix@colData)
gene_list <- gene_integration_matrix@elementMetadata$name
gene_list <- gene_list[duplicated(gene_list)]
gene_list <- unique(gene_list)
gene_list <- !(gene_integration_matrix@elementMetadata$name %in% gene_list)
gene_integration_matrix <- gene_integration_matrix[gene_list,]
table(duplicated(gene_integration_matrix@elementMetadata$name))
gene_list <- gene_integration_matrix@elementMetadata$name
gene_integration_matrix <- gene_integration_matrix@assays@data$GeneIntegrationMatrix
rownames(gene_integration_matrix) <- gene_list
gene_integration_matrix <- CreateSeuratObject(counts = gene_integration_matrix,project = 'temp',assay = 'RNA',meta.data = temp,min.cells = 0,min.features = 0)

#create pseudo-bulk matrix
gene_express_matrix <- AverageExpression(object = gene_express_matrix,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
gene_express_matrix <- log1p(gene_express_matrix$RNA)
gene_integration_matrix <- AverageExpression(object = gene_integration_matrix,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
gene_integration_matrix <- log1p(gene_integration_matrix$RNA)

#define gene list
gene_list <- VariableFeatures(macaque_RNA_seurat)
gene_list <- dplyr::intersect(gene_list,rownames(gene_express_matrix))
gene_list <- dplyr::intersect(gene_list,rownames(gene_integration_matrix))

temp <- matrix(nrow = dim(gene_express_matrix)[2],ncol = dim(gene_integration_matrix)[2])
colnames(temp) <- colnames(gene_integration_matrix)
rownames(temp) <- colnames(gene_express_matrix)
temp <- as.data.frame(temp)
for (i in rownames(temp)) {
  for (j in colnames(temp)) {
    temp[i,j] <- cor(x = gene_express_matrix[gene_list,i],y = gene_integration_matrix[gene_list,j],method = 'pearson')
  }
}

cell_list <- c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','Cycling','RG-3','RG-2','RG-1','OPC','Ependymal','Mic','End/Per')

pdf(file = './res/step_7_fig_220307/macaque_multiome_RNA_prediction_expression_cor_on_variable_gene.pdf',width = 6,height = 6)
Heatmap(matrix = temp[cell_list,cell_list],cluster_rows = FALSE,cluster_columns = FALSE,row_title = 'multiome express',column_title = 'prediction express',
        height = unit(4,'inches'),width = unit(4,'inches'),
        row_split = factor(c('In','In','Ex','Ex','Ex','Ex','NPC','NPC','NPC','NPC','NN','NN','NN','NN'),levels = c('In','Ex','NPC','NN')),
        column_split = factor(c('In','In','Ex','Ex','Ex','Ex','NPC','NPC','NPC','NPC','NN','NN','NN','NN'),levels = c('In','Ex','NPC','NN')),
        border = TRUE,name = 'Cor')
dev.off()

# marker gene validation --------------------------------------------------
#marker express
macaque_multiome_ArchR <- addImputeWeights(ArchRProj = macaque_multiome_ArchR,reducedDims = 'IterativeLSI')
marker_gene <- c('PECAM1', #End
                 'PDGFRB', #Per
                 'CX3CR1', #Mic
                 'DNAH11','GJA1','SPATA17','LRRC71','SPAG17','FOXJ1', #Ependymal
                 'SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','NPY','FGFR3','CD9','GPX3', #RG
                 'SOX10','OLIG2','EGFR', #OPC
                 'TOP2A','MKI67','CLSPN','AURKA', #Cyc
                 'EOMES','PPP1R17', #IP
                 'NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2', #Ex
                 'DLX5','GAD2','GAD1','DLX2', #In
                 'LHX6','SST', #InMGE
                 'SP8','NR2F2', #InCGE
                 'MEIS2','ETV1' #InPSB
)

p <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "GeneScoreMatrix", 
  name = marker_gene, 
  embedding = "UMAP",
  imputeWeights = getImputeWeights(macaque_multiome_ArchR)
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

#End checked
p1+p2+p$PECAM1+plot_layout(ncol = 3)

#Per checked
p1+p2+p$PDGFRB+plot_layout(ncol = 3)

#Mic checked
p1+p2+p$CX3CR1+plot_layout(ncol = 3)

#Ependymal not specific
p1+p2+p$FOXJ1+plot_layout(ncol = 3)

#RG checked
p1+p2+p$SOX9+plot_layout(ncol = 3)
p1+p2+p$PAX6+plot_layout(ncol = 3)
p1+p2+p$VIM+plot_layout(ncol = 3)

#OPC checked
p1+p2+p$SOX10+plot_layout(ncol = 3)
p1+p2+p$OLIG2+plot_layout(ncol = 3)
p1+p2+p$EGFR+plot_layout(ncol = 3)

#Cyc surprisingly checked!
p1+p2+p$TOP2A+plot_layout(ncol = 3)
p1+p2+p$MKI67+plot_layout(ncol = 3)

#IP checked
p1+p2+p$EOMES+plot_layout(ncol = 3)
p1+p2+p$PPP1R17+plot_layout(ncol = 3)

#Ex checked
p1+p2+p$NEUROD2+plot_layout(ncol = 3)
p1+p2+p$NEUROD6+plot_layout(ncol = 3)
p1+p2+p$TBR1+plot_layout(ncol = 3)
p1+p2+p$SATB2+plot_layout(ncol = 3)
p1+p2+p$FEZF2+plot_layout(ncol = 3)

#In checked
p1+p2+p$DLX5+plot_layout(ncol = 3)
p1+p2+p$DLX2+plot_layout(ncol = 3)
p1+p2+p$GAD2+plot_layout(ncol = 3)
p1+p2+p$GAD1+plot_layout(ncol = 3)

#InMGE checked
p1+p2+p$LHX6+plot_layout(ncol = 3)
p1+p2+p$SST+plot_layout(ncol = 3)

#InMGE checked!
p1+p2+p$NR2F2+plot_layout(ncol = 3)
p1+p2+p$SP8+plot_layout(ncol = 3)

#all checked except Ependymal, understandable

# call peak ---------------------------------------------------------------
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")
p4 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")

p1+p2+p3+p4+plot_layout(ncol = 2)

#call peak
macaque_multiome_ArchR <- addGroupCoverages(ArchRProj = macaque_multiome_ArchR, 
                                            groupBy = "cell_type", 
                                            minCells = 40, 
                                            maxCells = 500, 
                                            minReplicates = 2, 
                                            maxReplicates = 5, 
                                            sampleRatio = 0.8)


pathToMacs2 <- '/data/User/sunym/env/MACS2/bin/macs2'

macaque_multiome_ArchR <- addReproduciblePeakSet(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  pathToMacs2 = pathToMacs2,
  maxPeaks = 200000,
  genomeSize = 2.7e+09,
  force = TRUE
)

my_send_sms('call peak done!')

getPeakSet(macaque_multiome_ArchR)

#add peak matrix
macaque_multiome_ArchR <- addPeakMatrix(macaque_multiome_ArchR)
getAvailableMatrices(macaque_multiome_ArchR)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220307/", load = FALSE, overwrite = TRUE)

# marker peak and marker gene ---------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')

getAvailableMatrices(ArchRProj = macaque_multiome_ArchR)

#calculate marker peak
markersPeaks <- getMarkerFeatures(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment","log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE)
markerList

peakset <- getPeakSet(macaque_multiome_ArchR)
names(peakset) <- paste(as.character(peakset@seqnames),as.character(peakset@ranges),sep = '-')
table(duplicated(names(peakset)))

#annotate marker peak
for (i in names(markerList)) {
  names(markerList[[i]]) <- paste(as.character(markerList[[i]]@seqnames),as.character(markerList[[i]]@ranges),sep = '-')
  if(sum(names(markerList[[i]]) %in% names(peakset)) < length(names(markerList[[i]]))){
    stop('peak name match error!')
  } else{
    markerList[[i]]$nearestGene <- peakset[names(markerList[[i]])]$nearestGene
  }
  print(paste(as.character(i),'done!',sep = ' '))
}

#get peak express matrix
express_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'PeakMatrix')

meta_data <- express_matrix@colData
meta_data <- as.data.frame(meta_data)
feature_list <- express_matrix@rowRanges
feature_list <- paste(as.character(feature_list@seqnames),as.character(feature_list@ranges),sep = '-')
feature_list[1:10]
table(duplicated(feature_list))

express_matrix <- express_matrix@assays@data$PeakMatrix
rownames(express_matrix) <- feature_list
temp <- Seurat::CreateSeuratObject(counts = express_matrix,project = 'temp',assay = 'RNA',meta.data = meta_data)
temp <- AverageExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts',verbose = TRUE)
temp <- temp$RNA
temp[1:3,1:3]
express_matrix <- temp

markerList <- markerList[c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','RG-3','OPC')]
feature_list <- base::lapply(names(markerList),FUN = function(x){
  temp_marker <- markerList[[x]]
  temp_marker <- paste(as.character(temp_marker@seqnames),as.character(temp_marker@ranges),sep = '-')
  return(temp_marker)
})
for (i in 1:length(feature_list)) {
  names(feature_list[[i]]) <- rep(names(markerList)[i],times = length(feature_list[[i]]))
}
feature_list <- unlist(feature_list)

table(feature_list %in% rownames(express_matrix))

express_matrix <- express_matrix[feature_list,]
express_matrix <- log1p(express_matrix)
express_matrix <- t(scale(t(express_matrix)))
express_matrix <- express_matrix[,names(markerList)]
express_matrix <- t(express_matrix)

pdf(file = './res/step_7_fig_220307/macaque_multiome_ArchR_marker_peak_heatmap.pdf',width = 12,height = 6)
Heatmap(matrix = express_matrix,cluster_columns = FALSE,cluster_rows = FALSE,show_row_names = TRUE,show_column_names = FALSE,
        name = 'z score',width = unit(10,'inches'),height = unit(5,'inches'),
        column_split = factor(names(feature_list),levels = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','RG-3','OPC')),
        border = TRUE,row_split = factor(c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','RG-3','OPC'),
                                         levels = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','RG-3','OPC')))
dev.off()

# peak co-accessibility ----------------------------------------------------
macaque_multiome_ArchR <- addCoAccessibility(
  ArchRProj = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.5,
  resolution = 1,
  returnLoops = FALSE
)
cA

markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getCoAccessibility(macaque_multiome_ArchR,resolution = 2000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG-2','RG-1')
)

#RG
grid::grid.newpage()
grid::grid.draw(p$SOX9)

grid::grid.newpage()
grid::grid.draw(p$HOPX)

#OPC
grid::grid.newpage()
grid::grid.draw(p$OLIG2)

grid::grid.newpage()
grid::grid.draw(p$EGFR)

#IP
grid::grid.newpage()
grid::grid.draw(p$EOMES)

grid::grid.newpage()
grid::grid.draw(p$PPP1R17)

#Ex
grid::grid.newpage()
grid::grid.draw(p$NEUROD2)

grid::grid.newpage()
grid::grid.draw(p$NEUROD6)

grid::grid.newpage()
grid::grid.draw(p$SATB2)

grid::grid.newpage()
grid::grid.draw(p$SATB2)

#In
grid::grid.newpage()
grid::grid.draw(p$DLX5)

grid::grid.newpage()
grid::grid.draw(p$GAD2)

grid::grid.newpage()
grid::grid.draw(p$SST)

grid::grid.newpage()
grid::grid.draw(p$SP8)

grid::grid.newpage()
grid::grid.draw(p$LHX6)

# peak and gene express correlation ---------------------------------------
macaque_multiome_ArchR <- addPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI",
  useMatrix = 'GeneExpressionMatrix'
)

p2g <- getPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g

markerGenes  <- c(
  'SOX9','PAX6','VIM','HOPX', #RG
  'OLIG2','EGFR', #OPC
  'EOMES','PPP1R17', #IP
  'NEUROD2','NEUROD6','SATB2','FEZF2', #Ex
  'DLX5','GAD2','SST','SP8','LHX6' #In
)

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = markerGenes, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','OPC','RG-3','RG-2','RG-1')
)

grid::grid.newpage()
grid::grid.draw(p$SOX9)

grid::grid.newpage()
grid::grid.draw(p$PAX6)

grid::grid.newpage()
grid::grid.draw(p$HOPX)

grid::grid.newpage()
grid::grid.draw(p$OLIG2)

grid::grid.newpage()
grid::grid.draw(p$EGFR)

grid::grid.newpage()
grid::grid.draw(p$EOMES)

grid::grid.newpage()
grid::grid.draw(p$NEUROD2)

grid::grid.newpage()
grid::grid.draw(p$NEUROD6)

grid::grid.newpage()
grid::grid.draw(p$DLX5)

p <- plotPeak2GeneHeatmap(ArchRProj = macaque_multiome_ArchR, groupBy = "cell_type", k = 10)
ComplexHeatmap::draw(p, heatmap_legend_side = "bot", annotation_legend_side = "bot", width = unit(2,'inches'), heigh = unit(8,'inches'))

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220307/", load = FALSE, overwrite = TRUE)


# cell type marker for FISH -----------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220307/Save-ArchR-Project.rds')
macaque_integration_Seurat <- readRDS(file = '../processed_data/220305_summary/macaque_integration_Seurat_220307.rds')

#dimplot
p1 <- my_dimplot(embedding = macaque_integration_Seurat@reductions$harmonyUMAP@cell.embeddings,meta_data = macaque_integration_Seurat@meta.data,group.by = 'cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1.5),
        axis.title = element_blank()) + labs(title = 'macaque RNA')

p2 <- my_dimplot(embedding = macaque_multiome_ArchR@embeddings$UMAP$df,meta_data = macaque_multiome_ArchR@cellColData,group.by = 'Gex_cell_type',label = TRUE,repel = TRUE) + 
  theme(aspect.ratio = 1,
        legend.position = 'none',
        axis.line = element_blank(),
        panel.background = element_rect(fill = NA,color = 'black',size = 1.5),
        axis.title = element_blank()) + labs(title = 'macaque multiome ATAC')

pdf(file = './res/step_7_fig_220307/macaque_RNA_multiome_dimplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#marker peak
markersPeaks <- getMarkerFeatures(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "cell_type",
  bias = c("TSSEnrichment","log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1",returnGR = TRUE)
markerList

## RG-1 --------------------------------------------------------------------
RG_1_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-1',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_1_marker_list <- RG_1_marker_list[order(RG_1_marker_list$pct.1,decreasing = TRUE),]
RG_1_marker_list <- RG_1_marker_list[RG_1_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_1_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

RG_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-1',ident.2 = c('RG-2','RG-3'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_subset_marker <- RG_subset_marker[order(RG_subset_marker$pct.1,decreasing = TRUE),]
RG_subset_marker <- RG_subset_marker[RG_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect RG-1 marker
gene_list <- dplyr::intersect(rownames(RG_1_marker_list),rownames(RG_subset_marker))
VlnPlot(object = macaque_integration_Seurat,features = gene_list[1:12],assay = 'RNA',group.by = 'cell_type',pt.size = 0,slot = 'data')

saveRDS(gene_list,file = './res/step_7_fig_220307/macaque_integration_Seurat_RG_1_marker_list.rds')

#integrate gene express and peak
p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = gene_list, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

grid::grid.newpage()
grid::grid.draw(p[[gene_list[43]]])

#PLEKHA7 MOB3B MMS22L GPR156 CLU TULP3 AMMECR1 SLC15A2 SFRP1 FSTL1
#SFRP1 WNT

pdf(file = './res/step_7_fig_220307/macaque_integration_Seurat_RG_1_marker_featureplot.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = c('PLEKHA7','MOB3B','MMS22L','GPR156','CLU','TULP3','AMMECR1','SLC15A2','SFRP1','FSTL1'),object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.2,0.99),print = TRUE,point.size = 0.5)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('PLEKHA7','MOB3B','MMS22L','GPR156','CLU','TULP3','AMMECR1','SLC15A2','SFRP1','FSTL1'), 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

for(i in c('PLEKHA7','MOB3B','MMS22L','GPR156','CLU','TULP3','AMMECR1','SLC15A2','SFRP1','FSTL1')){
  char <- paste0('./res/step_7_fig_220307/macaque_integration_Seurat_RG_1_',i,'_trackplot.pdf')
  pdf(file = char,width = 12,height = 6)
  grid::grid.newpage()
  grid::grid.draw(p[[i]])
  dev.off()
}

## RG-2 --------------------------------------------------------------------
RG_2_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-2',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_2_marker_list <- RG_2_marker_list[order(RG_2_marker_list$pct.1,decreasing = TRUE),]
RG_2_marker_list <- RG_2_marker_list[RG_2_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_2_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

RG_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-2',ident.2 = c('RG-1','RG-3'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_subset_marker <- RG_subset_marker[order(RG_subset_marker$pct.1,decreasing = TRUE),]
RG_subset_marker <- RG_subset_marker[RG_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect RG-2 marker
gene_list <- dplyr::intersect(rownames(RG_2_marker_list),rownames(RG_subset_marker))
saveRDS(gene_list,file = './res/step_7_fig_220307/macaque_integration_Seurat_RG_2_marker_list.rds')

#integrate gene express and peak
p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = gene_list, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

grid::grid.newpage()
grid::grid.draw(p[[gene_list[57]]])

#ITPR2 PREX1 CAMK2G ARHGEF4 SHISA6 CA8 ITGA6 TIMP3 RAB31 AQP4
pdf(file = './res/step_7_fig_220307/macaque_integration_Seurat_RG_2_marker_featureplot.pdf',width = 12,height = 8)
geneSetAveragePlot(genes = c('ITPR2','PREX1','CAMK2G','ARHGEF4','SHISA6','CA8','ITGA6','TIMP3','RAB31','AQP4'),object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.2,0.99),print = TRUE,point.size = 0.5)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('ITPR2','PREX1','CAMK2G','ARHGEF4','SHISA6','CA8','ITGA6','TIMP3','RAB31','AQP4'), 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

for(i in c('ITPR2','PREX1','CAMK2G','ARHGEF4','SHISA6','CA8','ITGA6','TIMP3','RAB31','AQP4')){
  char <- paste0('./res/step_7_fig_220307/macaque_integration_Seurat_RG_2_',i,'_trackplot.pdf')
  pdf(file = char,width = 12,height = 6)
  grid::grid.newpage()
  grid::grid.draw(p[[i]])
  dev.off()
}

## RG-3 --------------------------------------------------------------------
RG_3_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-3',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_3_marker_list <- RG_3_marker_list[order(RG_3_marker_list$pct.1,decreasing = TRUE),]
RG_3_marker_list <- RG_3_marker_list[RG_3_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_3_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

RG_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'RG-3',ident.2 = c('RG-1','RG-2'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
RG_subset_marker <- RG_subset_marker[order(RG_subset_marker$pct.1,decreasing = TRUE),]
RG_subset_marker <- RG_subset_marker[RG_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(RG_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect RG-3 marker
gene_list <- dplyr::intersect(rownames(RG_3_marker_list),rownames(RG_subset_marker))
saveRDS(gene_list,file = './res/step_7_fig_220307/macaque_integration_Seurat_RG_3_marker_list.rds')

#integrate gene express and peak
p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = gene_list, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

grid::grid.newpage()
grid::grid.draw(p[[gene_list[49]]])

#SH3D19 CDK6 MAP3K1 CEP83 PCDH15
pdf(file = './res/step_7_fig_220307/macaque_integration_Seurat_RG_3_marker_featureplot.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('SH3D19','CDK6','MAP3K1','CEP83'),object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.2,0.99),print = TRUE,point.size = 0.5)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('SH3D19','CDK6','MAP3K1','CEP83'), 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

for(i in c('SH3D19','CDK6','MAP3K1','CEP83')){
  char <- paste0('./res/step_7_fig_220307/macaque_integration_Seurat_RG_3_',i,'_trackplot.pdf')
  pdf(file = char,width = 12,height = 6)
  grid::grid.newpage()
  grid::grid.draw(p[[i]])
  dev.off()
}

## Ex-1 --------------------------------------------------------------------
Ex_1_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-1',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
Ex_1_marker_list <- Ex_1_marker_list[order(Ex_1_marker_list$pct.1,decreasing = TRUE),]
Ex_1_marker_list <- Ex_1_marker_list[Ex_1_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(Ex_1_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

Ex_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-1',ident.2 = c('Ex-3','Ex-2'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
Ex_subset_marker <- Ex_subset_marker[order(Ex_subset_marker$pct.1,decreasing = TRUE),]
Ex_subset_marker <- Ex_subset_marker[Ex_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(Ex_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect Ex-1 marker
gene_list <- dplyr::intersect(rownames(Ex_1_marker_list),rownames(Ex_subset_marker))
saveRDS(gene_list,file = './res/step_7_fig_220307/macaque_integration_Seurat_Ex_1_marker_list.rds')

#integrate gene express and peak
p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = gene_list, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

grid::grid.newpage()
grid::grid.draw(p[[gene_list[47]]])

#LRP8 PRSS12
pdf(file = './res/step_7_fig_220307/macaque_integration_Seurat_Ex_1_marker_featureplot.pdf',width = 12,height = 3)
geneSetAveragePlot(genes = c('LRP8','PRSS12'),object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.2,0.99),print = TRUE,point.size = 0.5)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('LRP8','PRSS12'), 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

for(i in c('LRP8','PRSS12')){
  char <- paste0('./res/step_7_fig_220307/macaque_integration_Seurat_Ex_1_',i,'_trackplot.pdf')
  pdf(file = char,width = 12,height = 6)
  grid::grid.newpage()
  grid::grid.draw(p[[i]])
  dev.off()
}

## Ex-2 --------------------------------------------------------------------
Ex_2_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-2',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
Ex_2_marker_list <- Ex_2_marker_list[order(Ex_2_marker_list$pct.1,decreasing = TRUE),]
Ex_2_marker_list <- Ex_2_marker_list[Ex_2_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(Ex_2_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

Ex_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-2',ident.2 = c('Ex-3','Ex-1'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
Ex_subset_marker <- Ex_subset_marker[order(Ex_subset_marker$pct.1,decreasing = TRUE),]
Ex_subset_marker <- Ex_subset_marker[Ex_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(Ex_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect Ex_2 marker
gene_list <- dplyr::intersect(rownames(Ex_2_marker_list),rownames(Ex_subset_marker))
saveRDS(gene_list,file = './res/step_7_fig_220307/macaque_integration_Seurat_Ex_2_marker_list.rds')

#integrate gene express and peak
p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = gene_list, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

grid::grid.newpage()
grid::grid.draw(p[[gene_list[98]]])

#DOK5 NWD2 KLHL1 GRIA4 CCDC85A ARHGAP20 ARFGEF3 LRFN2 KIF16B SSX2IP CDH8 ACTN2 POSTN CA10 MYO5B GLIS3
pdf(file = './res/step_7_fig_220307/macaque_integration_Seurat_Ex_2_marker_featureplot.pdf',width = 12,height = 12)
geneSetAveragePlot(genes = c('DOK5','NWD2','KLHL1','GRIA4','CCDC85A','ARHGAP20','ARFGEF3','LRFN2','KIF16B','SSX2IP','CDH8','ACTN2','POSTN','CA10','MYO5B','GLIS3'),object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.2,0.99),print = TRUE,point.size = 0.5)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('DOK5','NWD2','KLHL1','GRIA4','CCDC85A','ARHGAP20','ARFGEF3','LRFN2','KIF16B','SSX2IP','CDH8','ACTN2','POSTN','CA10','MYO5B','GLIS3'), 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

for(i in c('DOK5','NWD2','KLHL1','GRIA4','CCDC85A','ARHGAP20','ARFGEF3','LRFN2','KIF16B','SSX2IP','CDH8','ACTN2','POSTN','CA10','MYO5B','GLIS3')){
  char <- paste0('./res/step_7_fig_220307/macaque_integration_Seurat_Ex_2_',i,'_trackplot.pdf')
  pdf(file = char,width = 12,height = 6)
  grid::grid.newpage()
  grid::grid.draw(p[[i]])
  dev.off()
}

## Ex-3 --------------------------------------------------------------------
Ex_3_marker_list <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-3',group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
Ex_3_marker_list <- Ex_3_marker_list[order(Ex_3_marker_list$pct.1,decreasing = TRUE),]
Ex_3_marker_list <- Ex_3_marker_list[Ex_3_marker_list$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(Ex_3_marker_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

Ex_subset_marker <- FindMarkers(object = macaque_integration_Seurat,ident.1 = 'Ex-3',ident.2 = c('Ex-2','Ex-1'),group.by = 'cell_type',
                                assay = 'RNA',slot = 'data',logfc.threshold = 0.5,test.use = 'bimod',min.pct = 0.4,
                                verbose = TRUE,only.pos = TRUE)
Ex_subset_marker <- Ex_subset_marker[order(Ex_subset_marker$pct.1,decreasing = TRUE),]
Ex_subset_marker <- Ex_subset_marker[Ex_subset_marker$pct.2 < 0.4,]
VlnPlot(object = macaque_integration_Seurat,features = rownames(Ex_subset_marker)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type')

#intersect Ex-3 marker
gene_list <- dplyr::intersect(rownames(Ex_3_marker_list),rownames(Ex_subset_marker))
saveRDS(gene_list,file = './res/step_7_fig_220307/macaque_integration_Seurat_Ex_3_marker_list.rds')

#integrate gene express and peak
p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = gene_list, 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

grid::grid.newpage()
grid::grid.draw(p[[gene_list[88]]])

#HS3ST4 KIAA1217 TSHZ3 TMEM178A OSBPL10 RALGPS2 RXFP1
pdf(file = './res/step_7_fig_220307/macaque_integration_Seurat_Ex_3_marker_featureplot.pdf',width = 12,height = 6)
geneSetAveragePlot(genes = c('HS3ST4','KIAA1217','TSHZ3','TMEM178A','OSBPL10','RALGPS2','RXFP1'),object = macaque_integration_Seurat,object.class = 'seurat',assay = 'RNA',embedding = 'harmonyUMAP',plot.type = 'panels',
                   scaled = FALSE,color.palette = c('lightgrey','blue'),aspectratio = 1,trim = c(0.2,0.99),print = TRUE,point.size = 0.5)
dev.off()

p <- plotBrowserTrack(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "cell_type", 
  geneSymbol = c('HS3ST4','KIAA1217','TSHZ3','TMEM178A','OSBPL10','RALGPS2','RXFP1'), 
  upstream = 50000,
  downstream = 50000,
  loops = getPeak2GeneLinks(macaque_multiome_ArchR,resolution = 1000),
  useGroups = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-3','RG-2','RG-1','Cycling','OPC','Mic','End/Per')
)

for(i in c('HS3ST4','KIAA1217','TSHZ3','TMEM178A','OSBPL10','RALGPS2','RXFP1')){
  char <- paste0('./res/step_7_fig_220307/macaque_integration_Seurat_Ex_3_',i,'_trackplot.pdf')
  pdf(file = char,width = 12,height = 6)
  grid::grid.newpage()
  grid::grid.draw(p[[i]])
  dev.off()
}
