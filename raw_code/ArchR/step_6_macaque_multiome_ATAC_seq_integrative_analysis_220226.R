#####################################################################################
## Project: macaque fetal Brain multiome ATAC_seq                                  ##
## Script Purpose: macaque multiome ATAC_seq integrative analysis                  ##
## Data: 2022.02.26                                                                ##
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
.libPaths(c('/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.3/',
            '/data/User/sunym/software/R_lib/R_4.1.3/'))
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
library(harmony,lib.loc = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.3/')
library(BSgenome.Mmulatta.UCSC.rheMac10)
library(rtracklayer)
library(org.Mmu.eg.db)
library(clusterProfiler)
library(scibet)
library(viridis)
library(networkD3)
library(ggpointdensity)
library(ggvenn)

#source list
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')
#initialize ArchR
addArchRThreads(threads = 5)

# label transfer using macaque_RNA_seurat ---------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(
  ArchRProj = macaque_multiome_ArchR, 
  colorBy = "cellColData", 
  name = 'DoubletEnrichment', 
  embedding = "UMAP",
  imputeWeights = NULL
)

p1+p2+p3+plot_layout(ncol = 3)

macaque_RNA_seurat <- readRDS(file = '../processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_220226.rds')
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$tech == 'scRNA']
macaque_RNA_seurat <- macaque_RNA_seurat[,!(macaque_RNA_seurat$cell_type %in% c('Cyc-S','Cyc-G2M'))]
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

# #unconstrained integration
# macaque_multiome_ArchR <- addImputeWeights(macaque_multiome_ArchR)
# macaque_multiome_ArchR <- addGeneIntegrationMatrix(
#   ArchRProj = macaque_multiome_ArchR, 
#   useMatrix = "GeneScoreMatrix",
#   matrixName = "GeneIntegrationMatrix",
#   reducedDims = "IterativeLSI",
#   seRNA = macaque_RNA_seurat,
#   addToArrow = FALSE,
#   groupRNA = "cell_type",
#   nameCell = "predictedCell_Un",
#   nameGroup = "predictedGroup_Un",
#   nameScore = "predictedScore_Un",
#   force = TRUE
# )
# 
# my_send_sms("intergration done!")
# 
# p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
# p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup_Un", embedding = "UMAP")
# 
# p1+p2+plot_layout(ncol = 2)

#constrained integration
group_list <- SimpleList(
  non_neuron <- SimpleList(
    RNA = colnames(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('End','Ependymal','Mic','OPC','Per')]),
    ATAC = rownames(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Clusters %in% c('C3','C1','C2'),])
  ),
  In_neuron <- SimpleList(
    RNA = colnames(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('InCGE','InMGE')]),
    ATAC = rownames(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Clusters %in% c('C20','C21','C22','C18','C19','C15','C16','C17'),])
  ),
  Ex_neuron <- SimpleList(
    RNA = colnames(macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('RG-1','RG-2','RG-3','IP','Ex-1','Ex-2','Ex-3','Ex-4')]),
    ATAC = rownames(macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Clusters %in% c('C7','C4','C5','C6','C8','C11','C10','C12','C14','C13','C9'),])
  )
)

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

p1+p2+plot_layout(ncol = 2)

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220216/", load = FALSE, overwrite = TRUE)

#compare with multiome annotation
macaque_RNA_seurat <- readRDS(file = '../processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_220226.rds')
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$tech == 'multiome']

macaque_multiome_ArchR@cellColData[colnames(macaque_RNA_seurat),'cell_type'] <- macaque_RNA_seurat$cell_type

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "cell_type", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")

p1+p2+plot_layout(ncol = 2)

scibet::Confusion_heatmap(ori = macaque_multiome_ArchR$cell_type,prd = macaque_multiome_ArchR$predictedGroup)

#sankey plot
confusion_matrix <- my_confusion_matrix(ori = paste(macaque_multiome_ArchR$cell_type,'RNA',sep = '.'),prd = paste(macaque_multiome_ArchR$predictedGroup,'ATAC',sep = '.'))
confusion_matrix <- confusion_matrix[confusion_matrix$value > 40,]
nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

p <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                   Source = "IDori", Target = "IDprd",
                   Value = "value", NodeID = "name", 
                   sinksRight=FALSE, nodeWidth= 20, fontSize=13, nodePadding=20)


# add gene expression matrix ----------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')

seRNA <- import10xFeatureMatrix(input = c('/data/User/sunym/project/Brain/data/Multiome/A50A/outs/filtered_feature_bc_matrix.h5',
                                          '/data/User/sunym/project/Brain/data/Multiome/A50B/outs/filtered_feature_bc_matrix.h5',
                                          '/data/User/sunym/project/Brain/data/Multiome/A82A/outs/filtered_feature_bc_matrix.h5',
                                          '/data/User/sunym/project/Brain/data/Multiome/A82B/outs/filtered_feature_bc_matrix.h5'),
                                names = c('A50A','A50B','A82A','A82B'))
#something wrong!

temp <- rowRanges(seRNA[[1]])
name_list <- seqlevels(temp)
name_list <- paste0('chr',name_list)
names(name_list) <- seqlevels(temp)
temp <- renameSeqlevels(x = temp,value = name_list)

seRNA <- cbind(assay(seRNA[[1]]),assay(seRNA[[2]]),assay(seRNA[[3]]),assay(seRNA[[4]]))
seRNA <- seRNA[,rownames(macaque_multiome_ArchR@cellColData)]
seRNA <- SummarizedExperiment(assays = list(counts=seRNA),rowRanges = temp)

macaque_multiome_ArchR <- addGeneExpressionMatrix(input = macaque_multiome_ArchR,seRNA = seRNA,verbose = TRUE,force = TRUE,threads = 1)

my_send_sms('add gene expression matrix done!')

#save data
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220216/", load = FALSE, overwrite = TRUE)

# compare multiome and label transfer -------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')
macaque_RNA_seurat <- readRDS(file = '../processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_220226.rds')
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$tech == 'multiome']
macaque_multiome_ArchR@cellColData[colnames(macaque_RNA_seurat),'RNA_cell_type'] <- macaque_RNA_seurat$cell_type

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "RNA_cell_type", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")

pdf(file = './res/step_6_fig_220227/macaque_multiome_RNA_cell_type_label_transfer_dimplot.pdf',width = 12,height = 7)
p1+p2+plot_layout(ncol = 2)
dev.off()

temp <- macaque_multiome_ArchR$predictedGroup
temp[temp %in% c('End','Per','Mic')] <- 'End/Per/Mic'

pdf(file = './res/step_6_fig_220227/macaque_multiome_RNA_cell_type_label_transfer_confusion_heatmap.pdf',width = 8,height = 6)
scibet::Confusion_heatmap(ori = macaque_multiome_ArchR$RNA_cell_type,prd = temp) + 
  theme(aspect.ratio = 13/17,
        axis.title = element_text(size = 14,face = 'bold')) + 
  xlab('RNA cell type') + ylab('ATAC prediction')
dev.off()


## cell type cor ----------------------------------------------------------------
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
gene_express_matrix <- AverageExpression(object = gene_express_matrix,assays = 'RNA',return.seurat = FALSE,group.by = 'RNA_cell_type',slot = 'counts',verbose = TRUE)
gene_express_matrix <- log1p(gene_express_matrix$RNA)
gene_integration_matrix <- AverageExpression(object = gene_integration_matrix,assays = 'RNA',return.seurat = FALSE,group.by = 'RNA_cell_type',slot = 'counts',verbose = TRUE)
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

cell_list <- c('InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','RG-3','RG-2','RG-1','OPC','Ependymal','Mic','Per','End')

pdf(file = './res/step_6_fig_220227/macaque_multiome_RNA_prediction_expression_cor_on_variable_gene.pdf',width = 6,height = 6)
Heatmap(matrix = temp[cell_list,cell_list],cluster_rows = FALSE,cluster_columns = FALSE,row_title = 'multiome express',column_title = 'prediction express',
        height = unit(4,'inches'),width = unit(4,'inches'),
        row_split = factor(c('In','In','Ex','Ex','Ex','Ex','Ex','NPC','NPC','NPC','NPC','NPC','NN','NN','NN','NN','NN'),levels = c('In','Ex','NPC','NN')),
        column_split = factor(c('In','In','Ex','Ex','Ex','Ex','Ex','NPC','NPC','NPC','NPC','NPC','NN','NN','NN','NN','NN'),levels = c('In','Ex','NPC','NN')),
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
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "RNA_cell_type", embedding = "UMAP")

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
#suggest:aggregate Cyc-S and Cyc-G2M


# call peak ---------------------------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")

p1+p2+plot_layout(ncol = 2)

#add metadata
macaque_integration_Seurat <- readRDS(file = '../processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_annotated_220302.rds')
macaque_multiome_ArchR$Gex_cell_type <- macaque_integration_Seurat@meta.data[rownames(macaque_multiome_ArchR@cellColData),"cell_type"]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

macaque_multiome_ArchR@cellColData[macaque_multiome_ArchR$Gex_cell_type %in% c('Cyc-S','Cyc-G2M'),"Gex_cell_type"] <- 'Cyc'

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#call peak
macaque_multiome_ArchR <- addGroupCoverages(ArchRProj = macaque_multiome_ArchR, 
                                            groupBy = "Gex_cell_type", 
                                            minCells = 40, 
                                            maxCells = 500, 
                                            minReplicates = 2, 
                                            maxReplicates = 5, 
                                            sampleRatio = 0.8)


pathToMacs2 <- '/data/User/sunym/env/MACS2/bin/macs2'

macaque_multiome_ArchR <- addReproduciblePeakSet(
  ArchRProj = macaque_multiome_ArchR, 
  groupBy = "Gex_cell_type", 
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
saveArchRProject(ArchRProj = macaque_multiome_ArchR, outputDirectory = "./processed_data/macaque_multiome_ArchR_220216/", load = FALSE, overwrite = TRUE)

# marker peak and marker gene ---------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')
getAvailableMatrices(ArchRProj = macaque_multiome_ArchR)

#calculate marker peak
markersPeaks <- getMarkerFeatures(
  ArchRProj = macaque_multiome_ArchR, 
  useMatrix = "PeakMatrix", 
  groupBy = "Gex_cell_type",
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
temp <- AverageExpression(object = temp,assays = 'RNA',return.seurat = FALSE,group.by = 'Gex_cell_type',slot = 'counts',verbose = TRUE)
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

Heatmap(matrix = express_matrix,cluster_columns = FALSE,cluster_rows = FALSE,show_row_names = TRUE,show_column_names = FALSE,
        name = 'z score',width = unit(10,'inches'),height = unit(5,'inches'),
        column_split = factor(names(feature_list),levels = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','RG-3','OPC')),
        border = TRUE,row_split = factor(c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','RG-3','OPC'),
                                         levels = c('InCGE','InMGE','Ex-3','Ex-2','Ex-1','IP','RG-1','RG-2','RG-3','OPC')))

# intersect with liuyt cell list ------------------------------------------
#load data
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')
cell_list <- readRDS(file = './res/step_6_fig_220227/liuyt_multiome_cell_list.rds')
cell_list <- paste0('A',cell_list)

temp <- list(sunym=rownames(macaque_multiome_ArchR@cellColData),liuyt=cell_list)
ggvenn(temp,columns = c('sunym','liuyt'))

#filter data
Cycling_cell_list <- rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Gex_cell_type == 'Cyc']
Cycling_cell_list <- Cycling_cell_list[!(Cycling_cell_list %in% cell_list)]
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(rownames(macaque_multiome_ArchR@cellColData) %in% Cycling_cell_list)]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#not satisfying
# filter In progenitor ----------------------------------------------------
macaque_multiome_ArchR <- readRDS(file = './processed_data/macaque_multiome_ArchR_220216/Save-ArchR-Project.rds')
p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#filter cluster 20 and 18
table(macaque_multiome_ArchR[macaque_multiome_ArchR$Clusters %in% c('C20','C18')]$Gex_cell_type)
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(macaque_multiome_ArchR$Clusters %in% c('C20','C18'))]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "predictedGroup", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#redo LSI, cluster and umap
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
  dimsToUse = 1:15, 
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

#cluster 20 seems weird

#umap
macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:15, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#filter all other cells in In
filter_cell_list <- rownames(macaque_multiome_ArchR@cellColData)[macaque_multiome_ArchR$Clusters %in% c('C9','C10','C11','C12','C13') & !(macaque_multiome_ArchR$Gex_cell_type %in% c('InCGE','InMGE'))]
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(rownames(macaque_multiome_ArchR@cellColData) %in% filter_cell_list)]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#cluster 3 seems weird
table(macaque_multiome_ArchR[macaque_multiome_ArchR$Clusters == 'C3']$Gex_cell_type)
macaque_multiome_ArchR <- macaque_multiome_ArchR[!(macaque_multiome_ArchR$Clusters %in% c('C3'))]

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#redo
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

macaque_multiome_ArchR <- addClusters(
  input = macaque_multiome_ArchR,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  maxClusters = 30, 
  dimsToUse = 1:22, 
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

macaque_multiome_ArchR <- addUMAP(
  ArchRProj = macaque_multiome_ArchR, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 50, 
  minDist = 0.6, 
  metric = "cosine", 
  dimsToUse = 1:22, 
  force = TRUE
)

p1 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
p3 <- plotEmbedding(ArchRProj = macaque_multiome_ArchR, colorBy = "cellColData", name = "Gex_cell_type", embedding = "UMAP")

p1+p2+p3+plot_layout(ncol = 3)

#save data
cell_list <- rownames(macaque_multiome_ArchR@cellColData)
saveRDS(cell_list,file = './res/step_6_fig_220227/macaque_multiome_ArchR_cell_list_220302_2253.rds')

# scRNA-seq filter --------------------------------------------------------
macaque_integration_Seurat <- readRDS(file = '../processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_annotated_220302.rds')
macaque_RNA_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'scRNA']
macaque_multiome_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']

#multiome ATAC filter
macaque_multiome_seurat <- macaque_multiome_seurat[,cell_list]

macaque_multiome_seurat <- my_process_seurat(object = macaque_multiome_seurat,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 3000,vars.to.regress = c('nCount_RNA','donor'),npcs = 50,preprocess = TRUE)
macaque_multiome_seurat <- my_process_seurat(object = macaque_multiome_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 32,resolution = c(0.5,1,1.5),group.by = 'old_cell_type',label = TRUE)

#harmony
DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

macaque_RNA_seurat <- FindVariableFeatures(object = macaque_RNA_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- VariableFeatures(macaque_RNA_seurat)
table(gene_list %in% rownames(macaque_multiome_seurat))

dim_used <- 30

macaque_integration_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_seurat),
                                                    assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                          multiome=c('nCount_RNA','donor')),
                                                    npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                    UMAP_dim = dim_used,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                                    yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_seurat),"cell_type"] <- macaque_multiome_seurat$old_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$old_cell_type

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 20 seems weird
DimPlot(macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cells.highlight = colnames(macaque_integration_Seurat)[macaque_integration_Seurat$seurat_clusters == '20'])
table(macaque_integration_Seurat[,macaque_integration_Seurat$seurat_clusters == '20']$cell_type)
VlnPlot(macaque_integration_Seurat,group.by = 'seurat_clusters',features = c('nCount_integration','nFeature_integration'),pt.size = 0,assay = 'integration',ncol = 1)
#cluster 20 seems low quality
table(macaque_integration_Seurat[,macaque_integration_Seurat$seurat_clusters == '20']$dataset)
#not tech specific

macaque_integration_Seurat <- macaque_integration_Seurat[,!(macaque_integration_Seurat$seurat_clusters %in% c('20'))]
p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 24 seems weird
DimPlot(macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cells.highlight = colnames(macaque_integration_Seurat)[macaque_integration_Seurat$seurat_clusters == '24'])

#save data
cell_list <- colnames(macaque_integration_Seurat)
saveRDS(cell_list,file = './res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220302_0047.rds')

#redo harmony
macaque_multiome_seurat <- macaque_multiome_seurat[,colnames(macaque_multiome_seurat) %in% cell_list]
macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_seurat) %in% cell_list]

DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

macaque_RNA_seurat <- FindVariableFeatures(object = macaque_RNA_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- VariableFeatures(macaque_RNA_seurat)
table(gene_list %in% rownames(macaque_multiome_seurat))

dim_used <- 26

macaque_integration_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_seurat),
                                                     assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                           multiome=c('nCount_RNA','donor')),
                                                     npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                     UMAP_dim = dim_used,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                                     yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_seurat),"cell_type"] <- macaque_multiome_seurat$old_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$old_cell_type

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 23 and 25 seems weird
DimPlot(macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cells.highlight = colnames(macaque_integration_Seurat)[macaque_integration_Seurat$seurat_clusters == '23'])
DimPlot(macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cells.highlight = colnames(macaque_integration_Seurat)[macaque_integration_Seurat$seurat_clusters == '25'])
table(macaque_integration_Seurat[,macaque_integration_Seurat$seurat_clusters == '23']$cell_type)
#filter cluster 23
macaque_integration_Seurat <- macaque_integration_Seurat[,!(macaque_integration_Seurat$seurat_clusters %in% c('23'))]

#save data
cell_list <- colnames(macaque_integration_Seurat)
saveRDS(cell_list,file = './res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220304_1548.rds')

#redo harmony
macaque_multiome_seurat <- macaque_multiome_seurat[,colnames(macaque_multiome_seurat) %in% cell_list]
macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_seurat) %in% cell_list]

DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

macaque_RNA_seurat <- FindVariableFeatures(object = macaque_RNA_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- VariableFeatures(macaque_RNA_seurat)
table(gene_list %in% rownames(macaque_multiome_seurat))

dim_used <- 25
#try 25
macaque_integration_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_seurat),
                                                     assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                           multiome=c('nCount_RNA','donor')),
                                                     npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                     UMAP_dim = dim_used,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                                     yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_seurat),"cell_type"] <- macaque_multiome_seurat$old_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$old_cell_type

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#cluster 23 seems weird
DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cells.highlight = colnames(macaque_integration_Seurat)[macaque_integration_Seurat$seurat_clusters == '23'])
table(macaque_integration_Seurat[,macaque_integration_Seurat$seurat_clusters == '23']$dataset)

VlnPlot(macaque_integration_Seurat,features = c('nCount_integration','nFeature_integration'),group.by = 'seurat_clusters',pt.size = 0)

macaque_integration_Seurat <- macaque_integration_Seurat[,!(macaque_integration_Seurat$seurat_clusters %in% c('23'))]
p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#save data
cell_list <- colnames(macaque_integration_Seurat)
saveRDS(cell_list,file = './res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220305_0844.rds')

#redo harmony
macaque_multiome_seurat <- macaque_multiome_seurat[,colnames(macaque_multiome_seurat) %in% cell_list]
macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_seurat) %in% cell_list]

DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

macaque_RNA_seurat <- FindVariableFeatures(object = macaque_RNA_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
macaque_multiome_seurat <- FindVariableFeatures(object = macaque_multiome_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- dplyr::intersect(VariableFeatures(macaque_RNA_seurat),VariableFeatures(macaque_multiome_seurat))
table(gene_list %in% rownames(macaque_multiome_seurat))

dim_used <- 25

macaque_integration_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_seurat),
                                                     assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                           multiome=c('nCount_RNA','donor')),
                                                     npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                     UMAP_dim = dim_used,resolution = c(0.5,1,2),kmeans_init_iter_max = 200,
                                                     yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_seurat),"cell_type"] <- macaque_multiome_seurat$old_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$old_cell_type

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#filter cluster 34
macaque_integration_Seurat <- macaque_integration_Seurat[,!(macaque_integration_Seurat$seurat_clusters %in% c('34'))]
p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

temp <- colnames(macaque_integration_Seurat)[macaque_integration_Seurat@reductions$umap@cell.embeddings[,1] > -5 & macaque_integration_Seurat@reductions$umap@cell.embeddings[,1] < 5 & macaque_integration_Seurat@reductions$umap@cell.embeddings[,2] > 1.5 & macaque_integration_Seurat@reductions$umap@cell.embeddings[,2] < 4]
DimPlot(object = macaque_integration_Seurat,cells.highlight = temp)

#filter weird cell
macaque_integration_Seurat <- macaque_integration_Seurat[,!(colnames(macaque_integration_Seurat) %in% temp)]
p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#save data
cell_list <- colnames(macaque_integration_Seurat)
saveRDS(cell_list,file = './res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220305_0939.rds')

#redo harmony
cell_list <- readRDS(file = './res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220305_0939.rds')
macaque_integration_Seurat <- readRDS(file = '../processed_data/220226_summary/macaque_RNA_multiome_harmony_integration_seurat_annotated_220302.rds')
macaque_RNA_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'scRNA']
macaque_multiome_seurat <- macaque_integration_Seurat[,macaque_integration_Seurat$tech == 'multiome']

macaque_multiome_seurat <- macaque_multiome_seurat[,colnames(macaque_multiome_seurat) %in% cell_list]
macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_seurat) %in% cell_list]

DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

macaque_RNA_seurat <- FindVariableFeatures(object = macaque_RNA_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- VariableFeatures(macaque_RNA_seurat)
table(gene_list %in% rownames(macaque_multiome_seurat))

dim_used <- 25

macaque_integration_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_seurat),
                                                     assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                           multiome=c('nCount_RNA','donor')),
                                                     npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                     UMAP_dim = dim_used,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                                     yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_seurat),"cell_type"] <- macaque_multiome_seurat$old_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$old_cell_type

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

temp <- colnames(macaque_integration_Seurat)[macaque_integration_Seurat$seurat_clusters %in% c('11') & macaque_integration_Seurat$cell_type %in% c('InCGE','InMGE')]
DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE,cells.highlight = temp)
DimPlot(object = macaque_integration_Seurat[,!(colnames(macaque_integration_Seurat) %in% temp)],group.by = 'cell_type',label = TRUE,repel = TRUE)

#save data
macaque_integration_Seurat <- macaque_integration_Seurat[,!(colnames(macaque_integration_Seurat) %in% temp)]
cell_list <- colnames(macaque_integration_Seurat)
saveRDS(cell_list,file = './res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220306_1630.rds')

#redo harmony
cell_list <- readRDS(file = './res/step_6_fig_220227/macaque_integration_Seurat_cell_list_220306_1630.rds')

macaque_multiome_seurat <- macaque_multiome_seurat[,colnames(macaque_multiome_seurat) %in% cell_list]
macaque_RNA_seurat <- macaque_RNA_seurat[,colnames(macaque_RNA_seurat) %in% cell_list]

DimPlot(macaque_RNA_seurat,group.by = 'old_cell_type',label = TRUE,repel = TRUE,reduction = 'harmonyUMAP')
macaque_multiome_seurat$tech <- 'multiome'
macaque_RNA_seurat$tech <- 'scRNA'

macaque_RNA_seurat <- FindVariableFeatures(object = macaque_RNA_seurat,assay = 'RNA',selection.method = 'vst',nfeatures = 3000,verbose = TRUE)
gene_list <- VariableFeatures(macaque_RNA_seurat)
table(gene_list %in% rownames(macaque_multiome_seurat))

dim_used <- 28
#try 28 31 32 
macaque_integration_Seurat <- my_harmony_integration(named_seurat_list = list(RNA=macaque_RNA_seurat,multiome=macaque_multiome_seurat),
                                                     assay = 'RNA',variable_feature = gene_list,var_to_regress_list = list(RNA=c('nCount_RNA','donor','batch'),
                                                                                                                           multiome=c('nCount_RNA','donor')),
                                                     npcs = 50,reference_loading = 'RNA',integration_var = 'tech',harmony_input_dim = dim_used,max.iter.harmony = 50,
                                                     UMAP_dim = dim_used,resolution = c(0.5,1,1.5),kmeans_init_iter_max = 200,
                                                     yiming_harmony = '/data/User/sunym/software/R_lib/yiming_harmony_R_4.1.2/')

my_send_sms('harmony done!')

macaque_integration_Seurat$cell_type <- NA
macaque_integration_Seurat@meta.data[colnames(macaque_multiome_seurat),"cell_type"] <- macaque_multiome_seurat$old_cell_type
macaque_integration_Seurat@meta.data[colnames(macaque_RNA_seurat),"cell_type"] <- macaque_RNA_seurat$old_cell_type

p1 <- DimPlot(object = macaque_integration_Seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_integration_Seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
p3 <- DimPlot(object = macaque_integration_Seurat,group.by = 'dataset',label = FALSE)
p1+p2+p3+plot_layout(ncol = 3)

#prepare to annotate
macaque_integration_Seurat <- NormalizeData(object = macaque_integration_Seurat)
VlnPlot(macaque_integration_Seurat,features = c('EGFR'),group.by = 'seurat_clusters',pt.size = 0,slot = 'data')
FeaturePlot(macaque_integration_Seurat,features = c('EGFR','EOMES','PPP1R17'))
