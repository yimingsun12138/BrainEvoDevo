#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: macaque progenitor cell atlas                                   ##
## Data: 2021.11.22                                                                ##
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
library(Rmisc)
library(Seurat)
library(dplyr)
library(patchwork)
library(scibet)
library(ggsci)
library(ggplot2)
library(cowplot)
library(ComplexHeatmap)
library(harmony)
library(monocle)
library(monocle3)
library(SingleCellExperiment)
library(clusterProfiler)
library(org.Mmu.eg.db)
library(topGO)
library(destiny)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('/data/User/sunym/back_up/twilio_send_messages.R')

#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,group.by = 'cell_type',label = TRUE,repel = TRUE)
DimPlot(macaque_RNA_seurat,group.by = 'seurat_clusters',label = TRUE,repel = TRUE)

# cell cycle related gene -------------------------------------------------
cell_cycle_gene <- getGOgeneSet(x = 'GO:0007049',OrgDb = 'org.Mmu.eg.db',ont = 'BP',keytype = 'SYMBOL')
table(cell_cycle_gene %in% VariableFeatures(macaque_RNA_seurat))

#macaque_modified_seurat umap
macaque_modified_seurat <- macaque_RNA_seurat@assays$RNA@counts
macaque_modified_seurat <- macaque_modified_seurat[!(rownames(macaque_modified_seurat) %in% cell_cycle_gene),]
dim(macaque_RNA_seurat)
dim(macaque_modified_seurat)
macaque_modified_seurat <- CreateSeuratObject(counts = macaque_modified_seurat,project = 'macaque',assay = 'RNA',meta.data = macaque_RNA_seurat@meta.data[,4:11],min.cells = 0,min.features = 0)

gene_list <- dplyr::intersect(x = rownames(macaque_modified_seurat),y = VariableFeatures(macaque_RNA_seurat))
macaque_modified_seurat <- NormalizeData(object = macaque_modified_seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)
VariableFeatures(macaque_modified_seurat) <- gene_list
macaque_modified_seurat <- ScaleData(object = macaque_modified_seurat,features = VariableFeatures(macaque_modified_seurat),assay = 'RNA',vars.to.regress = c("donor","batch","nCount_RNA"))

#projection
temp <- macaque_modified_seurat@assays$RNA@scale.data
temp <- projectMatrix_SeuratUMAP(X_scaled = temp,object = macaque_RNA_seurat,umap_model = NULL,assayUsed = 'RNA',missing_gene = TRUE)

macaque_RNA_seurat@reductions$UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_orig,assay = 'RNA')
macaque_modified_seurat@reductions$UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')

p1 <- DimPlot(object = macaque_RNA_seurat,group.by = 'cell_type',reduction = 'UMAP',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_modified_seurat,group.by = 'cell_type',reduction = 'UMAP',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)
#seems no obvious change...

#Cycling cell markers
gene_list <- FindMarkers(object = macaque_RNA_seurat,ident.1 = c('Cyc-S','Cyc-G2M'),group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE)
VlnPlot(object = macaque_RNA_seurat,features = rownames(gene_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')
gene_list <- rownames(gene_list[gene_list$pct.2 < 0.2,])

#topGO analysis
temp <- rownames(macaque_RNA_seurat)
temp <- c(temp %in% gene_list)
temp[temp == TRUE] <- 1
temp[temp == FALSE] <- 0
table(temp)
names(temp) <- rownames(macaque_RNA_seurat)
gene_list <- temp
gene_list <- factor(gene_list,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "Cycling enriched",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Mmu.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))

#gene gene list
gene_list <- allRes$GO.ID
gene_list <- base::lapply(X = gene_list,FUN = function(x){
  return(getGOgeneSet(x = x,OrgDb = 'org.Mmu.eg.db',ont = 'BP',keytype = 'SYMBOL'))
})
gene_list <- unique(unlist(gene_list))

#reprogection
macaque_modified_seurat <- macaque_RNA_seurat@assays$RNA@counts
macaque_modified_seurat <- CreateSeuratObject(counts = macaque_modified_seurat,project = 'macaque',assay = 'RNA',meta.data = macaque_RNA_seurat@meta.data[,4:11],min.cells = 0,min.features = 0)
gene_list <- VariableFeatures(macaque_RNA_seurat)[!(VariableFeatures(macaque_RNA_seurat) %in% gene_list)]
macaque_modified_seurat <- NormalizeData(object = macaque_modified_seurat,assay = 'RNA',normalization.method = 'LogNormalize',scale.factor = 10000)
VariableFeatures(macaque_modified_seurat) <- gene_list
macaque_modified_seurat <- ScaleData(object = macaque_modified_seurat,features = VariableFeatures(macaque_modified_seurat),assay = 'RNA',vars.to.regress = c("donor","batch","nCount_RNA"))

temp <- macaque_modified_seurat@assays$RNA@scale.data
temp <- projectMatrix_SeuratUMAP(X_scaled = temp,object = macaque_RNA_seurat,umap_model = NULL,assayUsed = 'RNA',missing_gene = TRUE)

macaque_RNA_seurat@reductions$UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_orig,assay = 'RNA')
macaque_modified_seurat@reductions$UMAP <- CreateDimReducObject(embeddings = temp$umapCoord_proj,assay = 'RNA')

p1 <- DimPlot(object = macaque_RNA_seurat,group.by = 'cell_type',reduction = 'UMAP',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = macaque_modified_seurat,group.by = 'cell_type',reduction = 'UMAP',label = TRUE,repel = TRUE)
p1+p2+plot_layout(ncol = 2)

#still no difference ... mother fuck!


# determing astrocyte -----------------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
DimPlot(macaque_RNA_seurat,reduction = 'umap',group.by = 'cell_type',label = TRUE,repel = TRUE)

#astrocyte marker 
gene_list <- FindMarkers(object = macaque_RNA_seurat,ident.1 = 'Astrocyte',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE)
pdf(file = './res/step_14_fig_211122/Astrocyte_marker_gene_vlnplot.pdf',width = 18,height = 8)
VlnPlot(object = macaque_RNA_seurat,features = rownames(gene_list)[1:12],pt.size = 0,assay = 'RNA',group.by = 'cell_type',slot = 'data')
dev.off()
rownames(gene_list)[1:10]

#it is ependymal cell!

dotplot_matrix <- my_dotplot(macaque_RNA_seurat,assay = 'RNA', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             Mic=c('CX3CR1'),
                                             Ependymal=c('DNAH11','GJA1','SPATA17','LRRC71','SPAG17'),
                                             RG=c('SOX9','PAX6','VIM','FAM107A','HOPX','MOXD1','FBXO32','CRYAB','NR4A1','FOXJ1','NPY','FGFR3','CD9','GPX3'),
                                             OPC=c('SOX10','OLIG2','EGFR'),
                                             Cyc=c('TOP2A','MKI67','CLSPN','AURKA'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NEUROD6','TBR1','SATB2','SLC17A7','FEZF2'),
                                             In=c('DLX5','GAD2','GAD1','DLX2'),
                                             MGE=c('LHX6','SST'),
                                             CGE=c('SP8','NR2F2'),
                                             PSB=c('MEIS2','ETV1')),
                             group.by = 'cell_type', cols = c('#2CA02CFF','white','#D62728FF'),
                             return_data_plot = TRUE)

dotplot_matrix$id <- factor(dotplot_matrix$id,levels = c('InPSB','InCGE','InMGE','Ex-4','Ex-3','Ex-2','Ex-1','IP','Cyc-G2M','Cyc-S','OPC','RG-3','RG-2','RG-1','Astrocyte','Mic','Per','End'))

pdf(file = './res/step_14_fig_211122/macaque_RNA_seurat_ependymal_dotplot.pdf',width = 16,height = 8)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',size = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom',
        panel.grid = element_line(color="grey",size = 0.1)) + 
  xlab('') + ylab('')
dev.off()

#load public data
express_matrix <- read.table(file = './data/public/Single_cell_atlas_of_early_human_brain_development_highlights_heterogeneity_of_human_neuroepithelial_cells_and_early_radial_glia/express_matrix_seurat_normalized.tsv')
express_matrix <- expm1(express_matrix)
meta_data <- read.table(file = './data/public/Single_cell_atlas_of_early_human_brain_development_highlights_heterogeneity_of_human_neuroepithelial_cells_and_early_radial_glia/meta.tsv',sep = '\t',header = TRUE,row.names = 1)
table(rownames(meta_data) == colnames(express_matrix))
human_Brain_seurat <- CreateSeuratObject(counts = express_matrix,project = 'human',assay = 'RNA',meta.data = meta_data,min.cells = 0,min.features = 0)
rm(express_matrix)

table(human_Brain_seurat$Cell.Type)
human_Brain_seurat <- human_Brain_seurat[,human_Brain_seurat$Area %in% c('Frontal cortex')]
human_Brain_seurat <- my_process_seurat(object = human_Brain_seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = c('nCount_RNA','Individual'),npcs = 50,preprocess = TRUE)
human_Brain_seurat <- my_process_seurat(object = human_Brain_seurat,assay = 'RNA',reduction.name = 'pca',preprocess = FALSE,dim_to_use = 10,resolution = 1,group.by = 'Cell.Type',label = TRUE)

NE_marker <- FindMarkers(object = human_Brain_seurat,ident.1 = 'Neuroepithelial',group.by = 'Cell.Type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE)
VlnPlot(macaque_RNA_seurat,features = rownames(NE_marker)[49:60],assay = 'RNA',group.by = 'cell_type',slot = 'data',pt.size = 0)
gene_list <- dplyr::intersect(rownames(gene_list),rownames(NE_marker))
VlnPlot(macaque_RNA_seurat,features = gene_list[61:78],assay = 'RNA',group.by = 'cell_type',slot = 'data',pt.size = 0)

temp <- c('ID4','CCDC181','LRRIQ1','AHI1','CNTLN','ANKRD26','RBM23','RABL2B','TRMT13','WDR6','FZD3','GNPAT','CC2D2A','KDM3A','CEP290','OBSL1')
table(temp %in% gene_list)
VlnPlot(macaque_RNA_seurat,features = temp,group.by = 'cell_type',pt.size = 0,assay = 'RNA',slot = 'data')
VlnPlot(human_Brain_seurat[,human_Brain_seurat$Cell.Type != 'IPC'],features = rownames(NE_marker)[1:12],group.by = 'Cell.Type',pt.size = 0,assay = 'RNA',slot = 'data')
#dead end! poor public data!

#calculate correlation
temp <- AverageExpression(object = human_Brain_seurat,assays = 'RNA',return.seurat = FALSE,group.by = 'Cell.Type',slot = 'data')
temp <- temp$RNA

express_matrix <- macaque_RNA_seurat@assays$RNA@counts
macaque_to_human_anno <- read.csv(file = './data/reference/Mmul_10_to_GRCh38.csv')
express_matrix <- My_Convert_Homology_Gene_ID(express_matrix = express_matrix,anno = macaque_to_human_anno,filter_anno = TRUE,workers = 6,future.globals.maxSize = 20^(1024)^3)
express_matrix <- CreateSeuratObject(counts = express_matrix,project = 'temp',assay = 'RNA',meta.data = macaque_RNA_seurat@meta.data,min.cells = 0,min.features = 0)
express_matrix <- express_matrix[,express_matrix$cell_type %in% c('Astrocyte','Cyc-G2M','Cyc-S','Ex-1','IP','RG-1','RG-2')]
express_matrix <- AverageExpression(object = express_matrix,assays = 'RNA',return.seurat = FALSE,group.by = 'cell_type',slot = 'counts')
express_matrix <- express_matrix$RNA

cor_matrix <- matrix(nrow = dim(express_matrix)[2],ncol = dim(temp)[2],data = NA)
rownames(cor_matrix) <- colnames(express_matrix)
colnames(cor_matrix) <- colnames(temp)
cor_matrix <- as.data.frame(cor_matrix)

gene_list <- dplyr::union(VariableFeatures(macaque_RNA_seurat),VariableFeatures(human_Brain_seurat))
gene_list <- gene_list[gene_list %in% rownames(temp) & gene_list %in% rownames(express_matrix)]
temp <- temp[gene_list,]
express_matrix <- express_matrix[gene_list,]

for(i in rownames(cor_matrix)){
  for (j in colnames(cor_matrix)) {
    cor_matrix[i,j] <- cor(express_matrix[,i],temp[,j],method = 'pearson')
  }
}

#dead end too! fuck
rm(list = c('temp','gene_list','express_matrix','cor_matrix'))
gc()

temp <- c('ID4','CCDC181','LRRIQ1','AHI1','CNTLN','ANKRD26','RBM23','RABL2B','TRMT13','WDR6','FZD3','GNPAT','CC2D2A','KDM3A','CEP290','OBSL1')
pdf(file = './res/step_14_fig_211122/Astrocyte_NE_cell_marker_vlnplot.pdf',width = 18,height = 8)
VlnPlot(macaque_RNA_seurat,features = temp,group.by = 'cell_type',pt.size = 0,assay = 'RNA',slot = 'data')
dev.off()

#go analysis for shared gene of astrocyte and NE
gene_list <- FindMarkers(object = macaque_RNA_seurat,ident.1 = 'Astrocyte',group.by = 'cell_type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE)
NE_marker <- FindMarkers(object = human_Brain_seurat,ident.1 = 'Neuroepithelial',group.by = 'Cell.Type',assay = 'RNA',slot = 'data',test.use = 'bimod',only.pos = TRUE)
gene_list <- dplyr::intersect(rownames(gene_list),rownames(NE_marker))

#topGO analysis
temp <- rownames(macaque_RNA_seurat)
temp <- c(temp %in% gene_list)
temp[temp == TRUE] <- 1
temp[temp == FALSE] <- 0
table(temp)
names(temp) <- rownames(macaque_RNA_seurat)
gene_list <- temp
gene_list <- factor(gene_list,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "Cycling enriched",
                 ontology = "BP",
                 allGenes = gene_list,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Mmu.eg.db',ID = 'symbol')

resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)

allRes$sig <- -log10(as.numeric(allRes$classicFisher))
#also dead end

# progenitor development diffusion map ------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')

#diffusion map without cycling cell
macaque_modified_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('RG-1','RG-2','RG-3','IP','OPC','Astrocyte')]
macaque_modified_seurat <- my_process_seurat(object = macaque_modified_seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
express_matrix <- macaque_modified_seurat@assays$RNA@data[VariableFeatures(macaque_modified_seurat),]
express_matrix <- as.matrix(t(express_matrix))
sigma <- find_sigmas(express_matrix,verbose = TRUE)
optimal_sigma(sigma)
dm <- DiffusionMap(data = express_matrix,sigma = 27,k = 30,n_eigs = 30,density_norm = TRUE,distance = 'cosine')

plot(dm)
temp <- data.frame(DC_1 = dm$DC1,DC_2 = dm$DC2,DC_3 = dm$DC3,DC_4 = dm$DC4,DC_5 = dm$DC5)
temp$cell_type <- macaque_modified_seurat@meta.data[rownames(temp),"cell_type"]

pdf(file = './res/step_14_fig_211122/Astrocyte_no_cycling_diffusion_map.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x=DC_1,y=DC_2)) + 
  geom_point(aes(color = cell_type),size = 0.5,alpha = 0.5) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = c(0.8,0.75),
        plot.title = element_text(size = 12,face = 'bold',hjust = 0.5),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1,fill = NA,colour = 'Black')) + 
  labs(title = 'Diffusion map')
dev.off()

#diffusion map with cycling cells
macaque_modified_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('RG-1','RG-2','RG-3','IP','OPC','Astrocyte','Cyc-G2M','Cyc-S')]
macaque_modified_seurat <- my_process_seurat(object = macaque_modified_seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
express_matrix <- macaque_modified_seurat@assays$RNA@data[VariableFeatures(macaque_modified_seurat),]
express_matrix <- as.matrix(t(express_matrix))
sigma <- find_sigmas(express_matrix,verbose = TRUE)
optimal_sigma(sigma)
dm <- DiffusionMap(data = express_matrix,sigma = 23,k = 30,n_eigs = 30,density_norm = TRUE,distance = 'cosine')

temp <- data.frame(DC_1 = dm$DC1,DC_2 = dm$DC2,DC_3 = dm$DC3,DC_4 = dm$DC4,DC_5 = dm$DC5)
temp$cell_type <- macaque_modified_seurat@meta.data[rownames(temp),"cell_type"]

ggplot(data = temp,aes(x=DC_1,y=DC_2)) + 
  geom_point(aes(color = cell_type),size = 1,alpha = 1) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = c(0.8,0.75),
        plot.title = element_text(size = 12,face = 'bold',hjust = 0.5),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1,fill = NA,colour = 'Black')) + 
  labs(title = 'Diffusion map')

#diffusion map without cycling cells and Astrocyte
macaque_modified_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('RG-1','RG-2','RG-3','IP','OPC')]
macaque_modified_seurat <- my_process_seurat(object = macaque_modified_seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
express_matrix <- macaque_modified_seurat@assays$RNA@data[VariableFeatures(macaque_modified_seurat),]
express_matrix <- as.matrix(t(express_matrix))
sigma <- find_sigmas(express_matrix,verbose = TRUE)
optimal_sigma(sigma)
dm <- DiffusionMap(data = express_matrix,sigma = 'local',k = 30,n_eigs = 30,density_norm = TRUE,distance = 'cosine')

temp <- data.frame(DC_1 = dm$DC1,DC_2 = dm$DC2,DC_3 = dm$DC3,DC_4 = dm$DC4,DC_5 = dm$DC5)
temp$cell_type <- macaque_modified_seurat@meta.data[rownames(temp),"cell_type"]

pdf(file = './res/step_14_fig_211122/no_Astrocyte_no_Cycling_diffusion_map.pdf',width = 5,height = 5)
ggplot(data = temp,aes(x=DC_1,y=DC_2)) + 
  geom_point(aes(color = cell_type),size = 1,alpha = 1) + 
  theme_classic() + 
  theme(aspect.ratio = 1,
        legend.position = c(0.8,0.75),
        plot.title = element_text(size = 12,face = 'bold',hjust = 0.5),
        axis.line = element_blank(),
        panel.background = element_rect(size = 1,fill = NA,colour = 'Black')) + 
  labs(title = 'Diffusion map')
dev.off()

# progenitor development URD ----------------------------------------------
#load data
macaque_RNA_seurat <- readRDS(file = './processed_data/211014_summary/macaque_200919_210922_merged_RNA_seurat_211021.rds')
macaque_RNA_seurat <- macaque_RNA_seurat[,macaque_RNA_seurat$cell_type %in% c('RG-1','RG-2','RG-3','IP','OPC','Astrocyte')]
macaque_RNA_URD <- createURD(count.data = macaque_RNA_seurat@assays$RNA@counts,meta = macaque_RNA_seurat@meta.data,min.cells = 0,min.genes = 0,min.counts = 0,gene.max.cut = 7000)
macaque_RNA_URD@group.ids$cell_type <- as.character(macaque_RNA_URD@meta$cell_type)
head(macaque_RNA_URD@group.ids)

#variable genes
macaque_RNA_seurat <- my_process_seurat(object = macaque_RNA_seurat,assay = 'RNA',reduction.name = 'pca',nfeatures = 2000,vars.to.regress = c('nCount_RNA','donor','batch'),npcs = 50,preprocess = TRUE)
macaque_RNA_URD@var.genes <- VariableFeatures(macaque_RNA_seurat)

#calculate diffusion map
macaque_RNA_URD <- calcDM(object = macaque_RNA_URD,knn = 30,sigma.use = 27,distance = 'cosine',verbose = TRUE)
plotDimArray(macaque_RNA_URD, reduction.use = "dm", dims.to.plot = 1:8, outer.title = "Diffusion Map", label="cell_type", plot.title="", legend=F)

#calculate pseudotime
root.cells <- cellsInCluster(macaque_RNA_URD, "cell_type", "Astrocyte")
macaque_RNA_URD.floods <- floodPseudotime(macaque_RNA_URD, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=TRUE)
macaque_RNA_URD <- floodPseudotimeProcess(macaque_RNA_URD, macaque_RNA_URD.floods, floods.name="pseudotime")
pseudotimePlotStabilityOverall(macaque_RNA_URD)
plotDim(object = macaque_RNA_URD,label = 'pseudotime',reduction.use = 'dm')

#find tips
cell_list <- rownames(macaque_RNA_URD@meta)[macaque_RNA_URD@dm$DC1>0.02 & macaque_RNA_URD@meta$cell_type == 'OPC']
cell_list <- append(cell_list,rownames(macaque_RNA_URD@meta)[macaque_RNA_URD@dm$DC1< -0.025 & macaque_RNA_URD@meta$cell_type == 'IP'])
macaque_RNA_URD@group.ids$tip.clusters <- NA
macaque_RNA_URD@group.ids[cell_list,"tip.clusters"] <- macaque_RNA_URD@meta[cell_list,"cell_type"]
table(macaque_RNA_URD@group.ids$tip.clusters)
macaque_RNA_URD@group.ids$tip.clusters[macaque_RNA_URD@group.ids$tip.clusters == 'OPC'] <- '1'
macaque_RNA_URD@group.ids$tip.clusters[macaque_RNA_URD@group.ids$tip.clusters == 'IP'] <- '2'
table(macaque_RNA_URD@group.ids$tip.clusters)

plotDim(object = macaque_RNA_URD,label = 'tip.clusters',reduction.use = 'dm')

#Biased random walks
macaque_RNA_URD.ptlogistic <- pseudotimeDetermineLogistic(macaque_RNA_URD, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = TRUE)
macaque_RNA_URD.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(macaque_RNA_URD, "pseudotime", logistic.params=macaque_RNA_URD.ptlogistic))
macaque_RNA_URD.walks <- simulateRandomWalksFromTips(macaque_RNA_URD, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = macaque_RNA_URD.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 20000,verbose = TRUE)
macaque_RNA_URD <- processRandomWalksFromTips(macaque_RNA_URD, macaque_RNA_URD.walks, verbose = TRUE)

plotDim(macaque_RNA_URD, "visitfreq.log.1", plot.title="Visitation frequency from tip 1 (log10)", transitions.plot=10000,reduction.use = 'dm')
plotDim(macaque_RNA_URD, "visitfreq.log.2", plot.title="Visitation frequency from tip 2 (log10)", transitions.plot=10000,reduction.use = 'dm')

#Build tree
macaque_RNA_URD.tree <- loadTipCells(macaque_RNA_URD, "tip.clusters")
macaque_RNA_URD.tree <- buildTree(macaque_RNA_URD.tree, pseudotime = "pseudotime", tips.use=1:2, divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = TRUE, p.thresh=0.001)

plotTree(macaque_RNA_URD.tree, "cell_type", title="Developmental Stage")

#not a good choice

# progenitor development monocle2 -----------------------------------------
