#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: genes affected by human SNC                                     ##
## Data: 2022.11.12                                                                ##
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

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')
source('/content/script/twilio_send_messages.R')

#initialize ArchR
addArchRThreads(threads = 5)

# try IP human specific human SNC -----------------------------------------

#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')
DF_list <- readRDS(file = './res/step_78_fig_221111/DF_list/IP_DF_list.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#get EP link
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#explore p2g
p2g@metadata$peakSet
p2g@metadata$geneSet
p2g

#create correspondence
human_peakset <- p2g@metadata$peakSet
names(human_peakset) <- paste(human_peakset@seqnames,as.character(human_peakset@ranges),sep = '-')
idxATAC <- names(human_peakset)

human_geneset <- p2g@metadata$geneSet
idxRNA <- as.character(human_geneset$name)

p2g_data_frame <- data.frame(idxATAC = idxATAC[p2g$idxATAC],idxRNA = idxRNA[p2g$idxRNA])

#get the peak that overlap with human specific SNC
SNC_list <- DF_list[DF_list$group == 'human_specific',"human_SNC"]
names(human_SNC) <- human_SNC$human_SNC
SNC_list <- human_SNC[SNC_list]

peak_list <- countOverlaps(query = human_peakset,subject = SNC_list)
peak_list <- names(peak_list)[peak_list > 0]

gene_list <- p2g_data_frame[p2g_data_frame$idxATAC %in% peak_list,"idxRNA"]

#try GO
human_all_gene <- unique(human_geneset$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
pdf(file = './res/step_79_fig_221112/human_specific_SNC_in_IP_GO_BP_barplot.pdf',width = 5,height = 3)
ggplot(allRes[c(1:3,5:6),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab('human specific SNC in IP')
dev.off()

#save gene list 
gene_list <- unique(gene_list)
saveRDS(object = gene_list,file = './res/step_79_fig_221112/gene_regulated_by_IP_human_specific_SNC.rds')

## try to find the expression difference -----------------------------------

#try to use gene quantile

#macaque
temp <- macaque_multiome_Seurat@assays$RNA@counts
temp <- as.matrix(temp)
options(future.globals.maxSize = 200*(1024^3))
plan(multisession,workers = 5)
temp <- base::do.call(what = cbind,args = future.apply::future_lapply(X = colnames(temp),FUN = function(x){
  x <- temp[,x]
  x <- rank(x)/length(x)
  return(x)
}))
plan(sequential)

colnames(temp) <- colnames(macaque_multiome_Seurat@assays$RNA@counts)
rownames(temp) <- rownames(macaque_multiome_Seurat@assays$RNA@counts)
temp <- as(temp,'sparseMatrix')
gc()

saveRDS(object = temp,file = './res/step_79_fig_221112/macque_RNA_quantile_matrix.rds')

#human
temp <- Greenleaf_RNA_Seurat@assays$RNA@counts
temp <- as.matrix(temp)
options(future.globals.maxSize = 200*(1024^3))
plan(multisession,workers = 10)
temp <- base::do.call(what = cbind,args = future.apply::future_lapply(X = colnames(temp),FUN = function(x){
  x <- temp[,x]
  x <- rank(x)/length(x)
  return(x)
}))
plan(sequential)

colnames(temp) <- colnames(Greenleaf_RNA_Seurat@assays$RNA@counts)
rownames(temp) <- rownames(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp <- as(temp,'sparseMatrix')
gc()

saveRDS(object = temp,file = './res/step_79_fig_221112/Greenleaf_RNA_quantile_matrix.rds')

#add matrix to seurat object
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

#try featureplot
DefaultAssay(macaque_multiome_Seurat) <- 'quantile'
FeaturePlot(object = macaque_multiome_Seurat,features = c('GAD2'),pt.size = 0.1,slot = 'data')
VlnPlot(object = macaque_multiome_Seurat,features = c('SLC1A3','SOX9','PAX6','PTPRZ1','EGFR'),group.by = 'cell_type',pt.size = 0,assay = 'quantile',slot = 'data')

my_dotplot <- function (object, assay = NULL, features, cols = c("lightgrey","blue"), 
                        col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, scale = FALSE, 
                        group.by = NULL, scale.by = "radius", scale.min = NA, scale.max = NA, 
                        return_data_plot = FALSE, data_plot = NULL) {
  require(Seurat)
  #default setting modified from Seurat::Dotplot
  idents <- NULL
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, 
                       stop("'scale.by' must be either 'size' or 'radius'"))
  #all settled
  
  if(is.null(data_plot)){
    PercentAbove <- function(x, threshold) {
      return(length(x = x[x > threshold]) / length(x = x))
    }
    assay <- assay %||% DefaultAssay(object = object)
    DefaultAssay(object = object) <- assay
    feature.groups <- NULL
    if (is.list(features) | any(!is.na(names(features)))) {
      feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                          FUN = function(x) {
                                            return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                          }))
      if (any(is.na(x = feature.groups))) {
        warning("Some feature groups are unnamed.", call. = FALSE, 
                immediate. = TRUE)
      }
      features <- unlist(x = features)
      names(x = feature.groups) <- features
    }
    cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
    data.features <- FetchData(object = object, vars = features, cells = cells)
    data.features$id <- if (is.null(x = group.by)) {
      Idents(object = object)[cells, drop = TRUE]
    }
    else {
      object[[group.by, drop = TRUE]][cells, drop = TRUE]
    }
    if (!is.factor(x = data.features$id)) {
      data.features$id <- factor(x = data.features$id)
    }
    id.levels <- levels(x = data.features$id)
    data.features$id <- as.vector(x = data.features$id)
    
    data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
      data.use <- data.features[data.features$id == ident, 1:(ncol(x = data.features) - 1), drop = FALSE]
      avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
        return(mean(x = expm1(x = x)))
      })
      pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, threshold = 0.6)
      return(list(avg.exp = avg.exp, pct.exp = pct.exp))
    })
    
    names(x = data.plot) <- unique(x = data.features$id)
    
    data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
      data.use <- as.data.frame(x = data.plot[[x]])
      data.use$features.plot <- rownames(x = data.use)
      data.use$id <- x
      return(data.use)
    })
    data.plot <- do.call(what = "rbind", args = data.plot)
    if (!is.null(x = id.levels)) {
      data.plot$id <- factor(x = data.plot$id, levels = id.levels)
    }
    if (length(x = levels(x = data.plot$id)) == 1) {
      scale <- FALSE
      warning("Only one identity present, the expression values will be not scaled",
              call. = FALSE, immediate. = TRUE)
    }
    avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot),
                             FUN = function(x) {
                               data.use <- data.plot[data.plot$features.plot == x, "avg.exp"]
                               if (scale) {
                                 data.use <- scale(x = data.use,center = TRUE)
                                 data.use <- MinMax(data = data.use, min = col.min, max = col.max)
                               }
                               else {
                                 data.use <- log1p(x = data.use)
                               }
                               return(data.use)
                             })
    avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
    
    data.plot$avg.exp.scaled <- avg.exp.scaled
    data.plot$features.plot <- factor(x = data.plot$features.plot, levels = features)
    data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
    data.plot$pct.exp <- data.plot$pct.exp * 100
    
    color.by <- "avg.exp.scaled"
    if (!is.na(x = scale.min)) {
      data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
    }
    if (!is.na(x = scale.max)) {
      data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
    }
    if (!is.null(x = feature.groups)) {
      data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                         levels = unique(x = feature.groups))
    }
    if(return_data_plot){
      return(data.plot)
    } else{
      plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
        geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
        scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
        theme(axis.title.x = element_blank(), 
              axis.title.y = element_blank()) + 
        guides(size = guide_legend(title = "Percent Expressed")) + 
        labs(x = "Features", y = "Identity") + 
        theme_cowplot()
      if (!is.null(x = feature.groups)) {
        plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                                  space = "free_x", switch = "y") + 
          theme(panel.spacing = unit(x = 1, units = "lines"), strip.background = element_blank())
      }
      
      plot <- plot + scale_colour_gradientn(colours = cols)
      plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
      return(plot)
    }
  } else{
    data.plot <- data_plot
    color.by <- "avg.exp.scaled"
    plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", y = "id")) + 
      geom_point(mapping = aes_string(size = "pct.exp", color = color.by)) + 
      scale.func(range = c(0, dot.scale), limits = c(scale.min, scale.max)) + 
      theme(axis.title.x = element_blank(), 
            axis.title.y = element_blank()) + 
      guides(size = guide_legend(title = "Percent Expressed")) + 
      labs(x = "Features", y = "Identity") + 
      theme_cowplot()
    if (!is.null(x = data.plot$feature.groups)) {
      plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                                space = "free_x", switch = "y") + 
        theme(panel.spacing = unit(x = 1, units = "lines"), strip.background = element_blank())
    }
    
    plot <- plot + scale_colour_gradientn(colours = cols)
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
    return(plot)
  }
}

dotplot_matrix <- my_dotplot(macaque_multiome_Seurat,assay = 'quantile', 
                             col.max = 2.5, col.min = -2.5, scale = TRUE, 
                             features = list(RG=c('SLC1A3','SOX9','PAX6','PTPRZ1','EGFR'),
                                             Cyc=c('TOP2A','MKI67'),
                                             IP=c('EOMES','PPP1R17'),
                                             Ex=c('NEUROD2','NRP1','SATB2','STMN2','NEUROD6','LMO4'),
                                             SP=c('NR4A2','CRYM','CDH18'),
                                             In=c('GAD1','DLX5','ADARB2','LHX6'),
                                             OPC=c('SOX10','OLIG2'),
                                             End=c('CLDN5','PECAM1'),
                                             Per=c('PDGFRB'),
                                             VLMC = c('COL1A1','LUM'),
                                             Mic=c('CX3CR1')),
                             group.by = 'cell_type', cols = c('#3B4992FF','white','#EE0000FF'),
                             return_data_plot = TRUE)
dotplot_matrix$id <- factor(dotplot_matrix$id,levels = names(color_param$celltype))

pdf(file = './res/step_79_fig_221112/macaque_RNA_quantile_marker_dotplot.pdf',width = 13,height = 6)
my_dotplot(data_plot = dotplot_matrix,col.max = 2.5, col.min = -2.5,
           cols = c('#3B4992FF','white','#EE0000FF'),return_data_plot = FALSE) + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        panel.border = element_rect(fill = NA,colour = 'black',linewidth = 0.5),
        axis.line = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = 'grey',colour = 'black'),
        legend.position = 'bottom') + 
  xlab('') + ylab('')
dev.off()

#consist with RNA express

#macaque
temp_quantile <- rowSums(macaque_multiome_Seurat@assays$quantile@data)
temp_quantile <- temp_quantile/sum(temp_quantile)*10^6
temp_RNA <- rowSums(macaque_multiome_Seurat@assays$RNA@counts)
temp_RNA <- temp_RNA/sum(temp_RNA)*10^6
temp <- data.frame(RNA = temp_RNA,quantile = temp_quantile)

pdf(file = './res/step_79_fig_221112/macaque_RNA_quantile_vs_counts_dotplot.pdf',width = 3.5,height = 4)
ggplot(data = temp,aes(x = RNA,y = quantile)) + 
  geom_point(size = 0.1) + 
  xlim(c(0,500)) + ylim(c(26,47)) + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(face = 'bold',size = 12,hjust = 0.5)) + 
  labs(title = 'macaque RNA quantile vs counts')
dev.off()

#human
temp_quantile <- rowSums(Greenleaf_RNA_Seurat@assays$quantile@data)
temp_quantile <- temp_quantile/sum(temp_quantile)*10^6
temp_RNA <- rowSums(Greenleaf_RNA_Seurat@assays$RNA@counts)
temp_RNA <- temp_RNA/sum(temp_RNA)*10^6
temp <- data.frame(RNA = temp_RNA,quantile = temp_quantile)

pdf(file = './res/step_79_fig_221112/Greenleaf_RNA_quantile_vs_counts_dotplot.pdf',width = 3.5,height = 4)
ggplot(data = temp,aes(x = RNA,y = quantile)) + 
  geom_point(size = 0.1) + 
  xlim(c(0,500)) + ylim(c(27,55)) + 
  theme_bw() + 
  theme(aspect.ratio = 1,
        plot.title = element_text(face = 'bold',size = 12,hjust = 0.5)) + 
  labs(title = 'Greenleaf RNA quantile vs counts')
dev.off()

#expression quantile seems consist with RNA expression

#gene expression quantile change between human and macaque
gene_list <- unique(gene_list)
gene_list <- gene_list[gene_list %in% rownames(macaque_multiome_Seurat@assays$quantile@data)]

temp_human <- colnames(Greenleaf_RNA_Seurat@assays$quantile@data)[Greenleaf_RNA_Seurat$ReAnno_celltype == 'IP']
temp_macaque <- colnames(macaque_multiome_Seurat@assays$quantile@data)[macaque_multiome_Seurat$cell_type == 'IP']

quantile_matrix <- cbind(Greenleaf_RNA_Seurat@assays$quantile@data[gene_list,temp_human],
                         macaque_multiome_Seurat@assays$quantile@data[gene_list,temp_macaque])
quantile_matrix <- CreateSeuratObject(counts = quantile_matrix,assay = 'RNA',project = 'temp',min.cells = 0,min.features = 0)
quantile_matrix$species <- NA
quantile_matrix@meta.data[temp_human,"species"] <- 'human'
quantile_matrix@meta.data[temp_macaque,"species"] <- 'macaque'

VlnPlot(object = quantile_matrix,features = gene_list,pt.size = 0,assay = 'RNA',slot = 'counts',group.by = 'species')

#get filtered gene list
temp_human <- rowMeans(quantile_matrix@assays$RNA@counts[,quantile_matrix$species == 'human'])
temp_macaque <- rowMeans(quantile_matrix@assays$RNA@counts[,quantile_matrix$species == 'macaque'])
gene_list <- names(temp_human)[temp_human > temp_macaque]

pdf(file = './res/step_79_fig_221112/up_gene_by_IP_human_specific_SNC_vlnplot.pdf',width = 16,height = 16)
VlnPlot(object = quantile_matrix,features = gene_list,pt.size = 0,assay = 'RNA',slot = 'counts',
        group.by = 'species',cols = color_param$species[c('human','macaque')])
dev.off()

#show some case
DefaultAssay(macaque_multiome_Seurat) <- 'RNA'
DefaultAssay(Greenleaf_RNA_Seurat) <- 'RNA'
p1 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('CHMP4B'),pt.size = 0.1,slot = 'data') + theme(aspect.ratio = 1)
p2 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('CHMP4B'),pt.size = 0.1,slot = 'data') + theme(aspect.ratio = 1)

p1+p2+plot_layout(ncol = 2)

#save data
saveRDS(object = gene_list,file = './res/step_79_fig_221112/up_gene_by_IP_human_specific_SNC.rds')

## try using gene score ----------------------------------------------------
#load data
gene_list <- readRDS(file = './res/step_79_fig_221112/gene_regulated_by_IP_human_specific_SNC.rds')

#get human gene score
human_gene_score_matrix <- getMatrixFromProject(ArchRProj = Greenleaf_ATAC_ArchR,useMatrix = 'GeneScoreMatrix',verbose = TRUE)
temp <- human_gene_score_matrix@assays@data$GeneScoreMatrix
temp <- temp[,which(Greenleaf_ATAC_ArchR@cellColData[colnames(temp),"cell_type"] == 'IP')]
idx <- which(human_gene_score_matrix@elementMetadata$name %in% gene_list)
idx <- idx[!duplicated(human_gene_score_matrix@elementMetadata$name[idx])]

temp <- temp[idx,]
rownames(temp) <- human_gene_score_matrix@elementMetadata$name[idx]
human_gene_score_matrix <- temp
gc()

#get macaque gene score
macaque_gene_score_matrix <- getMatrixFromProject(ArchRProj = macaque_multiome_ArchR,useMatrix = 'GeneScoreMatrix',verbose = TRUE)
temp <- macaque_gene_score_matrix@assays@data$GeneScoreMatrix
temp <- temp[,which(macaque_multiome_ArchR@cellColData[colnames(temp),"cell_type"] == 'IP')]
idx <- which(macaque_gene_score_matrix@elementMetadata$name %in% gene_list)
idx <- idx[!duplicated(macaque_gene_score_matrix@elementMetadata$name[idx])]

temp <- temp[idx,]
rownames(temp) <- macaque_gene_score_matrix@elementMetadata$name[idx]
macaque_gene_score_matrix <- temp
gc()

#human gene score up regulated in human
dim(human_gene_score_matrix)
dim(macaque_gene_score_matrix)
temp <- dplyr::intersect(rownames(human_gene_score_matrix),rownames(macaque_gene_score_matrix))
human_gene_score_matrix <- human_gene_score_matrix[temp,]
macaque_gene_score_matrix <- macaque_gene_score_matrix[temp,]

gene_score_matrix <- cbind(human_gene_score_matrix,macaque_gene_score_matrix)
gene_score_matrix <- CreateSeuratObject(counts = gene_score_matrix,project = 'temp',assay = 'RNA',min.cells = 0,min.features = 0)
gene_score_matrix$species <- NA
gene_score_matrix@meta.data[colnames(human_gene_score_matrix),"species"] <- 'human'
gene_score_matrix@meta.data[colnames(macaque_gene_score_matrix),"species"] <- 'macaque'

VlnPlot(object = gene_score_matrix,features = temp,pt.size = 0,assay = 'RNA',slot = 'counts',group.by = 'species')

temp_human <- rowMeans(human_gene_score_matrix)
temp_macaque <- rowMeans(macaque_gene_score_matrix)
gene_list <- names(temp_human)[temp_human > temp_macaque]

#save data
saveRDS(object = gene_list,file = './res/step_79_fig_221112/up_gene_score_by_IP_human_specific_SNC.rds')

#overlap between gene score and gene quantile
# gene_quantile <- readRDS(file = './res/step_79_fig_221112/up_gene_by_IP_human_specific_SNC.rds')
# gene_score <- readRDS(file = './res/step_79_fig_221112/up_gene_score_by_IP_human_specific_SNC.rds')
# gene_list <- dplyr::intersect(gene_quantile,gene_score)

gene_list <- c('RNASEH2B','CHMP2A','CSTB','NR2F1','NEUROD6','GOLM1')
gene_anno <- getGeneAnnotation(ArchRProj = Greenleaf_ATAC_ArchR)
gene_anno$genes <- gene_anno$genes[gene_anno$genes$symbol %in% gene_list]

# p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,geneSymbol = gene_list,groupBy = 'cell_type',
#                       useGroups = c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4'),
#                       features = getPeakSet(Greenleaf_ATAC_ArchR),
#                       upstream = 100000,downstream = 100000,
#                       geneAnnotation = gene_anno,
#                       sizes = c(10,1.5,2,2))
# 
# grid::grid.newpage()
# grid::grid.draw(p$CELSR1)

#get loop
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g_loop <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 5000,
  returnLoops = TRUE
)

# temp <- which(p2g@metadata$geneSet$name %in% gene_list)
# temp <- which(p2g$idxRNA %in% temp)
# temp <- p2g$idxATAC[temp]
# temp <- unique(temp)
# temp <- p2g@metadata$peakSet[temp]
# names(temp) <- paste(temp@seqnames,as.character(temp@ranges),sep = '-')
# 
# SNC_list <- DF_list[DF_list$group == 'human_specific',"human_SNC"]
# names(human_SNC) <- human_SNC$human_SNC
# SNC_list <- human_SNC[SNC_list]
# 
# idx <- countOverlaps(query = temp,subject = SNC_list)
# temp <- temp[names(idx)[idx > 0]]

# names(p2g_loop$Peak2GeneLinks) <- paste(p2g_loop$Peak2GeneLinks@seqnames,as.character(p2g_loop$Peak2GeneLinks@ranges),sep = '-')
# p2g_loop_data_frame <- rtracklayer::as.data.frame(p2g_loop$Peak2GeneLinks)
# options(future.globals.maxSize = 200*(1024^3))
# plan(multisession,workers = 6)
# p2g_loop_data_frame <- unlist(future_lapply(X = names(p2g_loop$Peak2GeneLinks),FUN = function(x){
#   
#   #create start site and end site
#   temp_p2g <- p2g_loop_data_frame[x,,drop = FALSE]
#   
#   start_site <- data.frame(chrom = temp_p2g$seqnames,start = temp_p2g$start,end = temp_p2g$start)
#   start_site <- as(start_site,'GRanges')
#   end_site <- data.frame(chrom = temp_p2g$seqnames,start = temp_p2g$end,end = temp_p2g$end)
#   end_site <- as(end_site,'GRanges')
#   
#   #intersect with peak set
#   i <- countOverlaps(query = start_site,subject = temp)
#   j <- countOverlaps(query = end_site,subject = temp)
#   if(i > 0 | j > 0){
#     return(x)
#   }else{
#     return(NA)
#   }
# }))
# plan(sequential)
# p2g_loop_data_frame <- p2g_loop_data_frame[!is.na(p2g_loop_data_frame)]
# p2g_loop$Peak2GeneLinks <- p2g_loop$Peak2GeneLinks[p2g_loop_data_frame]

names(human_SNC) <- human_SNC$human_SNC
p <- plotBrowserTrack(ArchRProj = Greenleaf_ATAC_ArchR,geneSymbol = gene_list,groupBy = 'cell_type',
                      useGroups = c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4'),
                      features = human_SNC[DF_list[DF_list$group == 'human_specific',"human_SNC"]],
                      upstream = 200000,downstream = 200000,
                      geneAnnotation = gene_anno,
                      loops = p2g_loop,
                      sizes = c(10,1.5,1.5,1.5))

pdf(file = './res/step_79_fig_221112/human_IP_NEUROD6_trackplot.pdf',width = 7,height = 4.5)
grid::grid.newpage()
grid::grid.draw(p$NEUROD6)
dev.off()

#macaque trackplot
gene_anno <- getGeneAnnotation(ArchRProj = macaque_multiome_ArchR)
gene_anno$genes <- gene_anno$genes[gene_anno$genes$symbol %in% gene_list]

p2g <- getPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

p2g_loop <- getPeak2GeneLinks(
  ArchRProj = macaque_multiome_ArchR,
  corCutOff = 0.45,
  resolution = 8000,
  returnLoops = TRUE
)

names(macaque_SNC) <- macaque_SNC$human_SNC

p <- plotBrowserTrack(ArchRProj = macaque_multiome_ArchR,geneSymbol = gene_list,groupBy = 'cell_type',
                      useGroups = c('RG-1','RG-2','Cycling','IP','Ex-1','Ex-2','Ex-3','Ex-4'),
                      features = macaque_SNC[DF_list[DF_list$group == 'human_specific',"human_SNC"]],
                      upstream = 200000,downstream = 200000,
                      geneAnnotation = gene_anno,
                      loops = p2g_loop,
                      sizes = c(10,1.5,1.5,1.5))

pdf(file = './res/step_79_fig_221112/macaque_IP_NEUROD6_trackplot.pdf',width = 7,height = 4.5)
grid::grid.newpage()
grid::grid.draw(p$NEUROD6)
dev.off()

#gene expression
DefaultAssay(object = macaque_multiome_Seurat) <- 'RNA'
DefaultAssay(object = Greenleaf_RNA_Seurat) <- 'RNA'
Idents(macaque_multiome_Seurat) <- 'cell_type'
Idents(Greenleaf_RNA_Seurat) <- 'ReAnno_celltype'

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('NEUROD6'),pt.size = 0.1,slot = 'data',label = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + 
  labs(title = 'Human')

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('NEUROD6'),pt.size = 0.1,slot = 'data',label = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'Macaque')

pdf(file = './res/step_79_fig_221112/NEUROD6_gene_expression_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

#using gene quantile
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

DefaultAssay(object = macaque_multiome_Seurat) <- 'quantile'
DefaultAssay(object = Greenleaf_RNA_Seurat) <- 'quantile'
Idents(macaque_multiome_Seurat) <- 'cell_type'
Idents(Greenleaf_RNA_Seurat) <- 'ReAnno_celltype'

p1 <- FeaturePlot(object = Greenleaf_RNA_Seurat,features = c('NEUROD6'),pt.size = 0.1,slot = 'data',label = TRUE) + 
  theme(aspect.ratio = 1,legend.position = 'none') + 
  labs(title = 'Human')

p2 <- FeaturePlot(object = macaque_multiome_Seurat,features = c('NEUROD6'),pt.size = 0.1,slot = 'data',label = TRUE) + 
  theme(aspect.ratio = 1) + 
  labs(title = 'Macaque')

pdf(file = './res/step_79_fig_221112/NEUROD6_gene_quantile_featureplot.pdf',width = 12,height = 6)
p1+p2+plot_layout(ncol = 2)
dev.off()

# for loop generate all gene regulated by human SNC -----------------------
#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#load gene quantile matrix
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

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

human_peakset <- p2g@metadata$peakSet
names(human_peakset) <- paste(human_peakset@seqnames,as.character(human_peakset@ranges),sep = '-')
idxATAC <- names(human_peakset)

human_geneset <- p2g@metadata$geneSet
idxRNA <- as.character(human_geneset$name)

p2g_data_frame <- data.frame(idxATAC = idxATAC[p2g$idxATAC],idxRNA = idxRNA[p2g$idxRNA])

#generate cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]

#for loop
for (i in cell_type_list) {
  #modify cell type
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load DF list
  DF_list <- list.files(path = './res/step_78_fig_221111/DF_list')
  DF_list <- DF_list[grep(pattern = cell_type_dot,x = DF_list,fixed = TRUE)]
  DF_list <- paste('./res/step_78_fig_221111/DF_list',DF_list,sep = '/')
  DF_list <- readRDS(DF_list)
  
  #get gene
  SNC_list <- DF_list[DF_list$group == 'human_specific',"human_SNC"]
  names(human_SNC) <- human_SNC$human_SNC
  SNC_list <- human_SNC[SNC_list]
  
  peak_list <- countOverlaps(query = human_peakset,subject = SNC_list)
  peak_list <- names(peak_list)[peak_list > 0]
  
  gene_list <- p2g_data_frame[p2g_data_frame$idxATAC %in% peak_list,"idxRNA"]
  
  #save data
  gene_list <- unique(gene_list)
  char <- paste0('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/',cell_type_dot,'_gene_regulated_by_human_specific_SNC.rds')
  saveRDS(object = gene_list,file = char)
  
  char <- paste(cell_type_dot,'done!',sep = ' ')
  print(char)
}

# for loop generate all up gene regulated by human SNC --------------------
#load data
human_SNC <- readRDS(file = './res/step_78_fig_221111/human_SNC.rds')
macaque_SNC <- readRDS(file = './res/step_78_fig_221111/macaque_SNC.rds')

Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

Greenleaf_RNA_Seurat <- readRDS(file = './processed_data/221008_summary/Greenleaf_RNA_Seurat_human_symbol_220917.rds')
macaque_multiome_Seurat <- readRDS(file = './processed_data/221008_summary/macaque_multiome_Seurat_220802.rds')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#load gene quantile matrix
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

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

human_peakset <- p2g@metadata$peakSet
names(human_peakset) <- paste(human_peakset@seqnames,as.character(human_peakset@ranges),sep = '-')
idxATAC <- names(human_peakset)

human_geneset <- p2g@metadata$geneSet
idxRNA <- as.character(human_geneset$name)

p2g_data_frame <- data.frame(idxATAC = idxATAC[p2g$idxATAC],idxRNA = idxRNA[p2g$idxRNA])

#generate cell type list
cell_type_list <- names(color_param$celltype)
cell_type_list <- cell_type_list[cell_type_list %in% Greenleaf_ATAC_ArchR$cell_type & cell_type_list %in% macaque_multiome_ArchR$cell_type]

#for loop
for (i in cell_type_list) {
  #cell_type
  cell_type <- i
  cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)
  
  #load gene list
  gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC')
  gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
  gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
  gene_list <- readRDS(gene_list)
  gene_list <- gene_list[gene_list %in% rownames(macaque_multiome_Seurat@assays$quantile@data) & gene_list %in% rownames(Greenleaf_RNA_Seurat@assays$quantile@data)]
  
  char <- paste(cell_type,'gene length',length(gene_list),sep = ' ')
  print(char)
  
  if(length(gene_list) > 0){
    temp_human <- Greenleaf_RNA_Seurat@assays$quantile@data[gene_list,Greenleaf_RNA_Seurat$ReAnno_celltype == cell_type]
    temp_macaque <- macaque_multiome_Seurat@assays$quantile@data[gene_list,macaque_multiome_Seurat$cell_type == cell_type]
    temp_human <- rowMeans(temp_human)
    temp_macaque <- rowMeans(temp_macaque)
    
    #get filtered gene list
    gene_list <- gene_list[temp_human > temp_macaque]
    
    #save data
    char <- paste0('./res/step_79_fig_221112/up_gene_regulated_by_human_specific_SNC/',cell_type_dot,'_up_gene_regulated_by_human_specific_SNC.rds')
    saveRDS(object = gene_list,file = char)
    
    char <- paste(cell_type,'done!',sep = ' ')
    print(char)
  }
}

# GO analysis for some cell type ------------------------------------------

## RG-1 --------------------------------------------------------------------
cell_type <- 'RG-1'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1,5,9,14,17),]
allRes$Term <- as.character(allRes$Term)
allRes[2,"Term"] <- 'regulation of epithelial cell differentiation'
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## RG-2 --------------------------------------------------------------------
cell_type <- 'RG-2'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(2,3,15,19,21),]
allRes$Term <- as.character(allRes$Term)
allRes[2,"Term"] <- 'regulation of endothelial cell development'
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## IP --------------------------------------------------------------------
cell_type <- 'IP'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes[c(1,2,3,5,6),],aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## Ex-1 --------------------------------------------------------------------
cell_type <- 'Ex-1'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1,3,4,6,7),]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## Ex-2 --------------------------------------------------------------------
cell_type <- 'Ex-2'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1:5),]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## Ex-3 --------------------------------------------------------------------
cell_type <- 'Ex-3'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1,2,14,16,23),]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## Ex-4 --------------------------------------------------------------------
cell_type <- 'Ex-4'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1,2,4,6,13),]
allRes$Term <- as.character(allRes$Term)
allRes[4,"Term"] <- 'hindbrain radial glia guided cell migration'
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## OPC --------------------------------------------------------------------
cell_type <- 'OPC'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1,2,4,5,15),]
allRes$Term <- as.character(allRes$Term)
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## InMGE --------------------------------------------------------------------
cell_type <- 'InMGE'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1,4,5,7,11),]
allRes$Term <- as.character(allRes$Term)
allRes[1,"Term"] <- 'positive regulation of fibroblast proliferation'
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()

## InCGE --------------------------------------------------------------------
cell_type <- 'InCGE'
cell_type_dot <- gsub(pattern = '-',replacement = '.',x = cell_type,fixed = TRUE)

#load data
Greenleaf_ATAC_ArchR <- loadArchRProject(path = './processed_data/221008_summary/Greenleaf_ATAC_ArchR_221019/')
macaque_multiome_ArchR <- loadArchRProject(path = './processed_data/221008_summary/macaque_multiome_ArchR_221011/')

color_param <- readRDS(file = './data/parameter/shared_param/MetaValue_color_param_220923.rds')

#gene list 
gene_list <- list.files(path = './res/step_79_fig_221112/gene_regulated_by_human_specific_SNC/')
gene_list <- gene_list[grep(pattern = cell_type_dot,x = gene_list,fixed = TRUE)]
gene_list <- paste('./res/step_79_fig_221112/gene_regulated_by_human_specific_SNC',gene_list,sep = '/')
gene_list <- readRDS(gene_list)

#get p2g
p2g <- getPeak2GeneLinks(
  ArchRProj = Greenleaf_ATAC_ArchR,
  corCutOff = 0.45,
  resolution = 1,
  returnLoops = FALSE
)

#GO
human_all_gene <- unique(p2g@metadata$geneSet$name)
temp <- human_all_gene
human_all_gene <- c(human_all_gene %in% gene_list)
human_all_gene[human_all_gene == TRUE] <- '1'
human_all_gene[human_all_gene == FALSE] <- '0'
names(human_all_gene) <- temp
human_all_gene <- factor(x = human_all_gene,levels = c('0','1'))

GO_enrich <- new("topGOdata",
                 description = "human specific SNC regulated gene BP",
                 ontology = "BP",
                 allGenes = human_all_gene,
                 nodeSize = 10,annotationFun = annFUN.org,
                 mapping = 'org.Hs.eg.db',ID = 'symbol')
resultFisher <- runTest(GO_enrich, algorithm = "classic", statistic = "fisher")
allRes <- GenTable(GO_enrich, classicFisher = resultFisher, 
                   orderBy = "resultFisher", ranksOf = "resultFisher", topNodes = 50)
allRes$sig <- -log10(as.numeric(allRes$classicFisher))
allRes[is.na(allRes$sig),"sig"] <- 30
allRes <- allRes[!duplicated(allRes$Term),]
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#modify
allRes <- allRes[c(1,3,5,9,14),]
allRes$Term <- as.character(allRes$Term)
allRes[1,"Term"] <- 'regulation of endothelial cell development'
allRes$Term <- factor(allRes$Term,levels = rev(allRes$Term))

#plot GO enrichment
char <- paste0('./res/step_79_fig_221112/GO_barplot/',cell_type_dot,'_human_specific_SNC_GO_BP_barplot.pdf')
pdf(file = char,width = 6.5,height = 2.5)
ggplot(allRes,aes(x=Term,y=sig)) + 
  geom_bar(stat = 'identity',position = 'stack',width = 0.7,
           fill = 'grey',colour = 'black') + 
  coord_flip() + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        axis.line = element_blank(),
        panel.background = element_rect(fill = NULL,colour = 'black',linewidth = 1),
        axis.ticks.y = element_blank()) + 
  ylab('-log10(pvalue)') + xlab(paste('human specific SNC in',cell_type,sep = ' '))
dev.off()