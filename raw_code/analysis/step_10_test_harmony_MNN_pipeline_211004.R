#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: test harmony MNN pipeline                                       ##
## Data: 2021.10.04                                                                ##
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
library(rliger)
library(factoextra)
library(circlize)
library(ggsci)
library(SeuratWrappers)
library(magrittr)
library(harmony)
library(networkD3)
library(htmlwidgets)
library(SeuratData)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')

# public data test --------------------------------------------------------
SeuratData::InstalledData()
data("panc8")
panc8@meta.data[1:3,]

panc8_celseq <- panc8[,panc8$tech == 'celseq']
panc8_celseq2 <- panc8[,panc8$tech == 'celseq2']
panc8_fluidigmc1 <- panc8[,panc8$tech == 'fluidigmc1']
panc8_indrop <- panc8[,panc8$tech == 'indrop']
panc8_smartseq2 <- panc8[,panc8$tech == 'smartseq2']

panc8_indrop <- my_process_seurat(object = panc8_indrop,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list_panc8 <- VariableFeatures(panc8_indrop)

#test harmony default method
test <- panc8
test <- my_process_seurat(object = test,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
test <- harmony::RunHarmony(object = test,group.by.vars = 'tech')
test <- RunUMAP(object = test,reduction = 'harmony',dims = 1:20)
p1 <- DimPlot(object = test,group.by = 'celltype',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = test,group.by = 'tech',label = FALSE)
pdf(file = './res/step_10_fig_211004/panc8_full_data_default_harmony_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout(ncol = 2)
dev.off()
#with full data
test <- my_harmony_integration(named_seurat_list = list(celseq=panc8_celseq,celseq2=panc8_celseq2,fluidigmc1=panc8_fluidigmc1,
                                                        indrop=panc8_indrop,smartseq2=panc8_smartseq2),assay = 'RNA',
                               variable_feature = gene_list_panc8,
                               var_to_regress_list = list(celseq=NULL,celseq2=NULL,fluidigmc1=NULL,indrop=NULL,smartseq2=NULL),
                               npcs = 50,reference_loading = 'indrop',integration_var = 'tech',harmony_input_dim = 20,max.iter.harmony = 50,
                               UMAP_dim = 20,resolution = 1,kmeans_init_iter_max = 200)

DimPlot(test,group.by = 'dataset')
test@meta.data[,'cell_type'] <- NA
test@meta.data[,'cell_type'] <- as.character(panc8@meta.data[colnames(test),"celltype"])

p1 <- DimPlot(object = test,group.by = 'cell_type',label = TRUE,repel = TRUE)
p2 <- DimPlot(object = test,group.by = 'dataset',label = FALSE)
pdf(file = './res/step_10_fig_211004/panc8_full_data_modified_harmony_umap.pdf',width = 16,height = 8)
p1+p2+plot_layout(ncol = 2)
dev.off()

temp <- my_MNN_label_transfer(data = test[,test$dataset == 'indrop'],
                              query = test[,test$dataset == 'smartseq2'],
                              reference_var = 'cell_type',reduction = 'pca',
                              mnn = 30,knn = 50)

p1 <- DimPlot(test[,test$dataset == 'indrop'],group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1,legend.position = 'none') + labs(title = 'indrop cell type')
p2 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1,legend.position = 'none') + labs(title = 'smartseq2 cell type')
p3 <- DimPlot(temp,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'smartseq2 predict label')

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_integration.pdf',width = 21,height = 8)
p1+p2+p3+plot_layout()
dev.off()

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_label_transfer_mnn_30_knn_50.pdf',width = 7,height = 6)
scibet::Confusion_heatmap(ori = temp$cell_type,prd = temp$predict_label) + 
  theme(aspect.ratio = 12/13,
        axis.title = element_text(size = 15)) + 
  xlab('cell_type') + ylab('predict_label')
dev.off()

score_matrix <- temp@meta.data[,grepl(pattern = '_score$',fixed = FALSE,x = colnames(temp@meta.data))]
score_matrix <- base::lapply(colnames(temp),FUN = function(x){
  return(max(score_matrix[x,]))
})
score_matrix <- unlist(score_matrix)
temp$score <- score_matrix

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_label_transfer_mnn_30_knn_50_predict_score_hist.pdf',width = 5,height = 5)
ggplot(data = data.frame(predict_score=temp$score)) + 
  geom_histogram(aes(x=predict_score),bins = 12,fill='#d3d3d3',colour='black') + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'histgram of predict score')
dev.off()
## test knn ----------------------------------------------------------------
for (i in c(5,10,30,50,100,300,500)) {
  temp <- my_MNN_label_transfer(data = test[,test$dataset == 'indrop'],
                                query = test[,test$dataset == 'smartseq2'],
                                reference_var = 'cell_type',reduction = 'pca',
                                mnn = 30,knn = i,return_query = FALSE)
  
  temp <- base::lapply(colnames(panc8_smartseq2),FUN = function(x){
    return(colnames(temp)[which.max(temp[x,])])
  })
  temp <- unlist(temp)
  temp <- paste(temp,as.character(i),sep = '_')
  
  assign(x = paste0('knn_',as.character(i)),value = temp)
  gc()
  print(paste0('knn_',as.character(i),'_done!'))
}

gc()

confusion_matrix <- my_confusion_matrix(ori = knn_5,prd = knn_10)
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_10,prd = knn_30))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_30,prd = knn_50))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_50,prd = panc8_smartseq2$celltype))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = panc8_smartseq2$celltype,prd = knn_100))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_100,prd = knn_300))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = knn_300,prd = knn_500))

nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

saveWidget(p1,file = './res/step_10_fig_211004/panc8_indrop_smartseq2_label_transfer_knn_sankey.html')


## test mnn ----------------------------------------------------------------
for (i in c(5,10,30,50,100,300,500)) {
  temp <- my_MNN_label_transfer(data = test[,test$dataset == 'indrop'],
                                query = test[,test$dataset == 'smartseq2'],
                                reference_var = 'cell_type',reduction = 'pca',
                                mnn = i,knn = 50,return_query = FALSE)
  
  temp <- base::lapply(colnames(panc8_smartseq2),FUN = function(x){
    return(colnames(temp)[which.max(temp[x,])])
  })
  temp <- unlist(temp)
  temp <- paste(temp,as.character(i),sep = '_')
  
  assign(x = paste0('mnn_',as.character(i)),value = temp)
  gc()
  print(paste0('mnn_',as.character(i),'_done!'))
}

gc()

confusion_matrix <- my_confusion_matrix(ori = mnn_5,prd = mnn_10)
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_10,prd = mnn_30))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_30,prd = panc8_smartseq2$celltype))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = panc8_smartseq2$celltype,prd = mnn_50))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_50,prd = mnn_100))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_100,prd = mnn_300))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_300,prd = mnn_500))

nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

saveWidget(p1,file = './res/step_10_fig_211004/panc8_indrop_smartseq2_label_transfer_mnn_sankey.html')

# missing public data test ------------------------------------------------
#using smartseq2 as reference
panc8_smartseq2 <- my_process_seurat(object = panc8_smartseq2,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list_panc8 <- VariableFeatures(panc8_smartseq2)
panc8_smartseq2 <- panc8_smartseq2[,!(panc8_smartseq2$celltype %in% c('activated_stellate','quiescent_stellate','schwann'))]

test <- my_harmony_integration(named_seurat_list = list(celseq=panc8_celseq,celseq2=panc8_celseq2,fluidigmc1=panc8_fluidigmc1,
                                                        indrop=panc8_indrop,smartseq2=panc8_smartseq2),assay = 'RNA',
                               variable_feature = gene_list_panc8,
                               var_to_regress_list = list(celseq=NULL,celseq2=NULL,fluidigmc1=NULL,indrop=NULL,smartseq2=NULL),
                               npcs = 50,reference_loading = 'smartseq2',integration_var = 'tech',harmony_input_dim = 20,max.iter.harmony = 50,
                               UMAP_dim = 20,resolution = 1,kmeans_init_iter_max = 200)

test@meta.data[,'cell_type'] <- NA
test@meta.data[,'cell_type'] <- as.character(panc8@meta.data[colnames(test),"celltype"])
test <- test[,test$dataset %in% c('smartseq2','indrop')]

p1 <- DimPlot(test,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme_classic() +theme(aspect.ratio = 1)
p2 <- DimPlot(test,group.by = 'dataset',label = FALSE) + theme_classic() +theme(aspect.ratio = 1)

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_missing_cell_integration.pdf',width = 16,height = 8)
p1+p2+plot_layout()
dev.off()

#label transfer
temp <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                              query = test[,test$dataset == 'indrop'],
                              reference_var = 'cell_type',reduction = 'pca',
                              mnn = 80,knn = 300)

p1 <- DimPlot(object = test[,test$dataset == 'smartseq2'],group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'samrtseq2 cell type')
p2 <- DimPlot(temp,group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop cell type')
p3 <- DimPlot(temp,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop predict label')

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_missing_cell_mnn_label_transfer.pdf',width = 24,height = 8)
p1+p2+p3+plot_layout()
dev.off()

## test mnn ----------------------------------------------------------------
#small mnn results more locally results
for (i in c(5,10,30,80,100,300,500)) {
  temp <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                query = test[,test$dataset == 'indrop'],
                                reference_var = 'cell_type',reduction = 'pca',
                                mnn = i,knn = 300,return_query = FALSE)
  
  temp <- base::lapply(colnames(panc8_indrop),FUN = function(x){
    return(colnames(temp)[which.max(temp[x,])])
  })
  temp <- unlist(temp)
  temp <- paste(temp,as.character(i),sep = '_')
  
  assign(x = paste0('mnn_',as.character(i)),value = temp)
  gc()
  print(paste0('mnn_',as.character(i),'_done!'))
}

gc()

confusion_matrix <- my_confusion_matrix(ori = mnn_5,prd = mnn_10)
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_10,prd = mnn_30))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_30,prd = mnn_80))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_80,prd = panc8_indrop$celltype))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = panc8_indrop$celltype,prd = mnn_100))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_100,prd = mnn_300))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = mnn_300,prd = mnn_500))

nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

saveWidget(p1,file = './res/step_10_fig_211004/panc8_indrop_smartseq2_missing_cell_label_transfer_mnn_sankey.html')

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_missing_cell_mnn_label_transfer_confusion_heatmap.pdf',width = 5,height = 6)
scibet::Confusion_heatmap(ori = panc8_indrop$celltype,prd = mnn_80) + 
  theme(aspect.ratio = 11/13,
        axis.title = element_text(size = 14)) + 
  xlab('ori cell type') + ylab('predict label')
dev.off()
#mnn=30,knn=50
mnn_50 <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                query = test[,test$dataset == 'indrop'],
                                reference_var = 'cell_type',reduction = 'pca',
                                mnn = 50,knn = 300)
mnn_100 <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                 query = test[,test$dataset == 'indrop'],
                                 reference_var = 'cell_type',reduction = 'pca',
                                 mnn = 100,knn = 300)
mnn_80 <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                query = test[,test$dataset == 'indrop'],
                                reference_var = 'cell_type',reduction = 'pca',
                                mnn = 80,knn = 300)

gc()

p1 <- DimPlot(object = test[,test$dataset == 'smartseq2'],group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'samrtseq2 cell type')
p2 <- DimPlot(mnn_50,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop mnn 50')
p3 <- DimPlot(mnn_80,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop mnn 80')
p4 <- DimPlot(mnn_100,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop mnn 100')

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_missing_cell_mnn_30_50_300_label_transfer.pdf',width = 16,height = 12)
p1+p2+p3+p4+plot_layout(ncol = 2)
dev.off()

# benchmark iteration -----------------------------------------------------
panc8_celseq <- panc8[,panc8$tech == 'celseq']
panc8_celseq2 <- panc8[,panc8$tech == 'celseq2']
panc8_fluidigmc1 <- panc8[,panc8$tech == 'fluidigmc1']
panc8_indrop <- panc8[,panc8$tech == 'indrop']
panc8_smartseq2 <- panc8[,panc8$tech == 'smartseq2']

panc8_smartseq2 <- my_process_seurat(object = panc8_smartseq2,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list_panc8 <- VariableFeatures(panc8_smartseq2)
panc8_smartseq2 <- panc8_smartseq2[,!(panc8_smartseq2$celltype %in% c('activated_stellate','quiescent_stellate','schwann'))]

test <- my_harmony_integration(named_seurat_list = list(celseq=panc8_celseq,celseq2=panc8_celseq2,fluidigmc1=panc8_fluidigmc1,
                                                        indrop=panc8_indrop,smartseq2=panc8_smartseq2),assay = 'RNA',
                               variable_feature = gene_list_panc8,
                               var_to_regress_list = list(celseq=NULL,celseq2=NULL,fluidigmc1=NULL,indrop=NULL,smartseq2=NULL),
                               npcs = 50,reference_loading = 'smartseq2',integration_var = 'tech',harmony_input_dim = 20,max.iter.harmony = 50,
                               UMAP_dim = 20,resolution = 1,kmeans_init_iter_max = 200)

DimPlot(test,group.by = 'dataset')
test@meta.data[,'cell_type'] <- NA
test@meta.data[,'cell_type'] <- as.character(panc8@meta.data[colnames(test),"celltype"])


for (i in 0:6) {
  temp <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                query = test[,test$dataset == 'indrop'],
                                reference_var = 'cell_type',reduction = 'pca',
                                mnn = 80,knn = 300,Iteration = i,return_query = FALSE)
  temp <- base::lapply(colnames(panc8_indrop),FUN = function(x){
    return(colnames(temp)[which.max(temp[x,])])
  })
  temp <- unlist(temp)
  temp <- paste(temp,as.character(i),sep = '_')
  
  assign(x = paste0('Iteration_',as.character(i)),value = temp)
  gc()
  print(paste0('Iteration_',as.character(i),'_done!'))
}

gc()

confusion_matrix <- my_confusion_matrix(ori = Iteration_0,prd = Iteration_1)
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = Iteration_1,prd = Iteration_2))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = Iteration_2,prd = Iteration_3))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = Iteration_3,prd = Iteration_4))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = Iteration_4,prd = Iteration_5))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = Iteration_5,prd = panc8_indrop$celltype))
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = panc8_indrop$celltype,prd = Iteration_6))

nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

saveWidget(p1,file = './res/step_10_fig_211004/panc8_indrop_smartseq2_label_transfer_Iteration_sankey.html')


Iteration_0 <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                     query = test[,test$dataset == 'indrop'],
                                     reference_var = 'cell_type',reduction = 'pca',
                                     mnn = 50,knn = 300,Iteration = 0)

Iteration_3 <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                     query = test[,test$dataset == 'indrop'],
                                     reference_var = 'cell_type',reduction = 'pca',
                                     mnn = 50,knn = 300,Iteration = 3)

Iteration_5 <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                     query = test[,test$dataset == 'indrop'],
                                     reference_var = 'cell_type',reduction = 'pca',
                                     mnn = 50,knn = 300,Iteration = 5)
gc()

p1 <- DimPlot(object = test[,test$dataset == 'smartseq2'],group.by = 'cell_type',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'samrtseq2 cell type')
p2 <- DimPlot(Iteration_0,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop iteration 0')
p3 <- DimPlot(Iteration_3,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop iteration 3')
p4 <- DimPlot(Iteration_5,group.by = 'predict_label',label = TRUE,repel = TRUE) + theme(aspect.ratio = 1) + labs(title = 'indrop iteration 5')

pdf(file = './res/step_10_fig_211004/panc8_indrop_smartseq2_label_transfer_umap.pdf',width = 18,height = 12)
p1+p2+p3+p4+plot_layout()
dev.off()

# gaussian vs exp ---------------------------------------------------------
my_MNN_label_transfer_exp <- function(data,query,reference_var,reduction='pca',mnn=30,knn=50,Iteration=5,return_query=TRUE){
  
  require(Seurat)
  require(dplyr)
  
  #create anchor matrix
  anchor_matrix <- my_low_dim_anchor_construction(data = data,query = query,reference_var = reference_var,reduction = reduction,knn = mnn)
  
  #create knn graph
  nearest <- My_FindKNN(data = query@reductions[[reduction]]@cell.embeddings,
                        query = query@reductions[[reduction]]@cell.embeddings,
                        k = knn+1)
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1]
  nearest$nn.cell <- nearest$nn.cell[,-1]
  
  #modify knn graph
  nearest_matrix <- do.call(rbind,base::lapply(colnames(query),FUN = function(x){
    temp <- rep(0,length(colnames(query)))
    names(temp) <- colnames(query)
    dist_list <- nearest$nn.dists[x,]
    temp[nearest$nn.cell[x,]] <- exp(-dist_list)/sum(exp(-dist_list))
    return(temp)
  }))
  
  colnames(nearest_matrix) <- colnames(query)
  rownames(nearest_matrix) <- colnames(query)
  nearest_matrix <- as.matrix(nearest_matrix)
  
  #predict
  predict_matrix <- anchor_matrix
  if(Iteration){
    for (i in 1:Iteration) {
      predict_matrix <- nearest_matrix %*% as.matrix(predict_matrix[colnames(query),])
    }
  }
  rownames(predict_matrix) <- colnames(query)
  colnames(predict_matrix) <- colnames(anchor_matrix)
  gc()
  
  if(return_query){
    label <- base::lapply(colnames(query),FUN = function(x){
      return(which.max(predict_matrix[x,]))
    })
    label <- unlist(label)
    label <- colnames(predict_matrix)[label]
    query$predict_label <- label
    
    colnames(predict_matrix) <- paste(colnames(predict_matrix),'score',sep = '_')
    predict_matrix <- predict_matrix[colnames(query),]
    query@meta.data <- cbind(query@meta.data,predict_matrix)
    return(query)
  }else{
    return(predict_matrix)
  }
}

guassian <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                  query = test[,test$dataset == 'indrop'],
                                  reference_var = 'cell_type',reduction = 'pca',
                                  mnn = 80,knn = 300,Iteration = 5,return_query = FALSE)

guassian <- base::lapply(colnames(panc8_indrop),function(x){
  return(colnames(guassian)[which.max(guassian[x,])])
})
guassian <- unlist(guassian)
guassian <- paste(guassian,'gaussian',sep = '_')

exp <- my_MNN_label_transfer_exp(data = test[,test$dataset == 'smartseq2'],
                                 query = test[,test$dataset == 'indrop'],
                                 reference_var = 'cell_type',reduction = 'pca',
                                 mnn = 80,knn = 100,Iteration = 2,return_query = FALSE)

exp <- base::lapply(colnames(panc8_indrop),function(x){
  return(colnames(exp)[which.max(exp[x,])])
})
exp <- unlist(exp)
exp <- paste(exp,'exp',sep = '_')

confusion_matrix <- my_confusion_matrix(ori = guassian,prd = panc8_indrop$celltype)
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = panc8_indrop$celltype,prd = exp))

nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)



#see nearest neighthor weight
my_MNN_label_transfer <- function(data,query,reference_var,reduction='pca',mnn=30,knn=50,Iteration=5,return_query=TRUE){
  
  require(Seurat)
  require(dplyr)
  
  #create knn graph
  nearest <- My_FindKNN(data = query@reductions[[reduction]]@cell.embeddings,
                        query = query@reductions[[reduction]]@cell.embeddings,
                        k = knn+1)
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1]
  nearest$nn.cell <- nearest$nn.cell[,-1]
  
  #modify knn graph
  nearest_matrix <- do.call(rbind,base::lapply(colnames(query),FUN = function(x){
    temp <- rep(0,length(colnames(query)))
    names(temp) <- colnames(query)
    dist_list <- nearest$nn.dists[x,]
    if(var(dist_list) == 0){
      temp[nearest$nn.cell[x,]] <- rep(x=1/knn,times = knn)
    }else if(sum(exp(-(dist_list^2)/(2*var(dist_list)))) == 0){
      temp[nearest$nn.cell[x,]] <- rep(x=1/knn,times = knn)
    }else{
      temp[nearest$nn.cell[x,]] <- exp(-(dist_list^2)/(2*var(dist_list)))/sum(exp(-(dist_list^2)/(2*var(dist_list))))
    }
    return(temp)
  }))
  
  colnames(nearest_matrix) <- colnames(query)
  rownames(nearest_matrix) <- colnames(query)
  nearest_matrix <- as.matrix(nearest_matrix)
  
  #predict
  predict_matrix <- nearest_matrix
  if(Iteration){
    for (i in 1:Iteration) {
      predict_matrix <- nearest_matrix %*% as.matrix(predict_matrix[colnames(query),])
    }
  }
  return(predict_matrix)
}

for (i in 0:6) {
  temp <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                                query = test[,test$dataset == 'indrop'],
                                reference_var = 'cell_type',reduction = 'pca',
                                mnn = 80,knn = 300,Iteration = i,return_query = FALSE)
  temp <- base::lapply(colnames(panc8_indrop),FUN = function(x){
    return(max(temp[x,]))
  })
  temp <- unlist(temp)
  
  assign(x = paste0('Iteration_',as.character(i)),value = temp)
  gc()
  print(paste0('Iteration_',as.character(i),'_done!'))
}

gc()

p1 <- ggplot(data = data.frame(max_weight=Iteration_0)) + geom_histogram(aes(x=max_weight),binwidth = 0.02) + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Iteration 1')
p2 <- ggplot(data = data.frame(max_weight=Iteration_1)) + geom_histogram(aes(x=max_weight),binwidth = 0.02) + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Iteration 2')
p3 <- ggplot(data = data.frame(max_weight=Iteration_2)) + geom_histogram(aes(x=max_weight),binwidth = 0.02) + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Iteration 3')
p4 <- ggplot(data = data.frame(max_weight=Iteration_3)) + geom_histogram(aes(x=max_weight),binwidth = 0.02) + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Iteration 4')
p5 <- ggplot(data = data.frame(max_weight=Iteration_4)) + geom_histogram(aes(x=max_weight),binwidth = 0.02) + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Iteration 5')
p6 <- ggplot(data = data.frame(max_weight=Iteration_5)) + geom_histogram(aes(x=max_weight),binwidth = 0.02) + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'Iteration 6')

pdf(file = './res/step_10_fig_211004/gaussian_iteration_histgram.pdf',width = 12,height = 6)
p1+p2+p3+p4+p5+p6+plot_layout(ncol = 3)
dev.off()

#exp
my_MNN_label_transfer_exp <- function(data,query,reference_var,reduction='pca',mnn=30,knn=50,Iteration=5,return_query=TRUE){
  
  require(Seurat)
  require(dplyr)
  
  #create anchor matrix
  anchor_matrix <- my_low_dim_anchor_construction(data = data,query = query,reference_var = reference_var,reduction = reduction,knn = mnn)
  
  #create knn graph
  nearest <- My_FindKNN(data = query@reductions[[reduction]]@cell.embeddings,
                        query = query@reductions[[reduction]]@cell.embeddings,
                        k = knn+1)
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1]
  nearest$nn.cell <- nearest$nn.cell[,-1]
  
  #modify knn graph
  nearest_matrix <- do.call(rbind,base::lapply(colnames(query),FUN = function(x){
    temp <- rep(0,length(colnames(query)))
    names(temp) <- colnames(query)
    dist_list <- nearest$nn.dists[x,]
    temp[nearest$nn.cell[x,]] <- exp(-dist_list)/sum(exp(-dist_list))
    return(temp)
  }))
  
  colnames(nearest_matrix) <- colnames(query)
  rownames(nearest_matrix) <- colnames(query)
  nearest_matrix <- as.matrix(nearest_matrix)
  return(nearest_matrix)
}

exp <- my_MNN_label_transfer_exp(data = test[,test$dataset == 'smartseq2'],
                                 query = test[,test$dataset == 'indrop'],
                                 reference_var = 'cell_type',reduction = 'pca',
                                 mnn = 80,knn = 100,Iteration = 2,return_query = FALSE)

exp <- base::lapply(colnames(panc8_indrop),function(x){
  return(max(exp[x,]))
})
exp <- unlist(exp)

pdf(file = './res/step_10_fig_211004/exp_iteration_histgram.pdf',width = 6,height = 3)
ggplot(data = data.frame(max_weight=exp)) + geom_histogram(aes(x=max_weight),binwidth = 0.02) + 
  theme_classic() + 
  theme(aspect.ratio = 0.7,
        plot.title = element_text(hjust = 0.5)) + 
  labs(title = 'exp')
dev.off()

# compare with seurat label transfer --------------------------------------
SeuratData::InstalledData()
data("panc8")
panc8@meta.data[1:3,]

panc8_celseq <- panc8[,panc8$tech == 'celseq']
panc8_celseq2 <- panc8[,panc8$tech == 'celseq2']
panc8_fluidigmc1 <- panc8[,panc8$tech == 'fluidigmc1']
panc8_indrop <- panc8[,panc8$tech == 'indrop']
panc8_smartseq2 <- panc8[,panc8$tech == 'smartseq2']

panc8_smartseq2 <- my_process_seurat(object = panc8_smartseq2,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)
gene_list_panc8 <- VariableFeatures(panc8_smartseq2)
panc8_smartseq2 <- panc8_smartseq2[,!(panc8_smartseq2$celltype %in% c('activated_stellate','quiescent_stellate','schwann'))]
panc8_indrop <- my_process_seurat(object = panc8_indrop,assay = 'RNA',reduction.name = 'pca',variable.feature = NULL,nfeatures = 2000,vars.to.regress = NULL,npcs = 50,preprocess = TRUE)

anchor <- FindTransferAnchors(reference = panc8_smartseq2,query = panc8_indrop,dims = 1:30)
prediction_table <- TransferData(anchorset = anchor,refdata = panc8_smartseq2$celltype,dims = 1:30)
panc8_indrop <- AddMetaData(object = panc8_indrop,metadata = prediction_table)

scibet::Confusion_heatmap(ori = panc8_indrop$celltype,prd = panc8_indrop$predicted.id)

test <- my_harmony_integration(named_seurat_list = list(celseq=panc8_celseq,celseq2=panc8_celseq2,fluidigmc1=panc8_fluidigmc1,
                                                        indrop=panc8_indrop,smartseq2=panc8_smartseq2),assay = 'RNA',
                               variable_feature = gene_list_panc8,
                               var_to_regress_list = list(celseq=NULL,celseq2=NULL,fluidigmc1=NULL,indrop=NULL,smartseq2=NULL),
                               npcs = 50,reference_loading = 'smartseq2',integration_var = 'tech',harmony_input_dim = 20,max.iter.harmony = 50,
                               UMAP_dim = 20,resolution = 1,kmeans_init_iter_max = 200)

test@meta.data[,'cell_type'] <- NA
test@meta.data[,'cell_type'] <- as.character(panc8@meta.data[colnames(test),"celltype"])
test <- test[,test$dataset %in% c('smartseq2','indrop')]

temp <- my_MNN_label_transfer(data = test[,test$dataset == 'smartseq2'],
                              query = test[,test$dataset == 'indrop'],
                              reference_var = 'cell_type',reduction = 'pca',
                              mnn = 90,knn = 300)

confusion_matrix <- my_confusion_matrix(ori = paste(temp$predict_label,'my',sep = '_'),prd = temp$cell_type)
confusion_matrix <- rbind(confusion_matrix,my_confusion_matrix(ori = panc8_indrop$celltype,prd = paste(panc8_indrop$predicted.id,'seurat',sep = '_')))

nodes <- data.frame(name=c(as.character(confusion_matrix$ori), as.character(confusion_matrix$prd)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
confusion_matrix$IDori=match(confusion_matrix$ori, nodes$name)-1 
confusion_matrix$IDprd=match(confusion_matrix$prd, nodes$name)-1

# Make the Network
p1 <- sankeyNetwork(Links = confusion_matrix, Nodes = nodes,
                    Source = "IDori", Target = "IDprd",
                    Value = "value", NodeID = "name", 
                    sinksRight=FALSE, nodeWidth=40, fontSize=13, nodePadding=20)

saveWidget(p1,file = './res/step_10_fig_211004/panc8_indrop_smartseq2_missing_cell_label_transfer_compare_seurat_with_mnn.html')
