#####################################################################################
## Project: macaque SP brain                                                       ##
## Script Purpose: ITGA2 RNA_Scope dot number analysis                             ##
## Data: 2023.06.21                                                                ##
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
.libPaths('/content/data/sunym/software/R_lib/R_4.3.0/')
Sys.setenv(HDF5_USE_FILE_LOCKING=FALSE,RHDF5_USE_FILE_LOCKING=FALSE)

#library
library(Rmisc)
library(dplyr)
library(Seurat)
library(ggplot2)
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
library(OpenAI4R)
library(paletteer)
library(ggpattern)
library(ggrepel)
library(openxlsx)
library(readxl)

#my function
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/sc_multiomics.R')
source('https://raw.githubusercontent.com/yimingsun12138/source_list/main/genomics.R')

#initialize ArchR
addArchRThreads(threads = 5)

#initialize OpenAI
Auth_OpenAI(key = readLines('/content/script/openai_API_key'))
chat <- Init_chat_session(model = 'gpt-4-0613')

# load RNA-Scope data -----------------------------------------------------
RNA_Scope_stat <- './res/step_121_fig_230621/ITGA2.xlsx'
sheet_names <- excel_sheets(path = RNA_Scope_stat)
sheet_names <- gsub(pattern = '-',replacement = '_',x = sheet_names,fixed = TRUE)

for (i in sheet_names) {
  temp <- gsub(pattern = '_',replacement = '-',x = i,fixed = TRUE)
  temp <- read_excel(path = RNA_Scope_stat,sheet = temp)
  assign(x = i,value = as.data.frame(temp))
}

# 20x 4.889 equals 100um
# 63x 30.833 euqals 100um

# check 20x and 63x PAX6 intensity ----------------------------------------
temp <- human_ITGA2
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p1 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank()) + 
  labs(title = 'human ITGA2') + 
  coord_cartesian(ylim = c(0,0.15))

temp <- human_PAX6
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p2 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'human PAX6') + 
  coord_cartesian(ylim = c(0,0.15))

temp <- monkey_ITGA2
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p3 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'monkey ITGA2') + 
  coord_cartesian(ylim = c(0,0.15))

temp <- monkey_PAX6
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p4 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'monkey PAX6') + 
  coord_cartesian(ylim = c(0,0.15))

pdf(file = './res/step_121_fig_230621/ITGA2_raw_average_size.pdf',width = 8,height = 4)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

#normalize by scale rule
temp <- human_ITGA2
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (100/4.889)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (100/30.833)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p1 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank()) + 
  labs(title = 'human ITGA2') + 
  ylab('um^2') + 
  coord_cartesian(ylim = c(0,12.5))

temp <- human_PAX6
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (100/4.889)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (100/30.833)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p2 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'human PAX6') + 
  coord_cartesian(ylim = c(0,12.5))

temp <- monkey_ITGA2
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (100/4.889)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (100/30.833)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p3 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'monkey ITGA2') + 
  coord_cartesian(ylim = c(0,12.5))

temp <- monkey_PAX6
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (100/4.889)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (100/30.833)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p4 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'monkey PAX6') + 
  coord_cartesian(ylim = c(0,12.5))

pdf(file = './res/step_121_fig_230621/ITGA2_average_size_normalized_by_rule.pdf',width = 8,height = 4)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

#normalize by scale
temp <- human_ITGA2
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (1000/20)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (1000/63)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p1 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank()) + 
  labs(title = 'human ITGA2') + 
  ylab('unit^2') + 
  coord_cartesian(ylim = c(0,70))

temp <- human_PAX6
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (1000/20)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (1000/63)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p2 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'human PAX6') + 
  coord_cartesian(ylim = c(0,70))

temp <- monkey_ITGA2
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (1000/20)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (1000/63)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p3 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'none',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'monkey ITGA2') + 
  coord_cartesian(ylim = c(0,70))

temp <- monkey_PAX6
temp <- temp[which(temp$average_size != 'NaN'),]
temp$average_size <- as.numeric(temp$average_size)
temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (1000/20)^2
temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (1000/63)^2
temp$scale <- factor(as.character(temp$scale),levels = c('20','63'))
temp$region <- factor(temp$region,levels = c('VZ','OSVZ'))
p4 <- ggplot(data = temp,mapping = aes(x = region,y = average_size,fill = scale)) + 
  geom_boxplot(width = 0.5) + 
  scale_fill_manual(values = c('#1F77B4FF','#FF7F0EFF')) + 
  theme_bw() + 
  theme(aspect.ratio = 3,
        panel.grid = element_blank(),
        plot.title = element_text(hjust = 0.5),
        legend.position = 'right',
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) + 
  labs(title = 'monkey PAX6') + 
  coord_cartesian(ylim = c(0,70))

pdf(file = './res/step_121_fig_230621/ITGA2_average_size_normalized_by_scale.pdf',width = 8,height = 4)
p1+p2+p3+p4+plot_layout(ncol = 4)
dev.off()

# 20x/63x ratio
scale_by_rule <- base::do.call(what = base::rbind,args = base::lapply(X = c('human_ITGA2','human_PAX6','monkey_PAX6'),FUN = function(x){
  temp <- get(x)
  temp <- temp[which(temp$average_size != 'NaN'),]
  temp$average_size <- as.numeric(temp$average_size)
  temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (100/4.889)^2
  temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (100/30.833)^2
  
  VZ_ratio <- mean(temp[which(temp$region == 'VZ' & temp$scale == 20),"average_size"]) / mean(temp[which(temp$region == 'VZ' & temp$scale == 63),"average_size"])
  OSVZ_ratio <- mean(temp[which(temp$region == 'OSVZ' & temp$scale == 20),"average_size"]) / mean(temp[which(temp$region == 'OSVZ' & temp$scale == 63),"average_size"])
  return(data.frame(ratio = c(VZ_ratio,OSVZ_ratio),region = c('VZ','OSVZ'),data_set = x,by = 'rule'))
}))

scale_by_scale <- base::do.call(what = base::rbind,args = base::lapply(X = c('human_ITGA2','human_PAX6','monkey_PAX6'),FUN = function(x){
  temp <- get(x)
  temp <- temp[which(temp$average_size != 'NaN'),]
  temp$average_size <- as.numeric(temp$average_size)
  temp[which(temp$scale == 20),"average_size"] <- temp[which(temp$scale == 20),"average_size"] * (1000/20)^2
  temp[which(temp$scale == 63),"average_size"] <- temp[which(temp$scale == 63),"average_size"] * (1000/63)^2
  
  VZ_ratio <- mean(temp[which(temp$region == 'VZ' & temp$scale == 20),"average_size"]) / mean(temp[which(temp$region == 'VZ' & temp$scale == 63),"average_size"])
  OSVZ_ratio <- mean(temp[which(temp$region == 'OSVZ' & temp$scale == 20),"average_size"]) / mean(temp[which(temp$region == 'OSVZ' & temp$scale == 63),"average_size"])
  return(data.frame(ratio = c(VZ_ratio,OSVZ_ratio),region = c('VZ','OSVZ'),data_set = x,by = 'scope'))
}))

temp <- rbind(scale_by_rule,scale_by_scale)
temp$by <- factor(temp$by,levels = c('rule','scope'))

pdf(file = './res/step_121_fig_230621/20x_divide_63x_dot_size_ratio.pdf',width = 4,height = 4)
ggplot(data = temp,mapping = aes(x = by,y = ratio,fill = by)) + 
  geom_boxplot(width = 0.4) + 
  scale_y_continuous(breaks = seq(0,17,1)) + 
  theme_bw() + 
  theme(aspect.ratio = 2.5,
        panel.grid = element_blank(),
        axis.title.x = element_blank()) + 
  xlab('normalize by') + 
  ylab('20x / 63x\ndot size ratio') + 
  labs(fill = 'normalize by')
dev.off()

# use scope scale is more resonable

# dot number per DAPI region ----------------------------------------------
human_ITGA2$species <- 'human'
human_ITGA2$gene <- 'ITGA2'
human_ITGA2$DAPI <- human_DAPI$Area * human_DAPI$percent_Area / 100

human_PAX6$species <- 'human'
human_PAX6$gene <- 'PAX6'
human_PAX6$DAPI <- human_DAPI$Area * human_DAPI$percent_Area / 100

monkey_ITGA2$species <- 'monkey'
monkey_ITGA2$gene <- 'ITGA2'
monkey_ITGA2$DAPI <- monkey_DAPI$Area * monkey_DAPI$percent_Area / 100

monkey_PAX6$species <- 'monkey'
monkey_PAX6$gene <- 'PAX6'
monkey_PAX6$DAPI <- monkey_DAPI$Area * monkey_DAPI$percent_Area / 100

temp <- rbind(human_ITGA2,human_PAX6,monkey_ITGA2,monkey_PAX6)
temp[which(temp$scale == 20),"DAPI"] <- temp[which(temp$scale == 20),"DAPI"] * (1000/20)^2
temp[which(temp$scale == 63),"DAPI"] <- temp[which(temp$scale == 63),"DAPI"] * (1000/63)^2
temp$Count_Per_DAPI <- temp$Count / temp$DAPI

# pdf(file = './res/step_121_fig_230621/ITGA2_dot_count_per_DAPI_boxplot.pdf',width = 4,height = 4)
# ggplot(data = temp,mapping = aes(x = species,y = Count_Per_DAPI,fill = gene)) + 
#   geom_boxplot(width = 0.5,outlier.alpha = 0) + 
#   theme_classic() + 
#   scale_fill_manual(values = c("#208A42","#D51F26")) + 
#   theme(aspect.ratio = 1.5,
#         axis.title.x = element_blank(),
#         axis.text.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         legend.position = c(0.8,0.9)) + 
#   coord_cartesian(ylim = c(0,0.002)) + 
#   ylab('dot count per DAPI region')
# dev.off()

# bar_table <- base::do.call(what = rbind,args = base::lapply(X = c('human','monkey'),FUN = function(x){
#   temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('ITGA2','PAX6'),FUN = function(y){
#     mean_value <- mean(temp[which(temp$species == x & temp$gene == y),"Count_Per_DAPI"])
#     sd_value <- sd(temp[which(temp$species == x & temp$gene == y),"Count_Per_DAPI"])
#     return(data.frame(mean = mean_value,sd = sd_value,species = x,gene = y))
#   }))
# }))
# 
# pdf(file = './res/step_121_fig_230621/ITGA2_dot_count_per_DAPI_barplot.pdf',width = 4,height = 4)
# ggplot(data = bar_table,mapping = aes(x = species,y = mean,fill = gene)) + 
#   geom_bar(stat = 'identity',position = position_dodge(width = 0.7),width = 0.6) + 
#   geom_errorbar(mapping = aes(ymin = mean - 0.25*sd,ymax = mean + 0.25*sd),position = position_dodge(width = 0.7),width = 0.3,linewidth = 0.2) + 
#   theme_classic() + 
#   scale_fill_manual(values = c('#5fb257','#af2318')) + 
#   theme(aspect.ratio = 1.5,
#         axis.title.x = element_blank(),
#         axis.ticks.y = element_blank(),
#         axis.text.y = element_blank(),
#         legend.position = c(0.8,0.9)) + 
#   ylab('fluorescent dot count per DAPI region')
# dev.off()

bar_table <- base::do.call(what = rbind,args = base::lapply(X = c('human','monkey'),FUN = function(x){
  temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('ITGA2','PAX6'),FUN = function(y){
    temp_temp_bar <- base::do.call(what = rbind,args = base::lapply(X = c('VZ','OSVZ'),FUN = function(z){
      mean_value <- mean(temp[which(temp$species == x & temp$gene == y & temp$region == z),"Count_Per_DAPI"])
      sd_value <- sd(temp[which(temp$species == x & temp$gene == y & temp$region == z),"Count_Per_DAPI"])
      return(data.frame(mean = mean_value,sd = sd_value,species = x,gene = y,region = z))
    }))
  }))
}))

bar_table$region <- factor(bar_table$region,levels = c('VZ','OSVZ'))

pdf(file = './res/step_121_fig_230621/ITGA2_dot_count_per_DAPI_barplot.pdf',width = 4.5,height = 4)
ggplot(data = bar_table,mapping = aes(x = species,y = mean,fill = gene)) + 
  facet_wrap(facets = ~ region,ncol = 2) + 
  geom_bar(stat = 'identity',position = position_dodge(width = 0.7),width = 0.6) + 
  geom_errorbar(mapping = aes(ymin = mean - 0.25*sd,ymax = mean + 0.25*sd),position = position_dodge(width = 0.7),width = 0.3,linewidth = 0.2) + 
  theme_bw() + 
  scale_fill_manual(values = c('#5fb257','#af2318')) + 
  theme(aspect.ratio = 2,
        panel.grid = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank(),
        legend.position = 'right') + 
  ylab('fluorescent dot count per DAPI region')
dev.off()
