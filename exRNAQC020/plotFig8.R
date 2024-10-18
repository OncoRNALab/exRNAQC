

#' # Setup

if (rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

library(tidyverse)
library(ggplot2)
library(data.table)
library(ggrepel)
library(ggpubr)
library(cowplot)
library(gridExtra)
library(repr)
library(gplots)
library(pheatmap)
library(gridGraphics)
library(grid)
library(cowplot)
library(RColorBrewer)


#' # Import data
data = fread("20230427_data metrics_v2.txt")

data2 <- data %>%
  pivot_longer(c(5:6), 
               names_to = "RNAtype",
               values_to = "RNA",
               values_drop_na = TRUE) 

data2$time2 <- apply(data2[,5], 2, function(x) {
  minVal <- min(x, na.rm = T)
  maxVal <- max(x, na.rm = T)
  (x - minVal) / (maxVal - minVal) * 10
})

data2_miRNA <- data2 %>% 
  filter(RNAtype == "miRNA") %>% 
  filter(RNA > 0)


data2_mRNA <- data2 %>% 
  filter(RNAtype == "mRNA") %>% 
  filter(RNA > 0)



# Define colors and the number of bins for the color scale
numBins <- 15
custom_palette = c(colorRampPalette(brewer.pal(n = 7, name =  "RdBu"))(30))
colorPalette <- c("white", custom_palette[12:1])
breaksList <- seq(1,10)


# miRNA
## Normalize the data to a range of 0-10
data2_miRNA$RNA2 <- apply(data2_miRNA[,7], 2, function(x) {
  minVal <- min(x, na.rm = T)
  maxVal <- max(x, na.rm = T)
  (x - minVal) / (maxVal - minVal) * 10
})


data2_miRNA2 <- data2_miRNA[,c(3:4,8,9)]
data2_miRNA2 <- as.data.frame(data2_miRNA2)
rownames(data2_miRNA2) <- data2_miRNA$author2
colnames(data2_miRNA2) <- c("tube", "kit", "time", "miRNA")


## plot
data2_miRNA2$mRNA <- NA
data2_miRNA2 <- data2_miRNA2 %>% 
  rownames_to_column() %>% 
  rowwise() %>% 
  mutate(m = mean(c(tube, kit, time, miRNA), na.rm = T)) %>% 
  arrange(desc(m)) %>% 
  column_to_rownames("rowname") %>% 
  select(-m)

Fig8A <- pheatmap(data2_miRNA2 %>%  filter(rownames(.) != "exRNAQC"), scale = "none", cluster_rows = F, 
                  cluster_cols = F, color = colorPalette, show_colnames = T, na_col = "grey85", breaks = breaksList)

# mRNA
## Normalize the data to a range of 0-10
data2_mRNA$RNA2 <- apply(data2_mRNA[,7], 2, function(x) {
  minVal <- min(x, na.rm = T)
  maxVal <- max(x, na.rm = T)
  (x - minVal) / (maxVal - minVal) * 10
})

data2_mRNA2 <- data2_mRNA[,c(3:4,8,9)]
data2_mRNA2 <- as.data.frame(data2_mRNA2)
rownames(data2_mRNA2) <- data2_mRNA$author2
colnames(data2_mRNA2) <- c("tube", "kit", "time", "mRNA")

## plot
data2_mRNA2$miRNA <- NA
data2_mRNA2 <- data2_mRNA2 %>% 
  relocate(miRNA, .before = mRNA) %>% 
  rownames_to_column() %>% 
  rowwise() %>% 
  mutate(m = mean(c(tube, kit, time, mRNA), na.rm = T)) %>% 
  arrange(desc(m)) %>% 
  column_to_rownames("rowname") %>% 
  select(-m)

Fig8B <- pheatmap(data2_mRNA2, scale = "none", cluster_rows = F, cluster_cols = F, color = colorPalette, na_col = "grey85", breaks = breaksList)


#miRNA + mRNA
## combine data
exRNAQC_miRNA <- data2_miRNA2 %>% 
  filter(rownames(.) == "exRNAQC") %>% 
  as.data.frame()
colnames(exRNAQC_miRNA) <- c("tube", "kit", "time", "miRNA", "mRNA")

exRNAQC_mRNA <- data2_mRNA2 %>% 
  filter(rownames(.) == "exRNAQC") %>% 
  as.data.frame()
colnames(exRNAQC_mRNA) <- c("tube", "kit", "time", "miRNA", "mRNA")

exRNAQC_both <- cbind(exRNAQC_miRNA[,c(1:4)], exRNAQC_mRNA[,5])
colnames(exRNAQC_both) <- c("tube", "kit", "time", "miRNA", "mRNA")


## plot
breaksList = seq(0, 10, by = 1)
Fig8C <- pheatmap(exRNAQC_both, scale = "none", cluster_rows = F, cluster_cols = F, color = colorPalette, breaks = breaksList, show_colnames = F)

### combine plots into 1 figure
data2_miRNA3 <- 
  data2_miRNA2 %>%  filter(rownames(.) != "exRNAQC")
colnames(data2_miRNA3) <- c("tube", "kit", "time", "miRNA", "mRNA")

data2_mRNA3 <- 
  data2_mRNA2 %>%  filter(rownames(.) != "exRNAQC")
colnames(data2_mRNA3) <- c("tube", "kit", "time", "miRNA", "mRNA")

data2_all <- rbind(data2_miRNA3, exRNAQC_both)
data2_all2 <- rbind(data2_all, data2_mRNA3)

Fig8A <- pheatmap(data2_all2 %>%  filter(miRNA > -1) %>%  filter(rownames(.) != "exRNAQC"), scale = "none", cluster_rows = F, cluster_cols = F, color = colorPalette, breaks = breaksList, show_colnames = T, cellwidth=20, cellheight=20, na_col = "grey70")
Fig8C <- pheatmap(data2_all2 %>%  filter(rownames(.) == "exRNAQC"), scale = "none", cluster_rows = F, cluster_cols = F, color = colorPalette, breaks = breaksList, show_colnames = F, cellwidth=20, cellheight=20, na_col = "grey70")
Fig8B <- pheatmap(data2_all2 %>%  filter(mRNA > -1) %>%  filter(rownames(.) != "exRNAQC"), scale = "none", cluster_rows = F, cluster_cols = F, color = colorPalette, breaks = breaksList, cellwidth=20, cellheight=20, na_col = "grey70")


dev.off()
plot_list=list()
plot_list[['Fig8C']]=Fig8C[[4]]
plot_list[['Fig8A']]=Fig8A[[4]]
plot_list[['Fig8B']]=Fig8B[[4]]
cowplot::plot_grid(plot_list[['Fig8C']], plot_list[['Fig8A']], plot_list[['Fig8B']], rel_heights=c(0.16, 0.3), ncol = 1, align = "v")
