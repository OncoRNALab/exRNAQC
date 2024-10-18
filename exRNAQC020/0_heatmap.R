setwd("/exRNAQC020")
require(data.table)
require(tidyverse)
require(dplyr)
library(ggplot2)
library(ggrepel)
require(pheatmap)
require(cowplot)
require(gridExtra)
library(RColorBrewer)

data <- readxl::read_xls("./Deconv_data/FAC_binary_Tube_Timelapse_signif_pvalue_fromOThas.xls")

annotation <- data.frame(CT = colnames(data)[2:ncol(data)], TimeLapse = as.character(data[1,2:ncol(data)]))
annotation$CT <- gsub(annotation$CT,pattern = "\\.[0-9]+",replacement = "")

data <- data[2:nrow(data),]
data<- as.data.frame(data)
rownames(data) <- data$...1
data$...1 <- NULL

colnames(data) <- paste(annotation$CT, annotation$TimeLapse, sep = "_")
rownames(annotation) <- colnames(data)

tubes <- rownames(data)
data <- apply(data, 2, as.numeric)
rownames(data) <- tubes

annotation_colors = list(
  TimeLapse = c(`T0-T1`="#E06C36", `T0-T2`="#F2904B", `T1-T2`="#51615E"),
  CT = c(Bcells="#999999", CD4_Tcells="#E69F00", CD8_Tcells="#56B4E9", Monocytes="#009E73", Neutrophils="#F0E442", NKcells="#0072B2", otherCells="#D55E00")
  )

pdf("./BINARY_sig_pvalue_inblack_beta_reg2.pdf", height = 3)
  pheatmap::pheatmap(data, annotation = annotation, annotation_colors = annotation_colors, color = c("white", "black"), show_colnames = FALSE) #cluster_cols = FALSE
dev.off()