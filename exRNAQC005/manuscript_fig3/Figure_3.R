if (rstudioapi::isAvailable()) {
}

knitr::opts_chunk$set(warning=FALSE, message=FALSE)


library(tidyverse)
library(reshape2)
library(broom)
library(RColorBrewer)
library(scales)
library(ggrepel)
library(gridExtra)
library(pheatmap)
library(colorRamps)


#load data
FC_all_005 <- read_tsv("summary_fig_exRNAQC005.tsv")
FC_all_005$study_detail <- "mRNA"
FC_all_013 <- read_tsv("summary_fig_exRNAQC013.tsv")
FC_all_013$study_detail <- "small RNA"
FC_all <- rbind(FC_all_005, FC_all_013)

FC_all$variable <- gsub("miRNAfract_FC", "biotype", FC_all$variable )
FC_all$variable <- gsub("biotype_FC", "biotype", FC_all$variable )

FC_all$variable <- gsub("miRNAcount_FC", "gene/miRNA count",FC_all$variable )
FC_all$variable <- gsub("numbergenes_FC", "gene/miRNA count", FC_all$variable )

FC_all$variable <- gsub("EndovsRC_FC", "RNA concentration", FC_all$variable )
FC_all$variable <- gsub("EndoVsSeq_FC", "RNA concentration", FC_all$variable )


#colors
paletteLength <- 12
breaksList = seq.int(1, 4, by = 0.25)
# color=colorRampPalette(c("#0072B2", "white", "#E69F00"))(paletteLength)
# color=colorRampPalette(c("white", "red"))(paletteLength)
# color <- c("#FFFFFF", "#FFE9E9", "#FFD4D4", "#FFBFBF", "#FFAAAA", "#FF9494", "#FF7F7F", "#FF6A6A", "#FF5555", "#FF3F3F", "#FF0000", "#CC0000")
custom_palette = c(colorRampPalette(brewer.pal(n = 7, name =  "RdBu"))(20))
color <- c('white', custom_palette[11:1])


# #mRNA preservation + non-preservation 
FC_all_005b <- FC_all_005 %>% 
  arrange(desc(TubeType)) %>% 
  select(Tube, variable, value, TubeType) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  rowwise() %>% mutate(average=mean(c(ALC, biotype_FC, numbergenes_FC, EndoVsSeq_FC, hemolysis_FC))) %>% 
  group_by(TubeType) %>% arrange(average) %>%
  ungroup() %>% 
  select(-TubeType) %>% 
  column_to_rownames("Tube") %>% 
  as.matrix()


my_sample_col <- FC_all_005 %>%
  arrange(desc(TubeType)) %>%
  select(Tube, variable, value, TubeType) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  column_to_rownames("Tube") %>%
  select(TubeType)


Fig3_005 <- pheatmap(FC_all_005b,
                     scale = "none", cluster_rows = F,
                     cluster_cols = F, show_colnames = T,
                     cellwidth=20, cellheight=20,
                     color = color,
                     gaps_row = 5,
                     gaps_col = 5,
                     annotation_row = my_sample_col, 
                     breaks = breaksList,
                     main = "mRNA"
)



FC_all_013b <- FC_all_013 %>% 
  arrange(desc(TubeType)) %>% 
  select(Tube, variable, value, TubeType) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  rowwise() %>% mutate(average=mean(c(ALC, miRNAfract_FC, miRNAcount_FC, EndovsRC_FC, hemolysis_FC))) %>% 
  group_by(TubeType) %>% arrange(average, .by_group = T) %>%
  ungroup() %>% 
  select(-TubeType) %>% 
  column_to_rownames("Tube") %>% 
  as.matrix()

#change values >4 to 4
FC_all_013c <- FC_all_013b %>% 
  as.data.frame() %>% 
  mutate_each(funs(replace(., . > 4, 4)))

my_sample_col2 <- FC_all_013 %>%
  arrange(desc(TubeType)) %>%
  select(Tube, variable, value, TubeType) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  column_to_rownames("Tube") %>%
  select(TubeType)


Fig3_013 <- pheatmap(FC_all_013c,
                     scale = "none", cluster_rows = F,
                     cluster_cols = F, show_colnames = T,
                     cellwidth=20, cellheight=20,
                     color = color,
                     gaps_row = 5,
                     gaps_col = 5,
                     annotation_row = my_sample_col2, 
                     breaks = breaksList,
                     main = "miRNA"
)


#mRNA and miRNA
dev.off()
plot_list=list()
plot_list[['Fig3_005']]=Fig3_005[[4]]
plot_list[['Fig3_013']]=Fig3_013[[4]]
cowplot::plot_grid(plot_list[['Fig3_005']], plot_list[['Fig3_013']], ncol = 2, nrow = 1, align = "hv")
