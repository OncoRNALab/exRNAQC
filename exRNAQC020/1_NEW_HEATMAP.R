# "Repeated measures analyses: beta regression models - FRANCISCO NOV 2022"

library(nlme)
library(ggplot2)
library(multcomp)
library(emmeans)
library(glmmTMB)
library(car)

require(data.table)
require(tidyverse)
require(dplyr)
library(ggplot2)
library(ggrepel)
require(pheatmap)
require(cowplot)
require(gridExtra)
library(RColorBrewer)

set.seed(189672)


db<-read.csv("duplicates.csv")
str(db)

db$Tube<-as.factor(db$Tube)
db$RNAisolation<-as.factor(db$RNAisolation)
db$TimeLapse<-as.factor(db$TimeLapse)
db$donorID<-as.factor(db$donorID)

Tube.levels<-levels(db$Tube)
RNAisolation.levels<-levels(db$RNAisolation)
TimeLapse.levels<-levels(db$TimeLapse)

db2 = fread("EPIC.BRef_mixed_model_input.txt", data.table = FALSE)
db2$Tube<-as.factor(db2$Tube)
db2$BiologicalReplicate<-as.factor(db2$BiologicalReplicate)
db2$TimeLapse<-as.factor(db2$TimeLapse)

cell.types<-unique(db2$CT)

cnt<-1
p<-list()
for(ct in cell.types) {
  p1<-ggplot(data = db2[db2$CT==ct,], aes(x = TimeLapse, y = CT.proportion, group = Tube:BiologicalReplicate))
  p2<-p1+geom_line(aes(color=Tube))+geom_point() + facet_grid(. ~ Tube)+
    ggtitle(paste("Cell type: ",ct,sep=""))
  p[[cnt]]<-p2
  cnt<-cnt+1
}
p

dev.off()

####################################
# Based on the previous visualisations, O.Thas suggested to have different variances for the tubes. Separate analyses for the cell types.
# Beta regression models with random effects for all cell types

ALL_OUTPUT = lapply(cell.types, function(cell_type){
  
  m <- glmmTMB(CT.proportion ~ Tube*TimeLapse + (1|BiologicalReplicate),
             dispformula = ~Tube,
             family = beta_family(),
             data=db2[db2$CT==cell_type,])
  
  m.methods <- emmeans(m, ~TimeLapse|Tube)
  
  results = reshape2::melt(pairs(m.methods))
  results = results[,c(1,2,ncol(results))]
  colnames(results) = c("time_comparison", "tube", "p-value")
  results$CT = gsub("_","\\.", cell_type)
  results
  
})

ALL_OUTPUT = do.call(rbind.data.frame, ALL_OUTPUT)

ALL_OUTPUT = data.table::dcast.data.table(data.table(ALL_OUTPUT), formula = tube ~ CT + time_comparison, value.var = "p-value") %>% as.data.frame()
rownames(ALL_OUTPUT) = ALL_OUTPUT$tube
ALL_OUTPUT$tube = NULL

######################################################################################
## Re-order columns on heatmap: all T0-T1, T1-T2 and T0-T2
ALL_OUTPUT = ALL_OUTPUT[, c(grep("T0 - T1", colnames(ALL_OUTPUT)), grep("T1 - T2", colnames(ALL_OUTPUT)), grep("T0 - T2", colnames(ALL_OUTPUT)))]

annotation = data.frame(
  CT = colnames(ALL_OUTPUT) %>% strsplit(., "_") %>% lapply(., function(x) x[1]) %>% unlist(),
  time_comparison = colnames(ALL_OUTPUT) %>% strsplit(., "_") %>% lapply(., function(x) x[2]) %>% unlist()
)
rownames(annotation) = colnames(ALL_OUTPUT)

annotation_colors = list(
  time_comparison = c(`T0 - T1`="#E06C36", `T0 - T2`="#F2904B", `T1 - T2`="#51615E"),
  CT = c(Bcells="#999999", CD4.Tcells="#E69F00", CD8.Tcells="#56B4E9", Monocytes="#009E73", NKcells="#0072B2", Neutrophils="#F0E442", otherCells="#D55E00")
)

ALL_OUTPUT.backup = ALL_OUTPUT

## Transform p-values > 0.05 TO WHITE, gradient of p-values only for significant ones!
ALL_OUTPUT = ALL_OUTPUT.backup

#Reorder matrix & only keep T0-T1 and T0-T2
ALL_OUTPUT = ALL_OUTPUT[c("ACD-A","Citrate","EDTA","EDTA separator","serum","Biomatrica","DNA Streck","RNA Streck","PAXgene","Roche"),grep("T1 - T2",colnames(ALL_OUTPUT),invert=T, value = T)]

#breaks for legend
bk1 <- c(seq(0, 0.051, by = 0.005))

#Make custom color palette
custom_palette = c(colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(20))
custom_palette2 = c(custom_palette[1:10], 'white')

#make matrix of labels (to display "<0.01" instead of the actual value)
#labels_tab <- round(ALL_OUTPUT, 3)
#labels_tab[labels_tab < 0.001] <- "<0.001"
labels_tab <- ALL_OUTPUT
labels_tab[labels_tab > 0.05] <- ""
#use scientific notation
labels_tab_formatted <- apply(labels_tab, c(1, 2), function(x) {
  formatted <- sprintf("%.1e", as.numeric(x)) # Format to scientific notation
  gsub("NA","",gsub("e-0", "e-", gsub("e\\+0", "e+", formatted))) # Remove leading zero in exponent
})

pdf("./NON_BINARY_sig_pvalue_inblack_beta_reg_OCT2023.pdf", height = 4, width = 10)
  pheatmap::pheatmap(round(ALL_OUTPUT, 3), 
                     cluster_cols = FALSE, 
                     cluster_rows = TRUE,
                     annotation = annotation, 
                     annotation_colors = annotation_colors, 
                     breaks = bk1,
                     display_numbers = labels_tab, fontsize_number = 6, number_color="black",
                     na_col = "white",
                     gaps_row = 5,
                     color = custom_palette2,
                     show_colnames = FALSE,
                     angle_col = 45)
dev.off()

######################################################################################
## Make an extra heatmap where cell type proportions themselves are plotted
DT = dcast(db2, fun.aggregate = mean, by = c("BiologicalReplicate"), formula = Tube ~ TimeLapse + CT, value.var = "CT.proportion") %>% data.frame()
rownames(DT) = DT$Tube
DT$Tube = NULL

#Make custom color palette
custom_palette = c(colorRampPalette(brewer.pal(n = 7, name = "RdBu"))(20))
custom_palette2 = c(custom_palette[1:8], 'white') %>% rev()

annotation = data.frame(
  CT = colnames(DT) %>% strsplit(., "_") %>% lapply(., function(x) x[2]) %>% unlist(),
  time_comparison = colnames(DT) %>% strsplit(., "_") %>% lapply(., function(x) x[1]) %>% unlist()
)
rownames(annotation) = colnames(DT)

annotation_colors = list(
  CT = c(Bcells="#999999", CD4="#E69F00", CD8="#56B4E9", Monocytes="#009E73", NKcells="#0072B2", Neutrophils="#F0E442", otherCells="#D55E00"),
  time_comparison = c(T0="#ede2f0", T1="#7f617e", T2="#4a345b")
)


pdf("./NON_BINARY_mean.CT.prop_OCT2023.pdf", height = 4, width = 10)
pheatmap::pheatmap(round(DT, 3), 
                   cluster_cols = FALSE, 
                   cluster_rows = TRUE,
                   annotation = annotation, 
                   annotation_colors = annotation_colors, 
                   display_numbers = TRUE, fontsize_number = 6, number_color="black",
                   na_col = "white",
                   color = custom_palette2,
                   show_colnames = FALSE,
                   angle_col = 45)
dev.off()
