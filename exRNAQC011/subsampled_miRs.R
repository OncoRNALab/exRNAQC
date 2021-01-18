require(tidyverse)
require(here)

## After scp of subsampled miR counts to local:
path_to_output <- "/Users/almorlio/Dropbox (Vandesompele lab)/Personal/Projects/exRNAQC011_smallRNA/data_raw/"
path_to_annotation <- "/Users/almorlio/Dropbox (Vandesompele lab)/Basecamp/exRNAQC/exRNAQC_011 Select three kits with one input volume (small RNA)/Data analysis/Annotation exRNAQC011.csv"
path_to_rawdata <- "/Users/almorlio/Dropbox (Vandesompele lab)/Shares/Big_data_exRNAQC/exRNAQC011/exRNAQC011_01_02_04/test_subsperkit/"


robzscore <- function(data, measurevar){
  require(dplyr)
  m <- median(pull(data,paste(measurevar))) #median
  s <- stats::mad(pull(data, paste(measurevar))) #mad = median absolute deviation (already contains scaling factor of 1.4826!)
  
  robz <- (pull(data, paste(measurevar)) - m) / (s) #calculate absolute difference between every value and median, and divide by mad 
}


sample_annotation <- data.table::fread(path_to_annotation, header=T, data.table = F)
sample_annotation<-sample_annotation %>% mutate(SampleID=paste0(GraphKit,"_",TechnicalReplicate))
theme_point<-theme_classic()+theme(strip.background = element_blank())
theme_bar<-theme_classic()+theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5),strip.background = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())
theme_boxplot<-theme_classic()+theme(axis.text.x = element_text(angle=90),strip.background = element_blank(),axis.ticks.x = element_blank(),axis.title.x = element_blank(),axis.line.x = element_blank(),legend.position = "none")
color_panel<-c("#e35d6a","#5bb75b","#428bca","#e87810","#23496b","#ffbf00","#cc2028","#039748","pink","gray","darkgray")
color_panel1 <-  c("#039748","#039748","#ffbf00","#ffbf00","#e35d6a","#e35d6a","#5bb75b","#5bb75b","#428bca","#428bca","#23496b","#23496b","grey","grey","#cc2028","#cc2028","#e87810")
names(color_panel1) <- c("CCF1","CCF4","CIRC0.25","CIRC5","MAP2","MAP4","MAX0.1","MAX0.5","MIRA0.2","MIRA0.6","MIRV0.1","MIRV0.625","MIRVE0.1","MIRVE0.625","NUC0.3","NUC0.9","MIR0.2")

color_panel2 <-  c("#e35d6a","#5bb75b","#428bca","#e87810","#23496b","grey","#ffbf00","#cc2028","#039748")
names(color_panel2) <- c("MagnaPure","Maxwell","miRNeasySPAkit","miRNeasySPkit","mirVana","mirVanaE","Norgen","NucleoSpin","QIAamp")


color_panel_spikes <- c("#e35d6a","gray","#5bb75b","#428bca")
names(color_panel_spikes) <- c("LP1","MIMAT","RC1","RC2")
color_panel_spikes2 <- c("#23496b","#ffbf00","gray","#ffbf01","orange","#23496b")
names(color_panel_spikes2) <- c("LP","RC","MIMAT","RC1","RC2","LP1")

data_path <- here::here("data_raw/")
basecamp_path <- "/Users/almorlio/Dropbox (Vandesompele lab)/Basecamp/exRNAQC/exRNAQC_011 Select three kits with one input volume (small RNA)/Data analysis/"

full_nr <- format_format(big.mark = ",", decimal.mark = ".", scientific = FALSE)
source(here::here("scripts","HelperFunctions.R"))

z_score_df <- data.frame("GraphKit"= unique(sample_annotation$GraphKit)) 
data_summary <- data.frame("SampleID"= unique(sample_annotation$SampleID))
## Join the different files
miRNAs<-data.frame(MIMATID=character())
files<-list.files(path_to_rawdata,recursive=TRUE)
files<-files[grep(".*_miRs_subs.txt", files)]
for(i in 1:length(files)){
  name_sample<-gsub(".*/","",sub("-[0-9]*$","",sub("_miRs_subs.txt","",files[i])))
  tmp<-read.table(paste0(path_to_rawdata,files[i]), header=F, sep=" ")
  colnames(tmp)<-c("MIMATID",name_sample)
  miRNAs<-full_join(miRNAs,tmp,by="MIMATID")
}

miRNAs

write.table(miRNAs, "/Users/almorlio/Dropbox (Vandesompele lab)/Personal/Projects/exRNAQC011_smallRNA/data_raw/miRNAs_extrasubs.txt", quote = F, 
            na = "", row.names = F, sep="\t")
## Cut-off
### All plots: 95% SP elimination cut-off

miRNAs[is.na(miRNAs)] <- 0
double_positives<-miRNAs %>%
  mutate_if(is.numeric, funs(replace(., .< 1, 0))) #remove pseudocounts

#make table for the 15 GraphKits (kit+volume) containing the pairwise combinations of replicates + amount of SP/DP/DN (before and after filtering)
single_pos <- data.frame(GraphKit = rep(unique(sample_annotation$GraphKit),3) %>% sort(), 
                         Replicates = rep(c("RNA1-RNA2", "RNA1-RNA3","RNA2-RNA3"), 17), 
                         SP_no_filter = NA, DP_no_filter = NA, DN_no_filter = NA,   
                         SP_after_filter = NA, DP_after_filter = NA, DN_after_filter = NA,
                         filter_threshold = NA) 
single_pos$full_name <- paste(single_pos$GraphKit, single_pos$Replicates, sep="_")

#for every kit+volume combination: determine the pairwise cut-offs
for (UniqueKit in unique(sample_annotation$GraphKit)){
  sample_duplicates<-sample_annotation %>% filter(GraphKit==UniqueKit) %>%
    dplyr::select(UniqueID,GraphKit,TechnicalReplicate)
  
  if(nrow(sample_duplicates)>1){
    #print(UniqueKit)
    double_pos_sample<-double_positives %>%
      dplyr::select(MIMATID,sample_duplicates$UniqueID)
    
    samples_comb <- combn(sample_duplicates$UniqueID,2) #compare 2 of the 3 samples at a time
    for (n_col in 1:ncol(samples_comb)){
      #print(samples_comb[,n_col])
      samplename1 <- samples_comb[,n_col][1]
      sample1 <- sample_annotation[sample_annotation$UniqueID==samplename1,]$TechnicalReplicate
      samplename2 <- samples_comb[,n_col][2]
      sample2 <- sample_annotation[sample_annotation$UniqueID==samplename2,]$TechnicalReplicate
      varname <- paste0(UniqueKit,"_",sample1,"-",sample2)
      #print(varname)
      
      double_pos_subset <- double_pos_sample %>% dplyr::select(MIMATID, paste(samplename1), paste(samplename2))
      
      double_pos_subset$pos_type <- as.factor(
        ifelse(double_pos_subset[,paste(samplename1)] > 0 & 
                 double_pos_subset[,paste(samplename2)] > 0, "DP", #double positive
               ifelse(double_pos_subset[,paste(samplename1)] == 0 &
                        double_pos_subset[,paste(samplename2)] ==0 , "DN", #double negative
                      "SP"))) # else: single positive
      single_p <- double_pos_subset %>% 
        filter(pos_type=="SP") %>% 
        mutate(., max=pmax(get(samplename1), get(samplename2)))
      
      #Threshold that removes 95% of the single positives
      threshold <- round(as.numeric(paste(quantile(single_p$max,probs = 0.95, na.rm=TRUE)))+0.005, 2) #get the 95% quantile number, round UP to two decimal numbers
      
      double_pos_subset$colouring <- as.factor(ifelse(double_pos_subset[,paste(samplename1)] > threshold, "> cut-off", 
                                                      ifelse(double_pos_subset[,paste(samplename2)] > threshold, "> cut-off", "<= cut-off")))
      
      single_pos[single_pos$full_name==paste(varname),]$filter_threshold <- threshold
      
      #Single Positives
      single_pos[single_pos$full_name==paste(varname),]$SP_no_filter <- sum( 
        ((double_pos_subset[,paste(samplename1)] > 0) & (double_pos_subset[,paste(samplename2)] == 0)) | 
          ((double_pos_subset[,paste(samplename1)] == 0) & (double_pos_subset[,paste(samplename2)] > 0))
      )
      single_pos[single_pos$full_name==paste(varname),]$SP_after_filter <-  sum( 
        ((double_pos_subset[,paste(samplename1)] > threshold) & (double_pos_subset[,paste(samplename2)] == 0)) | 
          ((double_pos_subset[,paste(samplename1)] == 0) & (double_pos_subset[,paste(samplename2)] > threshold))
      )
      
      #Double Positives
      single_pos[single_pos$full_name==paste(varname),]$DP_no_filter <- sum((double_pos_subset[,paste(samplename1)] > 0) & 
                                                                              (double_pos_subset[,paste(samplename2)] > 0))
      single_pos[single_pos$full_name==paste(varname),]$DP_after_filter <- sum((double_pos_subset[,paste(samplename1)] > threshold) & 
                                                                                 (double_pos_subset[,paste(samplename2)] > threshold))
      
      #Double Negatives
      single_pos[single_pos$full_name==paste(varname),]$DN_no_filter <- sum((double_pos_subset[,paste(samplename1)] == 0) & 
                                                                              (double_pos_subset[,paste(samplename2)] == 0))
      single_pos[single_pos$full_name==paste(varname),]$DN_after_filter <- sum((double_pos_subset[,paste(samplename1)] <= threshold) & 
                                                                                 (double_pos_subset[,paste(samplename2)] <= threshold))
      
      #Calculate percentages of SP and DP remaining after filtering
      single_pos$remainingSP <- single_pos$SP_after_filter/single_pos$SP_no_filter
      single_pos$remainingDP <- single_pos$DP_after_filter/single_pos$DP_no_filter
      
      correlation_data <- double_pos_subset %>% dplyr::select(c(paste(samplename1),paste(samplename2))) %>% dplyr::filter_if(is.numeric,all_vars(.>0)) #take out genes that are not >0 in both samples
      
      p <- ggplot(double_pos_subset, aes(x=log(get(samplename1)+1,10), y=log(get(samplename2)+1,10), col=colouring)) +
        geom_point(alpha=0.7,size=0.5) +
        #geom_hline(yintercept = log(quantile(single_pos$max,probs = 0.95, na.rm=TRUE)+1,10)) +
        #geom_vline(xintercept = log(quantile(single_pos$max,probs = 0.95, na.rm=TRUE)+1,10)) +
        theme_classic() +
        scale_x_continuous(limits=c(0,5)) +
        scale_y_continuous(limits=c(0,5)) +
        theme(plot.title=element_text(size=9), plot.subtitle=element_text(size=7),
              legend.title=element_blank(), legend.text=element_text(size=8), legend.position="top",
              axis.title=element_text(size=8)) +
        #scale_color_manual(values=c("lightblue","black")) +
        scale_color_manual(values=color_panel[2:6]) +
        labs(title=paste(varname, "with cut-off of", threshold),
             subtitle=paste0("% remaining after filtering: ",
                             round(single_pos[single_pos$full_name==paste(varname),]$remainingDP*100,2), "% of DP (",  
                             round(single_pos[single_pos$full_name==paste(varname),]$remainingSP*100,2),"% of SP)\nSpearman r (genes>0): ",round(cor(correlation_data[,paste(samplename1)], correlation_data[,paste(samplename2)], method="spearman"),2)),
             x=paste0("log(",sample1,"+1,10)"), y=paste0("log(",sample2,"+1,10)")) +
        guides(colour = guide_legend(override.aes = list(alpha = 1)))
      print(p)
    }
  }
}

cutoff_kit <- single_pos %>% dplyr::group_by(GraphKit) %>% dplyr::summarise(median_th = median(filter_threshold))

spread(cutoff_kit, key = GraphKit, value=median_th)

p <- ggplot(single_pos, aes(x=GraphKit, y=filter_threshold, color=Replicates)) +
  geom_point() +
  theme_point +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))+
  scale_color_manual(values=color_panel[3:5]) +
  scale_y_continuous(limits = c(0,NA)) +
  labs(y="cut-off (counts)", x="", title="95% SP elimination cut-off",
       subtitle=paste("overall median cut-off:",median(single_pos$filter_threshold))) +
  geom_hline(yintercept = median(single_pos$filter_threshold), color="grey80")

print(p)

### number of miRNAs based on filtering
#### Original data
miRNAs_long <- miRNAs %>% gather(., key="UniqueID", value="counts", -"MIMATID") %>% #long format
  left_join(., dplyr::select(sample_annotation,c(UniqueID,GraphKit,SampleID)), by="UniqueID") #add kit column

#keep only protein coding miR with more than 0 counts
miRNAs_long <- miRNAs_long %>% filter(counts > 0)


number_miR_detected <- miRNAs_long %>% group_by(SampleID) %>% dplyr::summarize(miR_above0=n()) 
number_miR_detected <- full_join(number_miR_detected, 
                                 miRNAs_long %>% group_by(SampleID) %>%
                                   dplyr::summarize(total_est_counts_above0=sum(counts)), 
                                 by="SampleID")


#cutoff_kit <- single_pos %>% group_by(GraphKit) %>% dplyr::summarise(median_th = median(filter_threshold))

miRNAs_cutoff <- miRNAs %>% gather(., key="UniqueID", value="counts", -"MIMATID") %>% #long format
  left_join(., dplyr::select(sample_annotation,c(UniqueID,GraphKit,SampleID)), by="UniqueID") %>% #add kit column
  left_join(., cutoff_kit, by="GraphKit") #add the median cut-off for each kit

miRNAs_cutoff <- miRNAs_cutoff %>% 
  filter(counts > median_th) 

number_miR_detected <- full_join(number_miR_detected, 
                                 miRNAs_cutoff %>% group_by(SampleID) %>%
                                   dplyr::summarize(miR_aboveTh=n()),
                                 by="SampleID")
number_miR_detected <- full_join(number_miR_detected, 
                                 miRNAs_cutoff %>% group_by(SampleID) %>%
                                   dplyr::summarize(total_est_counts_aboveTh=sum(counts)), 
                                 by="SampleID")

number_miR_detected <- left_join(number_miR_detected, 
                                 dplyr::select(sample_annotation, c(GraphKit, RNAisolation, EluateV,SampleID, TechnicalReplicate)),
                                 by="SampleID")
#convert to the original format
miRNAs_cutoff <- dplyr::select(miRNAs_cutoff, MIMATID, UniqueID, counts) %>% 
  spread(., key=UniqueID, value=counts)

######keep only plots!
p1 <- ggplot(number_miR_detected,aes(x=GraphKit,y=miR_above0, col=RNAisolation))+
  geom_point(size=1.3) +
  theme_point+
  theme(panel.grid.major.x=element_line(linetype = "dashed",color="gray88"), axis.text.x=element_text(size=8)) +
  labs(x="", y="# miRNAs",title="Absolute number of miRNAs (unfiltered)") +
  scale_y_continuous(limits = c(0, NA))+
  scale_color_manual(values=color_panel2) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
#facet_wrap(~GraphKit, nrow = 1)

#p1

p3 <- ggplot(number_miR_detected,aes(x=GraphKit,y=miR_aboveTh, col=RNAisolation))+
  geom_point(size=1.3) +
  theme_point+
  theme(panel.grid.major.x=element_line(linetype = "dashed",color="gray88"), axis.text.x=element_text(size=8)) +
  labs(x="", y="# miRNAs",title="Absolute number of miRNAs (above cutoff)") +
  scale_y_continuous(limits = c(0, NA))+
  scale_color_manual(values=color_panel2) +
  theme(axis.text.x=element_text(angle=90, vjust=0.5, hjust=1))
#facet_wrap(~GraphKit, nrow = 1)

source(here::here("scripts","HelperFunctions.R"))
library(gridExtra)

grid_arrange_shared_legend(p1 + scale_y_continuous(limits=c(0,500)),
                           p3 + scale_y_continuous(limits=c(0,500)))


#pdf(file="/Users/almorlio/Dropbox (Vandesompele lab)/Personal/Projects/exRNAQC011_smallRNA/data_output/Kits_smallRNA_number_extrasubs.pdf", height=4, width=6, useDingbats = F)
p3 + coord_flip() + theme(axis.text.x = element_text(angle=0, hjust=0.5))
#dev.off()

tmp_summary <- number_miR_detected %>% 
  filter(GraphKit != "MAP2") %>% filter(GraphKit != "MAP4") %>% #remove MagnaPure kits
  group_by(GraphKit) %>% dplyr::summarize(score_median=median(miR_aboveTh), score_mean=mean(miR_aboveTh))
tmp <- cbind("GraphKit"=paste(tmp_summary$GraphKit),
             "gene_count_median"=paste(scale(tmp_summary$score_median, center=T, scale=T)),
             "gene_count_mean"=paste(scale(tmp_summary$score_mean, center=T, scale=T))) %>% 
  as.tibble(.) %>% 
  type_convert(.) #convert columns with all numbers to type "dbl"

z_score_df <- full_join(z_score_df, tmp, by="GraphKit")

data_summary <- full_join(data_summary, number_miR_detected %>% 
                            filter(GraphKit != "MAP2") %>% filter(GraphKit != "MAP4") %>%
                            dplyr::select(c(SampleID,GraphKit,miR_count_aboveTh = miR_aboveTh)),
                          by="SampleID")
number_miR_detected2 <- number_miR_detected %>% filter(GraphKit != "MAP2") %>% filter(GraphKit != "MAP4")
z_score_indiv <- cbind("SampleID"=paste(number_miR_detected2$SampleID),
                       "GraphKit"=paste(number_miR_detected2$GraphKit),
                       "gene_count"=paste(scale(number_miR_detected2$miR_aboveTh, center=T, scale=T))) %>% #z-score transformation
  as.tibble(.) %>% 
  type_convert(.) #convert columns with all numbers to type "dbl"
#variable correlation
#cor_z_scores <- cor(z_score_df %>% select(-"GraphKit"), method = "spearman")

#variable correlation
#cor_z_scores <- cor(z_score_df %>% select(-"GraphKit"), method = "spearman")
robust_z <- robzscore(number_miR_detected2, "miR_aboveTh")
#cor(robust_z, z_score_indiv$gene_count, method="spearman")
#test <- cbind(z_score_indiv, rob_z=robust_z)
z_score_indiv_robust <- cbind("SampleID"=paste(number_miR_detected2$SampleID),
                              "GraphKit"=paste(number_miR_detected2$GraphKit),
                              "gene_count"=robust_z) %>%
  as.tibble(.) %>% 
  type_convert(.) #convert columns with all numbers to type "dbl"


rm(list=grep("tmp|test|filter",ls(),value=T))

### ALC
## ALC{.tabset}
# - AlC (area-left-of-curve) calculation:
#   - Only take genes into account that reach threshold (95% SP elimination)
# - Replace all zero counts by NA
# - Determine log2 ratios of counts for each gene (each time between 2 replicates of the kit)
# - Take the absolute value of these log2 ratios and rank them (percentage rank) -> plot this
# - The faster the curve reaches 100% -> the smaller the ALC -> the better (indicates that the replicates give similar counts for each gene)
# - Summarizing ALC value 
# - = the average of the ALC values (after removing NAs)
# - Dot plot at the end gives an overview of these values 
# - **Score: lower ALC = better**
  
### All - At least 1 replicate > threshold (and other > 0)

double_positives<-miRNAs %>% replace(., is.na(.),0) %>% #first change NAs to 0 
  mutate_if(is.numeric, funs(round)) #round everything to nearest integer

double_positives<-double_positives[rowSums(double_positives %>% select_if(is.numeric) > 1, na.rm=T)>0,] # keep only the rows where at least one sample has more than 1

ALC_input_all <- data.frame()
for (duplicate_type in unique(sample_annotation$GraphKit)){
  sample_duplicates<-sample_annotation%>% 
    #filter(TechnicalReplicate != "RNA2") %>% #leave out bad pippinprep replicates
    filter(GraphKit==duplicate_type) %>%
    dplyr::select(UniqueID,GraphKit,TechnicalReplicate)
  cutoff95SP <- cutoff_kit[cutoff_kit$GraphKit==duplicate_type,]$median_th
  if(nrow(sample_duplicates)>1){
    #print(duplicate_type)
    #double_positives_sample<-double_positives %>% dplyr::select(MIMATID,sample_duplicates$UniqueID) # only keep the replicates of one type
    
    samples_comb <- combn(sample_duplicates$UniqueID,2) #compare 2 of the 3 samples at a time
    for (n_col in 1:ncol(samples_comb)) {
      #print(samples_comb[,n_col])
      nr_runA <- gsub("^[A-Z]+","",sample_annotation[sample_annotation$UniqueID==
                                                       samples_comb[1,n_col],]$TechnicalReplicate)
      nr_runB <- gsub("^[A-Z]+","",sample_annotation[sample_annotation$UniqueID==
                                                       samples_comb[2,n_col],]$TechnicalReplicate)
      varname <- paste0("Rep",nr_runA,"_",nr_runB) #make a name so you can backtrace which replicates are compared
      #print(varname)
      double_positives_2repl <- double_positives %>% dplyr::select(MIMATID, paste(samples_comb[1,n_col]), paste0(samples_comb[2,n_col])) %>%
        filter_if(is.numeric, all_vars(.>0)) %>% #keep only the 2 replicates of one kit and remove genes where not at least 1 of them has > 0 counts
        filter_if(is.numeric, any_vars(.>cutoff95SP)) # remove genes where neither of the replicates reach the threshold that removes 95% of SP in that kit
      correlation_samples<-double_positives_2repl %>%
        mutate(log2_ratio=abs(log(get(samples_comb[1,n_col]),2)-log(get(samples_comb[2,n_col]),2))) %>%
        dplyr::select(log2_ratio,MIMATID) #%>% drop_na()
      ALC_input<-correlation_samples %>% arrange(log2_ratio) %>% # order by log2 ratio and then make a rank (counter) and rescale this to 1 (perc_counter)
        #mutate(rank=percent_rank(log2_ratio)) %>% # this does not work: gives everything with log2ratio = 0 rank 0
        mutate(counter = seq(1:nrow(double_positives_2repl)), GraphKit=duplicate_type, Replicates=varname) %>%
        mutate(ReplID = paste0(GraphKit,"-",Replicates), perc_counter = counter/nrow(double_positives_2repl))
      
      ALC_input_all <- rbind(ALC_input_all, ALC_input)
    }
  }
}

max_ALC <-max(ALC_input_all$log2_ratio) # calculate the max ALC over everything (necessary for area calculation -> should always be the same in order to compare among kits)
ALC <- ALC_input_all %>% group_by(GraphKit,Replicates) %>% 
  dplyr::summarise(ALC_calc = sum(log2_ratio)/(max_ALC*length(MIMATID))) %>%
  #dplyr::summarise(ALC_calc = sum(log2_ratio)/(max_ALC)) %>%
  mutate(ReplID = paste0(GraphKit,"-",Replicates))

# for (replicates in unique(ALC_input_all$ReplID)) { # plot the ALC (= the colored part of the plot)
#   print(ggplot(ALC_input_all %>% dplyr::filter(ReplID==replicates) %>% 
#                  mutate(log2_ratio_resc = log2_ratio/max_ALC))+
#           geom_line(aes(x=log2_ratio_resc,y=perc_counter))+
#           #facet_wrap(~ReplID) +
#           geom_ribbon(aes(x=log2_ratio_resc,ymin=perc_counter,ymax=1), fill=color_panel1[gsub("-.*$","",replicates)])+
#           geom_hline(aes(yintercept = 1))+
#           theme_classic()+
#           scale_x_continuous(limits=c(0,1)) + 
#           scale_y_continuous(expand = c(0, 0)) +
#           theme(legend.position = "none") +
#           labs(title=paste(replicates),
#                subtitle=paste("ALC:", round(ALC[ALC$ReplID==replicates,]$ALC_calc,2)), #print ALC for this particular comparison!
#                y="rank percentile",x="rescaled log2 ratio"))
#   
# }
#ALC<-sum(ALC_input$log2_ratio)/(length(ALC_input$log2_ratio))

ALC_melt <- left_join(ALC, sample_annotation[,c("GraphKit","RNAisolation")], by="GraphKit") %>% unique()
ALCplot <- ggplot(ALC_melt %>% drop_na())+
  geom_point(aes(y=ALC_calc,x=GraphKit,color=RNAisolation),size=1.3)+
  #geom_text_repel(aes(y=value,x=GraphKit), nudge_x=0.1)+
  theme_point +
  labs(y="ALC",title="Pairwise ALCs (both > 0, at least 1 > cut-off)", subtitle="lower ALC = better")+
  #  theme(panel.grid.major.y=element_line(linetype = "dashed",color="lightgray"))+
  scale_color_manual(values=color_panel2) +
  scale_y_continuous(limits=c(0,NA))

#pdf(file="/Users/almorlio/Dropbox (Vandesompele lab)/Personal/Projects/exRNAQC011_smallRNA/data_output/Kits_smallRNA_ALC_extrasubs.pdf", height=4, width=6, useDingbats = F)
ALCplot + coord_flip()
#dev.off()
ALCplot + theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1))


tmp_summary <- ALC_melt %>% 
  filter(GraphKit != "MAP2") %>% filter(GraphKit != "MAP4") %>% #remove MagnaPure kits
  group_by(GraphKit) %>% dplyr::summarize(score_median=median(ALC_calc), score_mean=mean(ALC_calc)) %>% 
  mutate_if(is.numeric, funs(. * -1)) #make counts negative so that a higher metric value corresponds to a better performance

#calculate z score for the median value of every kit (already negative value!)
tmp <- cbind("GraphKit"=paste(tmp_summary$GraphKit),
             "reproducibility_ALCgenesSPno0_median"=paste(scale(tmp_summary$score_median, center=T, scale=T)),
             "reproducibility_ALCgenesSPno0_mean"=paste(scale(tmp_summary$score_mean, center=T, scale=T))) %>% 
  as.tibble(.) %>% 
  type_convert(.) #convert columns with all numbers to type "dbl"

#join the z values to the data frame
z_score_df <- full_join(z_score_df, tmp, by="GraphKit")

# add the column with the median ALC values to the data_summary df (reconvert to original sign by multiplying with -1 again)
data_summary <- full_join(data_summary, tmp_summary %>% 
                            filter(GraphKit != "MAP2") %>% filter(GraphKit != "MAP4") %>%
                            dplyr::select(c(GraphKit,ALCgenes_withSPno0_median = score_median, ALCgenes_withSPno0_mean = score_mean)),
                          by="GraphKit") %>% mutate(ALCgenes_withSPno0_median = ALCgenes_withSPno0_median*(-1), ALCgenes_withSPno0_mean = ALCgenes_withSPno0_mean*(-1))

# Make individual z-scores!
ALC_melt2 <- ALC_melt %>% mutate(counter = rep(seq(1,3),length(unique(GraphKit)))) %>%
  mutate(SampleID=paste0(GraphKit,"_RNA",counter))

ALC_melt3 <- ALC_melt2 %>% filter(GraphKit != "MAP2") %>% filter(GraphKit != "MAP4")

tmp <- cbind("SampleID"=paste(ALC_melt3$SampleID),
             "GraphKit"=paste(ALC_melt3$GraphKit),
             "ALC_SPno0"=paste(scale(ALC_melt3$ALC_calc*(-1), center=T, scale=T))) %>% #z-score transformation
  as.tibble(.) %>% 
  type_convert(.) #convert columns with all numbers to type "dbl"

z_score_indiv <- left_join(z_score_indiv, tmp, by=c("SampleID","GraphKit"))
#variable correlation
#cor_z_scores <- cor(z_score_df %>% select(-"GraphKit"), method = "spearman")


# Calculate robust z scores!
robust_z <- robzscore(ALC_melt3 %>% mutate(ALC_calc_neg=ALC_calc*(-1)), "ALC_calc_neg")
tmp <- cbind(GraphKit = paste(ALC_melt3$GraphKit), SampleID=paste(ALC_melt3$SampleID), ALC_SPno0=robust_z) %>% as.tibble(.) %>% type_convert(.)
#cor(robust_z, z_score_indiv$gene_count, method="spearman")
#test <- cbind(z_score_indiv, rob_z=robust_z)
z_score_indiv_robust <- left_join(z_score_indiv_robust, tmp, by=c("SampleID","GraphKit"))

rm(list=grep("tmp|test|filter|melt",ls(),value=T))



### Plot miR count against ALC to select kit

ggplot(z_score_indiv, aes(x=gene_count, y=ALC_SPno0, col=GraphKit)) + geom_point() + 
  theme_point + ggrepel::geom_text_repel(aes(label=GraphKit))  + 
  scale_colour_manual(values=color_panel1) + theme(legend.position="none")

z_score_indiv_med <- z_score_indiv %>% group_by(GraphKit) %>% summarise_if(is.numeric, median, na.rm=T)
#pdf("data_output/selection_kits_miRNA_regularz.pdf", height=5, width=6)
ggplot(z_score_indiv_med, aes(x=gene_count, y=ALC_SPno0, col=GraphKit)) + 
  geom_point() + 
  theme_point + 
  ggrepel::geom_text_repel(aes(label=GraphKit)) + scale_colour_manual(values=color_panel1) + 
  theme(legend.position="none") + 
  labs(x="sensitivity (miRNA count)", y="reproducibility (ALC)", subtitle="Regular z-scores (median per kit)")
#dev.off()


z_score_indiv_robust_med <- z_score_indiv_robust %>% group_by(GraphKit) %>% summarise_if(is.numeric, median, na.rm=T)
#pdf("data_output/selection_kits_miRNA_robustz.pdf", height=5, width=6)
ggplot(z_score_indiv_robust_med, aes(x=gene_count, y=ALC_SPno0, col=GraphKit)) + 
  geom_point() + 
  theme_point + 
  ggrepel::geom_text_repel(aes(label=GraphKit)) + 
  scale_colour_manual(values=color_panel1) + 
  theme(legend.position="none") + labs(x="sensitivity (miRNA count)", y="reproducibility (ALC)", subtitle="Robust z-scores (median per kit)")
