# Initialise variables ===================================================
cat("Initialising variables\n")
args			  <- commandArgs(trailingOnly = TRUE)
sample_name		  <- args[2]
path	  <- args[1]


# Define names of in- and output files =============================================
cat("Defining names of in- and output files\n")
file_read_length <- paste(path,"/"	, sample_name, "_read_length_new.txt"	, sep = "")
file_mirs		 <- paste(path,"/"	, sample_name, "_isomiRs.txt"		, sep = "")
file_contam		 <- paste(path,"/"	, sample_name, "_contam.txt"		, sep = "")
file_noannot	 <- paste(path,"/"	, sample_name, "_not_annotated.txt" , sep = "")
plot_name		 <- paste(path,"/"	, sample_name, ".pdf"				, sep = "")

# Check if input files are empty =========================================
cat("Checking if input files are empty\n")
if (file.info(file_read_length)$size == 0) {
	stop(paste0("File is empty: ", file_read_length))
} else{data_read_length <- read.table(file_read_length, header = T	, sep = "\t")}
if (file.info(file_mirs)$size == 0) {
	stop(paste0("File is empty: ", file_mirs))
} else {data_mirs		 <- read.table(file_mirs					, sep = "\t")}
if (file.info(file_contam)$size == 0) {
	stop(paste0("File is empty: ", file_contam))
} else {data_contam		 <- read.table(file_contam					, sep = "\t")}
if (file.info(file_noannot)$size == 0) {
	stop(paste0("File is empty: ", file_noannot))
}else{data_noannot	 <- read.table(file_noannot				 	, sep = "\t")}

# # Read input files =======================================================
# cat("Reading input files\n")
# data_read_length <- read.table(file_read_length, header = T	, sep = "\t")
# data_mirs		 <- read.table(file_mirs					, sep = "\t")
# data_contam		 <- read.table(file_contam					, sep = "\t")
# data_noannot	 <- read.table(file_noannot				 	, sep = "\t")

if(ncol(data_read_length) != 3){
	stop(paste0("Incorrect number of columns found in file ('", file_read_length, "')! Found ", ncol(data_read_length), ", should be 3!"))
}
if(ncol(data_mirs)		  != 7){
	stop(paste0("Incorrect number of columns found in file ('", file_mirs		, "')! Found ", ncol(data_mirs)		  , ", should be 7!"))
}
if(ncol(data_contam)	  != 5){
	stop(paste0("Incorrect number of columns found in file ('", file_contam		, "')! Found ", ncol(data_contam)	  , ", should be 5!"))
}
if(ncol(data_noannot)	  != 7){
	stop(paste0("Incorrect number of columns found in file ('", file_noannot	, "')! Found ", ncol(data_noannot)	  , ", should be 7!"))
}
rownames(data_mirs)	  <- as.character(data_mirs[,2])
rownames(data_contam) <- as.character(data_contam[,2])

data_noannot <- data_noannot[!data_noannot[,1]%in%rownames(data_mirs),]
data_noannot <- data_noannot[!data_noannot[,1]%in%rownames(data_contam),]

# Fill array with data for plot ==========================================
cat("Filling array with data for plot\n")
total_mir		<- c()
total_tRNA		<- c()
total_piRNA		<- c()
total_rRNA		<- c()
total_snRNA		<- c()
total_miscRNA	<- c()
total_multiple	<- c()
total_no_annot	<- c()

for(i in seq(15,40,1)){
	tmp				<- data_read_length[data_read_length$read_length == i,]
	tmp_mir			<- tmp			[  tmp[,1]			 %in% rownames(data_mirs)			,]
	tmp_no_mir		<- tmp			[! tmp[,1]			 %in% rownames(data_mirs)			,]
	tmp_contam		<- tmp_no_mir	[  tmp_no_mir[,1]	 %in% rownames(data_contam)			,]
	tmp_no_contam 	<- tmp_no_mir	[! tmp_no_mir[,1]	 %in% rownames(data_contam)			,]
	tmp_noannot		<- tmp_no_contam[  tmp_no_contam[,1] %in% as.character(data_noannot[,1]),]
#	tmp_nomap		<- tmp_no_contam[! tmp_no_contam[,1] %in% as.character(data_noannot[,1]),]
	
	type			<- as.character(data_contam[rownames(tmp_contam),4])
	tmp_contam		<- cbind(tmp_contam, type)
	
	total_mir		<- c(total_mir,		 sum(tmp_mir[,3]))
	total_tRNA		<- c(total_tRNA,	 sum(tmp_contam[tmp_contam[,"type"] %in% c("Mt_tRNA", "tRNA"),3]))
	total_rRNA		<- c(total_rRNA,	 sum(tmp_contam[tmp_contam[,"type"] %in% c("rRNA")			 ,3]))
	total_piRNA		<- c(total_piRNA,	 sum(tmp_contam[tmp_contam[,"type"] %in% c("piRNA")			 ,3]))
	total_snRNA		<- c(total_snRNA,	 sum(tmp_contam[tmp_contam[,"type"] %in% c("snoRNA", "snRNA"),3]))
	total_miscRNA	<- c(total_miscRNA,	 sum(tmp_contam[tmp_contam[,"type"] %in% c("misc_RNA")		 ,3]))
	total_multiple	<- c(total_multiple, sum(tmp_contam[tmp_contam[,"type"] %in% c("multiple")		 ,3]))
	total_no_annot	<- c(total_no_annot, sum(tmp_noannot[,3]))
}
total <- rbind(total_mir, total_tRNA, total_rRNA, total_piRNA, total_snRNA, total_miscRNA, total_multiple, total_no_annot)

# Create plot ============================================================
cat("Creating plot\n")
pdf(plot_name)
#par(mar = c(5,6,4,9), xpd=TRUE)

bar_colors = c("#96C123", "#E61D84", "#000098", "#7FCDE9", "#FFD700", "#FFA500", "#808080", "#000000", "white")
gene_types = c("miRNA", "tRNA *", "rRNA *", "piRNA", "sn(o)RNA *", "miscRNA", "multiple", "none", "* (processed) fragments")

barplot(total, names = seq(15,40,1), beside = F,las = 1, ylab = "", border = "white", 
		col = bar_colors, xlab = "Read length", cex.names = 0.7, main = sample_name)

tmp_legend = legend("topright", legend = gene_types, fill = bar_colors, border = "white", bty = "n",
					title = "Annotation", title.col = "black", cex = 1, text.col = 'white')
#text(tmp_legend$text$x, tmp_legend$text$y, lab = gene_types, cex = c(1, 1, 1, 1, 1, 1, 1, 1, 0.7), pos = 4, offset = c(0, 0))
text(tmp_legend$text$x, tmp_legend$text$y, lab = gene_types, pos=4)

title(ylab = "Read count", line = 4)
png_create <- dev.off()

if(length(warnings()) != 0){
	warnings()
}

cat("\nScript finished\n")
