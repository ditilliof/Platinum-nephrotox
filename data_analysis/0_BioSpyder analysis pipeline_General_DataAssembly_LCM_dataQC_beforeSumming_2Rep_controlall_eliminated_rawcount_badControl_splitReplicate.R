

## Data analysis pipeline BioSpyder datasets ##

# 0 - Data assembly, normalization, QCs

# Marije Niemeijer

# Version 5

## Steps:
# Loading of data
# Library size calculation and filtering
# Exluding specific samples for further analysis (NoCells, and others if necessary)
# Normalization of raw counts using DESEq2 package
# Quality control plots (library size and count distribution, replicate correlation, variance across means)
# DEGs determination: Log2FC and p-value calculation using DESEq2 package
# Filtering of DEGs (based on p-value/Log2FC)
# Preparation of txt files for IPA as input files
# Assessing gene expression for specific genes
# Principal component analysis (PCA plots, top25 genes in PCs, contribution of variables)
# Hierarchical clustering of samples and genes in heatmap
# Saving data

#############################################

# NOTES

# SOP for BioSpyder sample preparation and analysis: TOX_021_RNA-based_BioSpyder transcriptomics.docx
#     @ "J:\Workgroups\FWN\LACDR\TOX\1_SOP's\2016 SOP's team\Overview TOX SOP's\2) RNA-based\2) RNA-based_cDNA synthesis"

# Please consider carefully which samples you want to normalize together using Deseq2. 
#     For instance, normalizing multiple cell types at once or separatly might influence the normalization. 
#     When having low amount of replicates per group, analysis together is preferred to increase power.

# When you want to compare different samples to different controls, either: 
#     run script separately on subsets of data (not possible to compare normalized counts between separate analysis, only at FC level)
# or  normalize together and retrieve specific comparisons (comparison both at normalized counts and FC level possible)

# By default Log2FCs using deseq2 are calculated using ControlSample (defined at input user)

#############################################

rm(list=ls(all=TRUE))


###### INPUT USER ######
NoCells <- c()

MetaName <- "meta_data_BCL.txt"      #use template and save as txt file
CountThreshold <- 10                              #Exclude samples lower than this number for library size
VariableQCs <- "meanID"                                #Variable to check distribution Library size and counts (choose meanID, CELL_ID, TREATMENT, SAMPLE_ID, EXP_ID, TIMEPOINT or REPLICATE)
RemoveSamples <- NoCells                             #Give list of SAMPLE_IDs for samples which you would like to remove for further analysis (for example MAQC, NoCells samples, outliers)
ControlSample <- c('HepG2_HepG2_DMSO_conc0.2_TP8')        #Sample name (meanID, same for replicates) of control for log2FC and p-value calculations
                                                          #meanID: 'EXP_ID'_'CELL_ID'_'TREATMENT'_conc'CONCENTRATION'_TP'TIMEPOINT

NormMethod <- "CPM"                                    #Fill in desired normalization method prior to Deseq2 DEG calculation: "CPM" or "DESeq2"
Threshold_padj <- 0.1                                     #Filtering of genes with p-adj lower than threshold
Threshold_log2FC <- 0.1                                   #Filtering of genes with absolute log2FC higher than threshold (genes lower than -0.5, and genes higher than 0.5)
Filtering <- "padj"                                       #Fill in filtering method ("padj", "log2FC" or "padj&log2FC")
CheckGenes <- NA #Check expression of individual genes among samples (Give gene symbols)
                                                          #You can list as many genes as you like, however graph may become too full
Geneset <- NA                 #Fill in file name of specific gene list of interest (GeneSymbols, txt format). If not used, fill in NA.
                                                          

# Working directories #

CountTableDir <- "./data/in vivo" #Folder containing only raw count table files (csv or excel files)
inputDir = CountTableDir
outputDir_GR <- CountTableDir

scriptDir <- "./data_analysis/"

source(file.path(scriptDir, "A_BioSpyder analysis pipeline_Lib.R"))

## Changes Tracking of versions ##

#V2
#Adjusted comparison samples in metadata and count table
#Changed DEGs calculation method (multiple pair-wise comparison)
#Addition of FC vs p-value/basemean graphs
#Addition of loading vsn Rpackage
#Saving IPA txt input files with no quotes around character strings
#Including error message when making a heatmap with more than 4000 genes

#v3
#Adding extra loading of libraries (vsn, org.Hs.eg.db, AnnotationDbi)
#Changed function PathwayProbesCounts, ability to match gene list based on entrez ID or probe ID to subset count table
#Addition of Entrez gene ID in output df_DEGs table
#Improved description of functions and arguments
#Possible to fill in NA at Geneset at input user, if not interested in particular gene set
#Message with how many genes are taken along in size factor calculation for deseq2 normalization

#v4
#CPM and Deseq2 normalization
#Deseq2 DEG calculation based on CPM or Deseq2 normalized counts
#Adding read depth to meta dataframe and returning summary
#Library size filtering at 500.000 reads as default
#Fixed bug with replicate correlation plot

#v5
#Fixed bug with replicate correlation plot when having different number of replicates
#Add plot to evaluate log2 norm counts for specific genes for each replicate/sample

## Assembly of Data ##

######  Load data & library size filter ###### 
meta <- read.delim(paste0(inputDir, MetaName), stringsAsFactors = FALSE)
meta[meta$CONCENTRATION == 0,]$TIMEPOINT <- "control"
rownames(meta) <- meta[,1]
count_1 <- LoadCountData(CountTableDir)
count1 <- read_csv("./data/in vivo/Count_per_gene_per_sample_raw_BCL-SP0195.csv")
colnames(count_1) <- gsub("S_2003_repeat", "S_2005", colnames(count_1))
colnames(count_1) <- gsub("S_16002", "S_16005", colnames(count_1))
#count$probe <- rownames(count)
count_2 <- read.csv("./data/in vivo/Counts_per_gene_per_sample_raw_BCL-SP0208.csv")

count <- merge(count_1, count_2, by=0, all=TRUE) 
rownames(count) <- count$Row.names
count$Row.names <- NULL

meta$meanID <- paste0(meta$EXP_ID, "_",
                      meta$LOCATION_ID, "_",
                      meta$TREATMENT, "_conc",
                      meta$CONCENTRATION,"_TP", 
                      meta$TIMEPOINT)

meta$SAMPLENAME <- paste0(meta$EXP_ID, "_",
                          meta$LOCATION_ID, "_",
                          meta$TREATMENT, "_TP", 
                          meta$TIMEPOINT, "_conc",
                          meta$CONCENTRATION, "_rep",
                          meta$REPLICATE)

if(!all(rownames(meta) %in% colnames(count)) |
   !all(colnames(count) %in% rownames(meta))){
  warning("No identical sample names (meta & counts)")
  print(paste0("Unmatched samples in count table: ", paste0(c(colnames(count)[which(!colnames(count) %in% rownames(meta))]), collapse = ", ")))
  print(paste0("Unmatched samples in meta: ", paste0(c(rownames(meta)[which(!rownames(meta) %in% colnames(count))]), collapse = ", ")))
} 

count <- count[, rownames(meta)]

## Library size filter ##

meta$LIB_SIZE <- colSums(count, na.rm = TRUE)
meta$ReadDepth <- meta$LIB_SIZE / nrow(count)
summary(meta$PERCENTAGE_MAPPED)
summary(meta$ReadDepth)
meta_LowLibSize <- meta[which(meta$PERCENTAGE_MAPPED < CountThreshold), ]
CELLID_LowLibSize <- as.data.frame(table(meta_LowLibSize $CELL_ID))
meta_Filtered <- meta[which(meta$PERCENTAGE_MAPPED > CountThreshold), ] 
count_Filtered <- count[, rownames(meta_Filtered)]

## Check if CountTable contains NAs ##

count_NA <- count_Filtered[rowSums(is.na(count_Filtered)) > 0, colSums(is.na(count_Filtered)) > 0]
if(nrow(count_NA) == 0){
  print("No NAs in counts table")
} else {
  warning("NAs present in counts table")
  print(table(count_NA)); print(dim(count_NA))
}

count_Filtered <- na.omit(count_Filtered)

## Excluding specific samples for further analysis

#meta_Filtered <- meta_Filtered[which(!grepl(paste0(RemoveSamples, collapse = "|"), meta_Filtered$SAMPLE_ID)),]

count_Filtered <- count_Filtered[, rownames(meta_Filtered)]

## Number of genes with at least one zero count among samples

NrGenes_zero <- sum(apply(count_Filtered, MARGIN = 1, function(x) sum(x == 0)) > 0)
message(paste0(NrGenes_zero, "/", nrow(count_Filtered), " genes have at least one zero count among samples and will not be used for size factor calculation"))

##correlation check rawcount###
require(PerformanceAnalytics)
require(Hmisc)


meta_filetered_whole <- meta_Filtered[meta_Filtered$LOCATION_ID == "WHOLE",]

meta_filetered_CPT <- meta_Filtered[meta_Filtered$LOCATION_ID == "CPT",]

meta_filetered_PPT <- meta[meta$LOCATION_ID == "PPT",]

#meta_filetered_GM <- meta_Filtered[meta_Filtered$LOCATION_ID == "GM",]


count_Filtered_whole <- count_Filtered[, rownames(meta_filetered_whole)]

for (i in 1 : length(unique(meta_filetered_whole$TIMEPOINT))){
  CON <- unique(meta_filetered_whole$TIMEPOINT)[i]
  
  meta_filtered_sub <- meta_filetered_whole[meta_filetered_whole$TIMEPOINT == CON, ]
  count_filtered_sub <- count_Filtered_whole[, rownames(meta_filtered_sub)]
  
  res2 <- rcorr(as.matrix(count_filtered_sub))
  rrSquare <- as.data.frame(res2[["r"]])
  
  colBreaks <- seq(-1,1, length.out = 100)
  
  my.colors  <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  pdf(paste0(file = "//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2021/LCM_Part2/graph/correlation_matrix_replicate_",CON,"_whole_.pdf"), height =15, width = 15)
  print(pheatmap(rrSquare, breaks = colBreaks, color = my.colors))
  dev.off()
}


count_Filtered_CPT <- count_Filtered[, rownames(meta_filetered_CPT)]

for (i in 1 : length(unique(meta_filetered_CPT$TIMEPOINT))){
  CON <- unique(meta_filetered_CPT$TIMEPOINT)[i]
  
  meta_filtered_sub <- meta_filetered_CPT[meta_filetered_CPT$TIMEPOINT == CON, ]
  count_filtered_sub <- count_Filtered_CPT[, rownames(meta_filtered_sub)]
  
  res2 <- rcorr(as.matrix(count_filtered_sub))
  rrSquare <- as.data.frame(res2[["r"]])
  
  colBreaks <- seq(-1,1, length.out = 100)
  
  my.colors  <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  pdf(paste0(file = "//VUW/PerSonal$/HomeS/W/wijayalS/DeSktop/ExperimentS/2021/LCM_Part2/graph/correlation_matrix_replicate_",CON,"_CPT_.pdf"), height =15, width = 15)
  print(pheatmap(rrSquare, breaks = colBreaks, color = my.colors))
  dev.off()
}

count_Filtered_PPT <- count[, rownames(meta_filetered_PPT)]
#meta_filetered_PPT <- meta_filetered_PPT[!meta_filetered_PPT$TIMEPOINT == "288",]

Low_corr_PPT <- NULL
for (i in 1 : length(unique(meta_filetered_PPT$TIMEPOINT))){
  CON <- unique(meta_filetered_PPT$TIMEPOINT)[8]
  
  meta_filtered_sub <- meta_filetered_PPT[meta_filetered_PPT$TIMEPOINT == CON, ]
  count_filtered_sub <- count_Filtered_PPT[, rownames(meta_filtered_sub)]
  
  res2 <- rcorr(as.matrix(count_filtered_sub))
  rrSquare <- as.data.frame(res2[["r"]])
  threshold <- 0.9
  rrSquare_low <- rrSquare
  diag(rrSquare_low) <- 0
  ok <- as.data.frame(apply(abs(rrSquare_low) > threshold, 1, any))
  Low_corr_PPT <- rbind(Low_corr_PPT, ok)
  my.colors  <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)
  
  pdf(paste0(file = "//VUW/PerSonal$/HomeS/W/wijayalS/DeSktop/ExperimentS/2021/LCM_Part2/graph/correlation_matrix_replicate_",CON,"_PPT_.pdf"), height =15, width = 15)
  print(pheatmap(rrSquare, breaks = colBreaks, color = my.colors))
  dev.off()
}
###diStribution###

count_Filtered_LCM <- count_Filtered[, !grepl("Slide", colnameS(count_Filtered))]
meta_Filtered_LCM <- meta_Filtered[!grepl("Slide", meta_Filtered$SAMPLE_ID),]
meta_Filtered_clean <- meta_Filtered[meta_Filtered$SAMPLE_ID %in% colnames(Count_Filtered_Clean),]

CPM <- apply(Count_Filtered_Clean, 2, function(x) (x/sum(x))*1000000)

dds <- DESeqDataSetFromMatrix(countData = Count_Filtered_Clean,
                              colData = meta_Filtered_clean,
                              design = ~ meanID)

if (NormMethod == "DESeq2"){ 
  dds <- estimateSizeFactors(dds) 
} else if (NormMethod == "CPM"){
  sizeFactors(dds) <- colSums(Count_Filtered_Clean)/1000000
} else { 
  warning("Incorrect input for 'NormMethod', please double check above!")
}

Norm <- as.data.frame(counts(dds, normalized = TRUE)) #Normalized by Size factorS
log2Norm <- log2(Norm + 1) #Log2 normalization of Size factor corrected countS


# CountS & normalized countS diStribution amoung CELL IDS #
count_long <- melt(as.matrix(Count_Filtered_Clean))
colnames(count_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
count_long <- left_join(count_long, meta_Filtered[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

log2Norm_long <- melt(as.matrix(log2Norm))
colnames(log2Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
log2Norm_long <- left_join(log2Norm_long, meta_Filtered[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])



# Counts & normalized counts distribution amoung CELL IDs #
for(i in 1:length(unique(count_long$CONCENTRATION))){
  DOSE <- unique(count_long$CONCENTRATION)[i] 
  count_long_sub <- count_long[count_long$CONCENTRATION == DOSE,]
  
  count_long_sub$VariableQCs = count_long_sub[,which(names(count_long_sub)==VariableQCs)]
  pdf(paste0(outputDir_GR, "Distribution-Counts_aftercorrCleanup_version2_", VariableQCs,"_", DOSE,".pdf"), width = 6, height = 4)
  print(ggplot(count_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS+1)) +
          geom_boxplot(size = 0.3, outlier.size = 0.5) +
          scale_y_log10(limits = c(1, max(count_long_sub$COUNTS))) +
          theme_classic() +
          theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
                axis.title.x = element_text(size=12, vjust = 0.25),
                axis.title.y = element_text(size=12, vjust = 1),
                axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
                axis.text.y = element_text(size=12)) +
          ggtitle("Distribution counts") + ylab('counts') + xlab(VariableQCs))
  dev.off()
  
}

for(i in 1:length(unique(log2Norm_long$CONCENTRATION))){
  DOSE <- unique(log2Norm_long$CONCENTRATION)[i] 
  log2Norm_long_sub <- log2Norm_long[log2Norm_long$CONCENTRATION == DOSE,]
  
  #for(u in 1:length(unique(log2Norm_long_sub$LOCATION_ID))){
  #LOC <- unique(log2Norm_long_sub$LOCATION_ID)[u]
  
  #log2Norm_long_sub_2 <- log2Norm_long_sub[log2Norm_long_sub$LOCATION_ID == LOC,]
  
  log2Norm_long_sub$VariableQCs = log2Norm_long_sub[,which(names(log2Norm_long_sub)==VariableQCs)] 
  pdf(paste0(outputDir_GR, "Distribution-CountsNorm_AfterCorrCleanup_version2_", DOSE,".pdf"), width = 6, height = 4)
  
  #print(plot(hist(log2Norm_long_sub_2$COUNTS, col = "green", breaks = 50)))
  #dev.off()
  
  #}
  #}
  
  print(ggplot(log2Norm_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS)) +
          geom_boxplot(size = 0.3, outlier.size = 0.5) +
          theme_classic() +
          theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
                axis.title.x = element_text(size=12, vjust = 0.25),
                axis.title.y = element_text(size=12, vjust = 1),
                axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
                axis.text.y = element_text(size=12)) +
          ggtitle("Distribution Normalized counts") + ylab('Normalized counts') + xlab(VariableQCs))
  dev.off()
  
}

##Eliminating the outlierS baSed on the correlation analySiS##

#Eliminated_samples <- c("S_4007_3_CPT", "S_5001_3_CPT", "S_2003_4_CPT", "S_2005_4_CPT", "S_1005_2_GM", "S_2005_1_GM", "S_5003_1_GM", "S_3003_4_GM",
 #                   "S_4003_1_GM", "S_5003_3_GM", "S_5001_4_PPT", "S_4005_1_PPT", "S_4001_3_PPT", "S_5003_2_PPT", "S_2001_1_PPT", "S_2007_3_PPT")

Eliminated_samples <- c("S_3003_PPT_1")

Count_Filtered_Clean <- t(count_Filtered)
Count_Filtered_Clean <- Count_Filtered_Clean[!rownames(Count_Filtered_Clean) %in% Eliminated_samples,]
Count_Filtered_Clean <- t(Count_Filtered_Clean)
Count_Filtered_Clean <- as.data.frame(Count_Filtered_Clean)



count_long <- melt(as.matrix(Count_Filtered_Clean))
colnames(count_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
count_long <- left_join(count_long, meta[,c("SAMPLE_ID", "SAMPLENAME", "CELL_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

count_long$sampleID <- count_long$SAMPLE_ID
count_long$sampleID <-  sub("(_[^_]+)_.*", "\\1", count_long$sampleID)
count_long$sampleID <- paste0(count_long$sampleID,"_", count_long$LOCATION_ID, sep = "")
count_long$sampleID <- gsub("Whole", "slide", count_long$sampleID)
count_long$SAMPLE_ID <- NULL
count_long$SAMPLENAME <- NULL
count_long <- aggregate(COUNTS ~ Probe_ID + sampleID + EXP_ID + TIMEPOINT + TREATMENT + CONCENTRATION + REPLICATE + meanID + LOCATION_ID , data = count_long, sum)
count_Filtered <- count_long[,c("Probe_ID", "sampleID", "COUNTS")]
count_Filtered <- spread(count_Filtered,sampleID, COUNTS)
rownames(count_Filtered) <- count_Filtered$Probe_ID
count_Filtered$Probe_ID <- NULL

names(count_long)[2] <- "SAMPLE_ID"
count_long$Probe_ID <- NULL
count_long$COUNTS <- NULL
meta_Filtered <- count_long[!duplicated(count_long$SAMPLE_ID),]
rownames(meta_Filtered) <- meta_Filtered$SAMPLE_ID

idx <- sapply(rownames(meta_Filtered), function(x) {
  which(colnames(count_Filtered) == x)
})
count_Filtered <- count_Filtered[,idx]

#transforming the data

count_Filtered <- count_Filtered+1


###base Mean calculation##
CPM <- apply(count_Filtered, 2, function(x) (x/sum(x))*1000000)

dds <- DESeqDataSetFromMatrix(countData = count_Filtered,
                              colData = meta_Filtered,
                              design = ~ meanID)

if (NormMethod == "DESeq2"){ 
  dds <- estimateSizeFactors(dds) 
} else if (NormMethod == "CPM"){
  sizeFactors(dds) <- colSums(count_Filtered)/1000000
} else { 
  warning("Incorrect input for 'NormMethod', please double check above!")
}

baseMean <-  as.data.frame(rowMeans(counts(dds, normalized=TRUE)))
summary(baseMean)
baseMean$Probe_ID <- rownames(baseMean)
Low_Expr_Gene <- baseMean[baseMean$`rowMeans(counts(dds, normalized = TRUE))` < 0.558,]

#df_DEGs_ver3_clean <- df_DEGs_ver3[!df_DEGs_ver3$Probe_ID %in% Low_Expr_Gene_list,]
#hist(df_DEGs_ver3_clean$baseMean)
#summary(df_DEGs_ver3_clean$baseMean)
##merging Replicate by summing the count##

count_Filtered$probe <- rownames(count_Filtered)
count_Filtered <- count_Filtered[!count_Filtered$probe %in% Low_Expr_Gene$Probe_ID,]
count_Filtered$probe <- NULL

#count_Filtered_1 <- count_Filtered

#count_Filtered_1 <- lapply(count_Filtered_1, as.integer)
#count_Filtered_1 <- as.data.frame(count_Filtered_1)

#rownames(count_Filtered_1) <- rownames(count_Filtered)
##choosing the best replicate##


#meta_Filtered_best <- meta_Filtered
#meta_Filtered <- NA

#for(i in 1:length(unique(meta_Filtered_best$SAMPLENAME))){
  #SAM <- unique(meta_Filtered_best$SAMPLENAME)[i]
  #.meta_best <- meta_Filtered_best[meta_Filtered_best$SAMPLENAME == SAM,]
  #.meta_best <- .meta_best[order(.meta_best$LIB_SIZE, decreasing = TRUE),]
  #.meta_best <- .meta_best[!duplicated(.meta_best$SAMPLENAME),]
  #meta_Filtered <- rbind(.meta_best, meta_Filtered)
#}

#meta_Filtered <- meta_Filtered[c(1:79),]

## Normalization & DEG calculation ##
#count_Filtered <- count[, rownames(meta_Filtered)]

###Normalized it manually###


conc <- colnames(count_Filtered)
conc <- gsub("[0-9]", "", sapply(strsplit(conc, '_', perl = T), '[', 3))
conc <- unique(conc)
paths <- conc

calc_norm_gf <- function(paths){
  toMatch <- c(paths)
  count_sub_com <- count_Filtered[, grepl(paste0(toMatch), names(count_Filtered))]
  
  metafiltered_sub <- meta_Filtered[meta_Filtered$SAMPLE_ID %in% names(count_sub_com),]
  
  dds <- DESeqDataSetFromMatrix(countData = count_sub_com,
                                colData = metafiltered_sub,
                                design = ~ meanID)
  
  if (NormMethod == "DESeq2"){ 
    dds <- estimateSizeFactors(dds) 
  } else if (NormMethod == "CPM"){
    sizeFactors(dds) <- colSums(count_sub_com)/1000000
  } else { 
    warning("Incorrect input for 'NormMethod', please double check above!")
  }
  
  
  dds_compound <- DESeq(dds)
}

dds <- lapply(paths, calc_norm_gf)

##whole normalization##

#CPM <- apply(count_Filtered, 2, function(x) (x/sum(x))*1000000)

#dds <- DESeqDataSetFromMatrix(countData = count_Filtered,
                  #            colData = meta_Filtered,
                          #    design = ~ meanID)
#if (NormMethod == "DESeq2"){ 
 # dds <- estimateSizeFactors(dds) 
#} else if (NormMethod == "CPM"){
  #sizeFactors(dds) <- colSums(count_Filtered)/1000000
#} else { 
 # warning("Incorrect input for 'NormMethod', please double check above!")
#}

#dds <- DESeq(dds)

Norm_1 <- as.data.frame(counts(dds[[1]], normalized = TRUE)) #Normalized by size factors
Norm_2 <- as.data.frame(counts(dds[[2]], normalized = TRUE)) #Normalized by size factors
Norm_3 <- as.data.frame(counts(dds[[3]], normalized = TRUE)) #Normalized by size factors
Norm <- cbind(Norm_1, Norm_2, Norm_3)
#Norm <- as.data.frame(counts(dds, normalized = TRUE))
log2Norm <- log2(Norm + 1) #Log2 normalization of size factor corrected counts

###retrieving the log2FC###

meta_filetered_whole <- meta_Filtered[meta_Filtered$LOCATION_ID == "WHOLE",]
ControlSample_whole <- "BCL_SP_0195_LCM_WHOLE_CIS_conc0_TPcontrol"

meta_filetered_CPT <- meta_Filtered[meta_Filtered$LOCATION_ID == "CPT",]
ControlSample_CPT <- "BCL_SP_0195_LCM_CPT_CIS_conc0_TPcontrol"

meta_filetered_PPT <- meta_Filtered[meta_Filtered$LOCATION_ID == "PPT",]
ControlSample_PPT <- "BCL_SP_0195_LCM_PPT_CIS_conc0_TPcontrol"



samples_whole <- unique(meta_filetered_whole$meanID[which(!meta_filetered_whole$meanID %in% ControlSample_whole)])
samples_CPT <- unique(meta_filetered_CPT$meanID[which(!meta_filetered_CPT$meanID %in% ControlSample_CPT)])
samples_PPT <- unique(meta_filetered_PPT$meanID[which(!meta_filetered_PPT$meanID %in% ControlSample_PPT)])

data_list <- c("samples_whole","samples_CPT","samples_PPT")

paths <- data_list

calc_summaries <- function( paths ){
  
  if (paths == "samples_whole") {
    df_DEGs <- c() 
    for(y in 1:length(samples_whole)){
      sample <- samples_whole[y]
      print(paste(y, sample, ControlSample_whole, sep = ", "))
      rslt_tmp <- as.data.frame(results(dds[[3]], contrast = c("meanID", sample, ControlSample_whole)))
      rslt_tmp$meanID <- sample
      rslt_tmp$ref_meanID <- ControlSample_whole
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      df_DEGs <- rbind(df_DEGs, rslt_tmp)
    }
    return(df_DEGs)
  }else if (paths == "samples_CPT"){
    df_DEGs <- c()
    for(y in 1:length(samples_CPT)){
      sample <- samples_CPT[y]
      print(paste(y, sample, ControlSample_CPT, sep = ", "))
      rslt_tmp <- as.data.frame(results(dds[[1]], contrast = c("meanID", sample, ControlSample_CPT)))
      rslt_tmp$meanID <- sample
      rslt_tmp$ref_meanID <- ControlSample_CPT
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      df_DEGs <- rbind(df_DEGs, rslt_tmp)
    }
    return(df_DEGs)
  } else if (paths == "samples_PPT"){
    df_DEGs <- c()  
    for(y in 1:length(samples_PPT)){
      sample <- samples_PPT[y]
      print(paste(y, sample, ControlSample_PPT, sep = ", "))
      rslt_tmp <- as.data.frame(results(dds[[2]], contrast = c("meanID", sample, ControlSample_PPT)))
      rslt_tmp$meanID <- sample
      rslt_tmp$ref_meanID <- ControlSample_PPT
      rslt_tmp$Probe_ID <- rownames(rslt_tmp)
      rownames(rslt_tmp) <- c()
      df_DEGs <- rbind(df_DEGs, rslt_tmp)
    }
    return(df_DEGs) 
  } else { stop ("variable not found")}
}

summary_lists <- lapply(paths, calc_summaries)
df_DEGs <- do.call(rbind,summary_lists)

###gene matching region rat##
#upload mouse gene Park, Jihwan, et al.#
Mouse_geneList <- read.csv("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2020/20200501_LCM_RealExperiment/Data/Mouse_GeneList.txt",
                             sep ="")

Gene_ortholog <- read.csv("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2020/20200501_LCM_RealExperiment/Data/RGD_ORTHOLOGS_NCBI.txt",
                           sep ="\t")

Gene_ortholog <- Gene_ortholog[,c(1:3,8:10)]
Gene_ortholog <- Gene_ortholog[!duplicated(Gene_ortholog$RAT_GENE_SYMBOL),]
Gene_ortholog <- Gene_ortholog[Gene_ortholog$MOUSE_ORTHOLOG_SYMBOL %in% Mouse_geneList$genes,]

names(Mouse_geneList)[1] <- "MOUSE_ORTHOLOG_SYMBOL"
Mouse_geneList_1 <- Mouse_geneList[Mouse_geneList$MOUSE_ORTHOLOG_SYMBOL %in% Gene_ortholog$MOUSE_ORTHOLOG_SYMBOL,]
Mouse_geneList_1 <- Mouse_geneList_1[!duplicated(Mouse_geneList_1$MOUSE_ORTHOLOG_SYMBOL),]

Gene_kidney <- merge(Gene_ortholog, Mouse_geneList_1, by = "MOUSE_ORTHOLOG_SYMBOL")
Gene_annotation <- Gene_kidney[,c(2,7)]
log2Norm$Gene_symbol <- rownames(log2Norm)
log2Norm$Gene_symbol<-gsub("(.*)_.*","\\1",log2Norm$Gene_symbol)
log2Norm <- log2Norm[log2Norm$Gene_symbol %in% Gene_annotation$RAT_GENE_SYMBOL,]
Gene_annotation <- Gene_annotation[Gene_annotation$RAT_GENE_SYMBOL %in% log2Norm$Gene_symbol,]
log2Norm <- log2Norm[!duplicated(log2Norm$Gene_symbol),]
log2Norm$Gene_symbol <- NULL


##heatmap_time##

##loop##
for (i in 1:length(unique(meta_Filtered$CONCENTRATION))){
  DOSE <- unique(meta_Filtered$CONCENTRATION)[i]
  meta_sub <- meta_Filtered[meta_Filtered$CONCENTRATION == DOSE,]
  log2Norm_sub <- log2Norm[,colnames(log2Norm) %in% meta_sub$SAMPLE_ID]
  
heatmap_final <- transpose(log2Norm_sub)
rownames(heatmap_final) <- colnames(log2Norm_sub)
colnames(heatmap_final) <- rownames(log2Norm_sub)

heatmap_final[1:10,1:10]


Gene_annotation <- Gene_kidney[,c(2,7)]
log2Norm_sub$Gene_symbol <- rownames(log2Norm_sub)
log2Norm_sub$Gene_symbol<-gsub("(.*)_.*","\\1",log2Norm_sub$Gene_symbol)
log2Norm_sub <- log2Norm_sub[log2Norm_sub$Gene_symbol %in% Gene_annotation$RAT_GENE_SYMBOL,]
Gene_annotation <- Gene_annotation[Gene_annotation$RAT_GENE_SYMBOL %in% log2Norm_sub$Gene_symbol,]
log2Norm_sub <- log2Norm_sub[!duplicated(log2Norm_sub$Gene_symbol),]
log2Norm_sub$Gene_symbol <- NULL

myColors = list()
require(RColorBrewer)
colFun <- colorRampPalette(c("red","green","blue","yellow", "purple", "brown"))
myColors$Annotation <- colFun(6)
names(myColors$Annotation) <- unique(Gene_annotation$Annotation)

colnames(heatmap_final)


# only cluster row based treatment (leave dose-range)
myColors$probe <- c("white","darkmagenta")
myColors$probe <- colFun(410)
names(myColors$probe) <- unique(Gene_annotation$probe)

head(heatmap_final)


# color anotations for columns. 
colnames(heatmap_final)


heatmap_final_color <- heatmap_final

paletteLength <- 300
my.colors <- colorRampPalette(c("white", "red"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
colBreaks <- c(seq(min(heatmap_final_color), 0, length.out=ceiling(paletteLength/2) + 1), 
               seq(max(heatmap_final_color)/paletteLength, max(heatmap_final_color), length.out=floor(paletteLength/2)))

#my.colors  <- colorRampPalette(c("blue", ="white", "red"), space="Lab")(n = 300)


pie(rep(1, length(my.colors)), labels = sprintf("%d (%s)", seq_along(my.colors), 
                                                my.colors), col = my.colors)


log2Norm_merge <- log2Norm_sub
log2Norm_merge$probe <- rownames(log2Norm_merge)
log2Norm_merge$RAT_GENE_SYMBOL <- rownames(log2Norm_merge)
log2Norm_merge$RAT_GENE_SYMBOL<-gsub("(.*)_.*","\\1",log2Norm_merge$RAT_GENE_SYMBOL)
log2Norm_merge <- log2Norm_merge[,c("probe", "RAT_GENE_SYMBOL")]
Gene_annotation <- merge(Gene_annotation, log2Norm_merge)


##clustering
#Sample

aggrtreats <- heatmap_final

dist_hm_data <- dist(aggrtreats, method = c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1])
clust_hm_data<-hclust(dist_hm_data, method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")[2])

pdf(paste0("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2021/LCM_Part2/graph/Sample_Clustering_Basedon_KidneyMarker_latest_new_20211109.pdf"),
    height =5, width = 25)
plot(clust_hm_data, hang = -1, xlab = 12)
dev.off()

tmp_1 <- table(rownames(heatmap_final))[clust_hm_data$order]

indl = alist()
for(i in seq_along(tmp_1)) {
  indl[[i]] <- which(rownames(heatmap_final) %in% names(tmp_1)[i])
  
}
ind_rows<- unlist(indl)
rownames(heatmap_final)[ind_rows]

#Gene

#aggrtreats <- transpose(heatmap_final)
#rownames(aggrtreats) <- colnames(heatmap_final)
#colnames(aggrtreats) <- rownames(heatmap_final)

#dist_hm_data <- dist(aggrtreats, method = c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1])
#clust_hm_data<-hclust(dist_hm_data, method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")[2])

#pdf(paste0("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2020/20200501_LCM_RealExperiment/Probe_Clustering_Basedon_KidneyMarker_clean10th_mergedSum_",DOSE,".pdf"),
 #   height =5, width = 70)
#plot(clust_hm_data, hang = -1, xlab = 12)
#dev.off()

##order alphabetical annot
#probe_list <- clust_hm_data$labels[clust_hm_data$order]  

Gene_annotation <- Gene_annotation[order(Gene_annotation$Annotation),]
probe_list <- Gene_annotation$probe
idx <- sapply(probe_list, function(x) {
  which(Gene_annotation$probe == x)
})
Gene_annotation <- Gene_annotation[idx,]
rownames(Gene_annotation) <- Gene_annotation$probe
ind_cols <- match(rownames(Gene_annotation), colnames(heatmap_final))
rownames(Gene_annotation) <- colnames(heatmap_final)[ind_cols]

if(nrow(Gene_annotation) != ncol(heatmap_final)){
  stop("check annCol: did you change feature number or sampled time points?")
}

Gene_annotation$probe <- NULL

dev.off()

pdf(paste0("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2021/LCM_Part2/graph/HeatmapSubset_KidneyMarker_probenotclus_new_20211109.pdf"), height =15, width = 10)
pheatmap(heatmap_final[ind_cols], cluster_cols = FALSE, cluster_rows = TRUE, cutree_rows = 4,
         clustering_method =  "ward.D2", border_color = "grey90",
         fontsize = 12, width = 10, height=10, legend = TRUE,
         annotation_col = Gene_annotation, annotation_colors = myColors,
         fontsize_col = 2, show_rownames = T, show_colnames = T 
)
dev.off()
}
### percetnage mapped plot###
for(i in 1:length(unique(meta$TIMEPOINT))){
  DOSE <- unique(meta$TIMEPOINT)[i] 
  
  meta_sub <- meta[meta$TIMEPOINT == DOSE,]
  meta_sub$Th3Percent <- meta_sub$PERCENTAGE_MAPPED > 10

pdf(paste0(outputDir_GR, "percentage_mapped_.pdf"), width = 10, height = 4)
print(ggplot(meta_sub, aes(x = SAMPLE_ID, y = PERCENTAGE_MAPPED)) +
        geom_hline(yintercept = 10, colour = 'grey75', size = 1) +
        geom_point(aes(reorder(SAMPLE_ID, PERCENTAGE_MAPPED),color = Th3Percent)) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12))+
        ggtitle("") + ylab('PERCENTAGE_MAPPED') + xlab(VariableQCs))
dev.off()


}

##hierarchical clustering###

log2NormFull <- log2Norm
meta_Filtered_full <- meta_Filtered
for(i in 1:length(unique(meta_Filtered_full$LOCATION_ID))){
  DOSE <- unique(meta_Filtered_full$LOCATION_ID)[i]
  meta_Filtered <- meta_Filtered_full[meta_Filtered_full$LOCATION_ID == DOSE,]
  sample <- meta_Filtered$SAMPLE_ID
  log2Norm <- log2NormFull[, sample]
  log2Norm_t <- transpose(log2Norm)
  colnames(log2Norm_t) <- rownames(log2Norm)
  rownames(log2Norm_t) <- colnames(log2Norm)
  aggrtreats <- log2Norm_t
  
  dist_hm_data <- dist(aggrtreats, method = c("euclidean", "maximum", "manhattan", "canberra", "binary" , "minkowski")[1])
  clust_hm_data<-hclust(dist_hm_data, method = c("ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median","centroid")[2])
  
  pdf(paste0("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2021/LCM_Part2/Replicate_2/SampleCLustering_allSample_latest_20211109_",DOSE,".pdf"),
      height =5, width = 20)
  plot(clust_hm_data, hang = -1, xlab = 12)
  dev.off()
  
}

##eliminating genes with total count < 5
Norm_long <- melt(as.matrix(Norm))
colnames(Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
Norm_long <- left_join(Norm_long, meta_Filtered[,c("SAMPLE_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID","LOCATION_ID")])
Norm_long$REPLICATE <- NULL

Norm_long_sum <- aggregate( COUNTS ~ Probe_ID + EXP_ID + TIMEPOINT + TREATMENT + CONCENTRATION  + meanID + LOCATION_ID, data = Norm_long, sum)

Norm_long_sum$gene <- Norm_long_sum$Probe_ID
Norm_long_sum$gene <-gsub("_.*","",Norm_long_sum$gene)

#elimnintating gene with norm count < 5
#Norm_long_sum <- Norm_long_sum[Norm_long_sum$COUNTS > 5,]

##TXG-MAPr_gene_Upload

path <- "/media/callegarog/GiuliaData_2017/"
jDrive <- "/vol/winshare/Public/"
if( Sys.info()[1] != "Linux" ) {
  path <- "H:/"
  jDrive <- "J:/"}

load(file.path(jDrive, 'Workgroups/FWN/LACDR/TOX/VDWATER/WGCNA/KIDNEY_RData/Kidney_ParsedModules_data.RData'))
modules   <<- unique(datKidney$module_definitions)
percentage_gene <- NA
for (i in 1:length(unique(Norm_long_sum$LOCATION_ID))){
  LOC <- unique(Norm_long_sum$LOCATION_ID)[i]
  norm_subset <- Norm_long_sum[Norm_long_sum$LOCATION_ID == LOC,]
  
  for(a in 1:length(unique(norm_subset$meanID))){
    MEAN <- unique(norm_subset$meanID)[a]
    norm_subset_2 <- norm_subset[norm_subset$meanID == MEAN,]
    
    norm_subset_2 <- norm_subset_2[norm_subset_2$gene %in% modules$rat.gene.symbol,]
    percentage_included_gene <- (nrow(norm_subset_2)/nrow(modules))*100
    percentage_gene_sub <- as.data.frame(cbind(unique(norm_subset_2$meanID), percentage_included_gene))
    percentage_gene <- rbind(percentage_gene_sub, percentage_gene)
    
}
}



#log2Norm_1 <- data.frame(log2(CPM + 1)) #Log2 normalization of size factor corrected counts

## Retrieving results table of desired comparisons ##

#df_DEGs <- c()
#samples <- unique(meta_Filtered$meanID[which(!meta_Filtered$meanID %in% ControlSample)])


#for(y in 1:length(samples)){
 # sample <- samples[y]
  
  #print(paste(y, sample, ControlSample, sep = ", "))
  
  #rslt_tmp <- as.data.frame(results(dds, contrast = c("meanID", sample, ControlSample)))
  
  #rslt_tmp$meanID <- sample
  #rslt_tmp$ref_meanID <- ControlSample
  #rslt_tmp$Probe_ID <- rownames(rslt_tmp)
  #rownames(rslt_tmp) <- c()
  #df_DEGs <- rbind(df_DEGs, rslt_tmp)
#}

#df_DEGs_1 <- read.csv("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2020/20200501_LCM_RealExperiment/Data/DEGsfiltered_vs_correctControl_clean_transformed_corrClean_ver5_normSplit.txt", sep="\t")

df_DEGs <- unique(left_join(df_DEGs, meta_Filtered[,c("EXP_ID","TIMEPOINT","TREATMENT","CONCENTRATION","meanID","LOCATION_ID" )], by = "meanID"))
df_DEGs <- separate(df_DEGs, c("Probe_ID"), into = c("GeneSymbol", "ProbeNr"), sep = "_", remove = FALSE)

path <- "/media/callegarog/GiuliaData_2017/"
jDrive <- "/vol/winshare/Public/"
if( Sys.info()[1] != "Linux" ) {
  path <- "H:/"
  jDrive <- "J:/"}

load(file.path(jDrive, 'Workgroups/FWN/LACDR/TOX/VDWATER/WGCNA/KIDNEY_RData/Kidney_ParsedModules_data.RData'))
modules   <- symbol_rat[!duplicated(symbol_rat$entrezgene),]
modules   <- symbol_rat[,c("entrezgene","external_gene_name"),]

names(modules)[2] <- "GeneSymbol"

df_DEGs_filtered <- df_DEGs[df_DEGs$GeneSymbol %in% modules$GeneSymbol,]
modules <- modules[modules$GeneSymbol %in% df_DEGs_filtered$GeneSymbol,]

df_DEGs_filtered_1 <- merge(df_DEGs_filtered, modules, by  = c("GeneSymbol"))
df_DEGs_filtered <- df_DEGs_filtered_1
#df_DEGs$Entrez_ID <- unname(mapIds(org.Hs.eg.db, keys= df_DEGs$GeneSymbol, column="ENTREZID", keytype="ALIAS", multiVals="first"))

## Quality control of data ##

# Dataframe for plots
meta_Filtered <- meta_Filtered_full
count_long <- melt(as.matrix(count_Filtered))
colnames(count_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
count_long <- left_join(count_long, meta_Filtered[,c("SAMPLE_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

Norm_long <- melt(as.matrix(Norm))
colnames(Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
Norm_long <- left_join(Norm_long, meta_Filtered[,c("SAMPLE_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

meta_Filtered$TIMEPOINT <- as.factor(meta_Filtered$TIMEPOINT)
meta_Filtered$REPLICATE <- as.factor(meta_Filtered$REPLICATE)


log2Norm_long <- melt(as.matrix(log2Norm))
colnames(log2Norm_long) <- c("Probe_ID", "SAMPLE_ID", "COUNTS")
log2Norm_long <- left_join(log2Norm_long, meta_Filtered[,c("SAMPLE_ID", "EXP_ID", "TIMEPOINT", "TREATMENT", "CONCENTRATION", "REPLICATE", "meanID", "LOCATION_ID")])

#log2Norm_mdm2<- log2Norm_long[log2Norm_long$Probe_ID == "Mdm2_9214",]


# Library size distribution #
for(i in 1:length(unique(meta$CONCENTRATION))){
 DOSE <- unique(meta$CONCENTRATION)[i] 
 for (u in 1 : length(unique(meta$LOCATION_ID))){
   LOC <- unique(meta$LOCATION_ID)[u] 
 
meta_sub <- meta[meta$CONCENTRATION == DOSE,]
meta_sub_1 <- meta_sub[meta_sub$LOCATION_ID == DOSE,]

meta_sub_1$VariableQCs = meta_sub_1[,which(names(meta_sub_1)==VariableQCs)]
pdf(paste0(outputDir_GR, "LibSize__complete_techRep_20211109.pdf"), width = 10, height = 5)
print(ggplot(meta_sub_1, aes(x = reorder(VariableQCs, LIB_SIZE, FUN = median), y = LIB_SIZE)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) + geom_hline(yintercept = 500000, colour = 'red', size = 0.5)+
        scale_y_log10(limits = c(1, max(meta_sub_1$LIB_SIZE))) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Library size distribution") + ylab('Library size') + xlab(VariableQCs))
dev.off()
}
}
# Counts & normalized counts distribution amoung CELL IDs #
for(i in 1:length(unique(count_long$CONCENTRATION))){
  DOSE <- unique(count_long$CONCENTRATION)[i] 
  for (u in 1 : length(unique(meta$LOCATION_ID))){
    LOC <- unique(meta$LOCATION_ID)[u] 
    
  count_long_sub <- count_long[count_long$CONCENTRATION == DOSE,]
  count_long_sub <- count_long[count_long$LOCATION_ID == LOC,]
  
count_long_sub$VariableQCs = count_long_sub[,which(names(count_long_sub)==VariableQCs)]
pdf(paste0(outputDir_GR, "Distribution-Counts_combinedControled_20211109.pdf"), width = 10, height = 5)
print(ggplot(count_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS+1)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        scale_y_log10(limits = c(1, max(count_long_sub$COUNTS))) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +
        ggtitle("Distribution counts") + ylab('counts') + xlab(VariableQCs))
dev.off()
}
}

for(i in 1:length(unique(log2Norm_long$CONCENTRATION))){
  DOSE <- unique(log2Norm_long$CONCENTRATION)[i] 
  for (u in 1 : length(unique(meta$LOCATION_ID))){
    LOC <- unique(meta$LOCATION_ID)[u] 
    
  log2Norm_long_sub <- log2Norm_long[log2Norm_long$CONCENTRATION == DOSE,]
  log2Norm_long_sub <- log2Norm_long[log2Norm_long$LOCATION_ID == LOC,]
log2Norm_long_sub$VariableQCs = log2Norm_long_sub[,which(names(log2Norm_long_sub)==VariableQCs)] 
pdf(paste0(outputDir_GR, "Distribution-CountsNorm_combinedControled_20211109.pdf"), width = 10, height = 5)
print(ggplot(log2Norm_long_sub, aes(x = reorder(VariableQCs, COUNTS, FUN = median), y = COUNTS)) +
        geom_boxplot(size = 0.3, outlier.size = 0.5) +
        theme_classic() +
        theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
              axis.title.x = element_text(size=12, vjust = 0.25),
              axis.title.y = element_text(size=12, vjust = 1),
              axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
              axis.text.y = element_text(size=12)) +        ggtitle("Distribution Normalized counts") + ylab('Normalized counts') + xlab(VariableQCs))
dev.off()
}
}
# SDs across means after different normalization methods

notAllZero <- (rowSums(counts(dds[1]))>0)

pdf(paste0(outputDir_GR, "SDvsMean_log2Norm_control360h.pdf"), width = 10, height = 8)
meanSdPlot(log2(counts(dds[1],normalized=TRUE)[notAllZero,] + 1))
dev.off()

#meta_Filtered$readClass <- NA
#meta_Filtered[meta_Filtered$LIB_SIZE > 100000,]$readClass <- "ReadAbove_100000"
#meta_Filtered[meta_Filtered$LIB_SIZE <= 100000,]$readClass <- "ReadBelow_100000"


# Replicate correlation #

meta_sub_1 <- meta_Filtered_full

repCors <- lapply(unique(meta_sub_1[, "meanID"]), function(expID) {
#expID <- "EUT_076_LCM_Whole_CIS_conc0_TP24"
  replicateIDs <- rownames(meta_sub_1)[which(meta_sub_1[, "meanID"] == expID)]
  if(length(replicateIDs)>1){
    experimentMean <- rowMeans(log2Norm[, replicateIDs])
    apply(log2Norm[, replicateIDs], 2, function(.expID) cor(.expID, experimentMean))
  }
})

repCors <- lapply(repCors,function(y){cbind("pearsonR" = y, "SAMPLE_ID" = names(y))})
repCors <- as.data.frame(do.call(rbind, repCors))
repCors$pearsonR <- as.numeric(as.character(repCors$pearsonR))
repCors$clr <- repCors$pearsonR > 0.8
repCors <- repCors[which(!is.na(repCors$pearsonR)),]
repCors <- left_join(repCors, meta_sub_1)
#repCors$CELL_ID <- as.factor(repCors$CELL_ID)
repCors$VariableQCs = repCors[,which(names(repCors)==VariableQCs)]
pdf(paste0(outputDir_GR, "RepCor_PearsonR_CorrClean_CombinedControled_log2Norm_20211109.pdf"), width = 10, height = 10)
print(ggplot(repCors, aes(x = reorder(VariableQCs, pearsonR, FUN = min), y = pearsonR)) +
  geom_hline(yintercept = 0.8, colour = 'grey75', size = 0.75) +
  geom_point(aes()) +
  theme_classic() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=12))+
  ggtitle("Replicate correlation") + ylab('PearsonR') + xlab(VariableQCs))
dev.off()


repCors_lowPR <- filter(repCors, clr == FALSE)

## Plot of dispersion estimate and MA-plot ##

pdf(paste0(outputDir_GR, "MA-plot_control360h.pdf"), width = 6, height = 4)
plotMA(results(dds), ylim=c(-2,2))
dev.off()

pdf(paste0(outputDir_GR, "Dispersion-Estimate_control360h.pdf"), width = 6, height = 4)
plotDispEsts(dds)
dev.off()

## Number of DEGs ##

df_nrDEGs <- ddply(df_DEGs, .(meanID, ref_meanID, TREATMENT, TIMEPOINT, EXP_ID, LOCATION_ID), summarize, 
                   NrDEGs = sum(padj < Threshold_padj , na.rm = TRUE),
                   NrDEGsFC = sum(padj < Threshold_padj & abs(log2FoldChange) > Threshold_log2FC , na.rm = TRUE))

if(Filtering == "padj"){
  df_DEGs_filtered <- subset(df_DEGs, padj < Threshold_padj)
} else {
  if(Filtering == "padj&log2FC"){
    df_DEGs_filtered <- subset(df_DEGs, padj < Threshold_padj & abs(log2FoldChange) > Threshold_log2FC)
  } else {
    if(Filtering == "log2FC"){
      df_DEGs_filtered <- subset(df_DEGs, abs(log2FoldChange) > Threshold_log2FC)
    } else {
      warning("Unrecognized input for Filtering variable. Check what has been filled in for Filtering at 'Input user' section")
    }
  }
}

df_DEGs_FCmatrix <- dcast(Count_batch2, Probe_ID ~ meanID, value.var = "log2FoldChange")
rownames(df_DEGs_FCmatrix) <- df_DEGs_FCmatrix[,1]
df_DEGs_FCmatrix <- df_DEGs_FCmatrix[,-1]

meta_Filtered_meanID <- unique(meta_Filtered[, c("CELL_ID", "TREATMENT", "CONCENTRATION", "TIMEPOINT", "EXP_ID", "meanID")])
rownames(meta_Filtered_meanID) <- meta_Filtered_meanID$meanID

## Separate text data frames for each sample for Ingenuity Pathway Analysis ##

for(i in 1:length(unique(df_DEGs$meanID))){
  df_DEGs_sub <- df_DEGs[which(df_DEGs$meanID == unique(df_DEGs$meanID)[i]), c("GeneSymbol", "log2FoldChange", "padj")]
  write.table(df_DEGs_sub, file = paste0(output_IPA, unique(df_DEGs$meanID)[i], ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
}

## Distribution baseMean ##

for(i in 1:length(unique(df_DEGs_filtered$CONCENTRATION))){
  DOSE <- unique(df_DEGs_filtered$CONCENTRATION)[i] 
  count_long_sub <- df_DEGs_filtered[df_DEGs_filtered$CONCENTRATION == DOSE,]

  count_long_sub$VariableQCs = count_long_sub[,which(names(count_long_sub)==VariableQCs)]
pdf(paste0(outputDir_GR, "Histogram_baseMean_combinedControl.pdf"), width = 30, height = 30)
print(ggplot(count_long_sub, aes(x = baseMean)) +
  geom_histogram(binwidth = 0.1) +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  ggtitle("Histogram baseMean") + ylab('Count') + xlab('Log10 baseMean') +
  facet_wrap(.~ VariableQCs))
dev.off()

}

## Distribution of log2FC and p-values ##
for(i in 1:length(unique(df_DEGs_filtered$CONCENTRATION))){
  DOSE <- unique(df_DEGs_filtered$CONCENTRATION)[i] 
  count_long_sub <- df_DEGs_filtered[df_DEGs_filtered$CONCENTRATION == DOSE,]

pdf(paste0(outputDir_GR, "Log2FC_vs_baseMean_CorrClean_combinedControl_20211109.pdf"), width = 30, height = 30)
print(ggplot(count_long_sub, aes(x = baseMean, y  = log2FoldChange)) +
  geom_point(aes(color = padj < Threshold_padj), size = 0.2, alpha = 0.3) +
  geom_hline(aes(yintercept = Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  geom_hline(aes(yintercept = -Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "red", `NA` = "darkgrey")) +
  scale_x_continuous(expand = c(0, 0), trans = "log10") +
  scale_y_continuous(expand = c(0, 0), limits = c(-20, 20)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  ggtitle("Log2FC vs baseMean") + ylab('Log2FC (vs control)') + xlab('Log10 baseMean') +
  facet_wrap(LOCATION_ID~ meanID))
dev.off()

}

df_DEGs_filtered_com <- df_DEGs_filtered[!df_DEGs_filtered$CONCENTRATION == 0,]
for(i in 1:length(unique(df_DEGs_filtered_com$TIMEPOINT))){
  DOSE <- unique(df_DEGs_filtered_com$TIMEPOINT)[i] 
  count_long_sub <- df_DEGs_filtered_com[df_DEGs_filtered_com$TIMEPOINT == DOSE,]
  
pdf(paste0(outputDir_GR, "padj_vs_log2FC_CombinedControl_20211109.pdf"), width = 30, height = 30)
print(ggplot(count_long_sub, aes(x = log2FoldChange, y  = -log10(padj))) +
  geom_point(aes(color = padj < Threshold_padj), size = 0.2, alpha = 0.3) +
  geom_vline(aes(xintercept = Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  geom_vline(aes(xintercept = -Threshold_log2FC), color = "dodgerblue", size = 0.6, alpha = 0.5) +
  scale_color_manual(values = c(`FALSE` = "darkgrey", `TRUE` = "red", `NA` = "darkgrey")) +
  scale_x_continuous(expand = c(0, 0), limits = c(-10, 10)) +
  scale_y_continuous(expand = c(0, 0), limits =c(0, 15)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  ggtitle("Adjusted p-value vs Log2FC") + ylab('-log10 adjusted p-value') + xlab('log2FC') +
  facet_wrap(LOCATION_ID ~ meanID))
dev.off()
}
## Check gene expression of individual genes 

# Log2FC
df_DEGs_GenesManual <- PathwayProbesCounts(Genelist = CheckGenes, CountTable = df_DEGs, ProbeCol = "Probe_ID", GeneID = "GeneSymbol", match = "Entrez")
pdf(paste0(outputDir_GR, "GenesManual_DESeq2Log2FC_", ControlSample, ".pdf"), width = 15, height = 15)
ggplot(df_DEGs_GenesManual, aes(x = meanID, y = log2FoldChange)) + 
  geom_bar(aes(fill = padj < Threshold_padj), stat="identity", position=position_dodge(), color = "black") +
  geom_errorbar(aes(ymin = log2FoldChange - lfcSE, ymax = log2FoldChange + lfcSE), 
                color = "black",  position=position_dodge(0.9), width = 0.5) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  facet_grid(Probe_ID ~ CELL_ID, scales = "free") +
  ggtitle("Log2FC") + ylab('Log2FC (vs control)') + xlab('Conditions')
dev.off()

#Log2 normalized counts
log2Norm_long_GenesManual <- PathwayProbesCounts(Genelist = CheckGenes, CountTable = log2Norm_long, ProbeCol = "Probe_ID", GeneID = "GeneSymbol", match = "Entrez")
log2Norm_long_GenesManual_sum <- setDT(log2Norm_long_GenesManual)[,.(meanLog2Count = mean(COUNTS), sdLog2Count = sd(COUNTS)), by = .(TREATMENT, CELL_ID, TIMEPOINT, Probe_ID, EXP_ID, meanID, CONCENTRATION)]
pdf(paste0(outputDir_GR, "GenesManual_Log2count_mean.pdf"), width = 15, height = 15)
ggplot(log2Norm_long_GenesManual_sum, aes(x = meanID, y = meanLog2Count)) + 
  geom_bar(stat="identity", position=position_dodge(), color = "black", fill = "steelblue") +
  geom_errorbar(aes(ymin = meanLog2Count, ymax = meanLog2Count + sdLog2Count), 
                color = "black",  position=position_dodge(0.9), width = 0.5) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  facet_grid(Probe_ID ~ CELL_ID + TREATMENT + TIMEPOINT, scales = "free") +
  ggtitle("Log2 normalized counts") + ylab('Log2 normalized counts') + xlab('Conditions')
dev.off()

pdf(paste0(outputDir_GR, "GenesManual_Log2count_reps.pdf"), width = 15, height = 15)
ggplot(log2Norm_long_GenesManual, aes(x = meanID, y = COUNTS)) + 
  geom_point(aes(color = factor(REPLICATE)), alpha = 0.8) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 20)) +
  theme_bw() +
  theme(plot.title = element_text(size=14, face="bold", vjust = 2, hjust = 0.5), 
        axis.title.x = element_text(size=12, vjust = 0.25),
        axis.title.y = element_text(size=12, vjust = 1),
        axis.text.x = element_text(size=8, angle=90, vjust=0.5, hjust=1),
        axis.text.y = element_text(size=8)) +
  facet_grid(Probe_ID ~ CELL_ID + TREATMENT + TIMEPOINT, scales = "free") +
  ggtitle("Log2 normalized counts") + ylab('Log2 normalized counts') + xlab('Conditions')
dev.off()

## Subsetting of data for specific geneset
if(!is.na(Geneset)){
  log2Norm_Genes <- PathwayProbesCounts(Genelist = Geneset, CountTable = log2Norm, ProbeCol = "rownames", GeneID = "Probes", match = "Entrez", Dir = GeneListDir)
  df_DEGs_FCmatrix_Genes <- PathwayProbesCounts(Genelist = Geneset, CountTable = df_DEGs_FCmatrix, ProbeCol = "rownames", GeneID = "Probes", match = "Entrez", Dir = GeneListDir)
}

## Principle component analysis ##
meta$percetangeMapped_cat <- NA
meta$LIB_SIZE <- NA
meta[meta$PERCENTAGE_MAPPED < 10,]$percetangeMapped_cat <- "mapped_below_10%"
meta[meta$PERCENTAGE_MAPPED >= 10,]$percetangeMapped_cat <- "mapped_Above_10%"
#meta_Filtered[meta_Filtered$Percentage_mapped < 10 & meta_Filtered$Percentage_mapped > 5,]$percetangeMapped_cat <- "mapped_between 5_10%"

meta$CELL_ID <- as.factor(meta$CELL_ID)
meta_Filtered
#meta_Filtered <- meta_Filtered_full
#for(i in 1:length(unique(meta_Filtered_full$CONCENTRATION))){
 # DOSE <- unique(meta_Filtered_full$CONCENTRATION)[i]
#meta_Filtered <- meta_Filtered_full[meta_Filtered_full$CONCENTRATION == DOSE,]
#sample <- meta_Filtered$SAMPLE_ID
#log2Norm <- log2NormFull[, sample]
  dev.off()
df_log2Norm_pca <- PCADF(log2Norm)
meta_Filtered$TIMEPOINT <- ordered(meta_Filtered$TIMEPOINT, levels = c("control", "1", "2", "4","24","72","120","192","240","288","360","480","672"))
PCAplots(df_log2Norm_pca, meta_Filtered, colourVar = "TIMEPOINT", Clusts = FALSE, FileID = "Raw_count_Filtered-AllSamples-AllGenes_th_new", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)

#}

PCAplots(df_log2Norm_pca, meta_Filtered, colourVar = "T", Clusts = FALSE, FileID = "Log2-AllSamples-AllGenes_20211109_", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)
pca_ALL_top25_log2Norm <- PCAcontrVar(df_log2Norm_pca, NrGenes = 25)
pcaCoV(df_log2Norm_pca, meta_Filtered, filename = "Log2-AllSamples-AllGenes")

if(!is.na(Geneset)){
df_log2Norm_Genes_pca <- PCADF(log2Norm_Genes)
PCAplots(df_log2Norm_Genes_pca, meta_Filtered, colourVar = "TREATMENT", Clusts = FALSE, FileID = "Log2-AllSamples-Geneset", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)
PCAplots(df_log2Norm_Genes_pca, meta_Filtered, colourVar = "CELL_ID", Clusts = FALSE, FileID = "Log2-AllSamples-Geneset", width = 10, height = 8, sizeP = 0.2, Dataclass = "log2", a = 0.7)
pca_ALL_top25_log2Norm_Genes <- PCAcontrVar(df_log2Norm_Genes_pca, NrGenes = 25)
pcaCoV(df_log2Norm_Genes_pca, meta_Filtered, filename = "Log2-AllSamples-Geneset")
}

df_log2FC_pca <- PCADF(na.omit(df_DEGs_FCmatrix))
PCAplots(df_log2FC_pca, meta_Filtered_meanID, colourVar = "TREATMENT", Clusts = FALSE, FileID = paste0("Log2FC-AllSamples-AllGenes_", ControlSample), width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
PCAplots(df_log2FC_pca, meta_Filtered_meanID, colourVar = "CELL_ID", Clusts = FALSE, FileID = paste0("Log2FC-AllSamples-AllGenes_", ControlSample), width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
pca_ALL_top25genes_log2FC <- PCAcontrVar(df_log2FC_pca, NrGenes = 25)
pcaCoV(df_log2FC_pca, meta_Filtered_meanID, filename = "Log2FC-AllSamples-AllGenes")

if(!is.na(Geneset)){
df_log2FC_Genes_pca <- PCADF(na.omit(df_DEGs_FCmatrix_Genes))
PCAplots(df_log2FC_Genes_pca, meta_Filtered_meanID, colourVar = "TREATMENT", Clusts = FALSE, FileID = paste0("Log2FC-AllSamples-Geneset_", ControlSample), width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
PCAplots(df_log2FC_Genes_pca, meta_Filtered_meanID, colourVar = "CELL_ID", Clusts = FALSE, FileID = paste0("Log2FC-AllSamples-Geneset_", ControlSample), width = 10, height = 8, sizeP = 0.2, Dataclass = "log2FC", a = 0.7)
pca_ALL_top25genes_log2FC_Genes <- PCAcontrVar(df_log2FC_Genes_pca, NrGenes = 25)
pcaCoV(df_log2FC_Genes_pca, meta_Filtered_meanID, filename = "Log2FC-AllSamples-Geneset")
}

## Heatmap ##

Heatmap(Data = log2Norm, Meta = meta_Filtered, FC = FALSE, heightGene = 0.2, widthGene = 7, FileName = "Log2-AllSamples-AllGenes")
Heatmap(Data = df_DEGs_FCmatrix, Meta = meta_Filtered_meanID, FC = TRUE, heightGene = 0.2, widthGene = 7, FileName = "Log2FC-AllSamples-AllGenes")

if(!is.na(Geneset)){
  Heatmap(Data = log2Norm_Genes, Meta = meta_Filtered, FC = FALSE, heightGene = 7, widthGene = 7, FileName = paste0("Log2FC-AllSamples-Geneset_", ControlSample), GeneName = TRUE)
  Heatmap(Data = df_DEGs_FCmatrix_Genes, Meta = meta_Filtered_meanID, FC = TRUE, heightGene = 7, widthGene = 7, FileName = paste0("Log2FC-AllSamples-Geneset_", ControlSample), GeneName = TRUE)
}

## Saving data ##	

setwd(inputDir)

write.table(meta_Filtered, file = 'meta_Filtered_20211109.txt', sep = '\t')
write.table(meta_Filtered_meanID, file = 'meta_Filtered_meanID.txt', sep = '\t')
write.table(meta, file = 'meta.txt', sep = '\t')
write.table(count, file = 'count.txt', sep = '\t')
write.table(count_Filtered, file = 'count_Filtered.txt', sep = '\t')

if(NormMethod == "DESeq2"){
  write.table(CPM, file = 'count_Filtered_CPMnorm_combinedControl_LCM2_20210911.txt', sep = '\t')
  write.table(Norm, file = 'count_Filtered_CPMnorm_combinedControl_LCM2_20210911.txt', sep = '\t')
  write.table(log2Norm, file = 'count_Filtered_log2norm_combinedControl_LCM2_20210911.txt', sep = '\t')
  write.table(Norm_long, file = 'count_Filtered_CPMnorm_combinedControl_LCM2_20210911_long.txt', sep = '\t')
  write.table(log2Norm_long, file = 'count_Filtered_log2norm_combinedControl_LCM2_20210911_long.txt', sep = '\t')
}
if(NormMethod == "CPM"){
  write.table(Norm, file = 'count_Filtered_CPMnorm_ver6.txt', sep = '\t')
  write.table(log2Norm, file = 'count_Filtered_log2CPMnorm_BestLibSize_ver6_norm.txt', sep = '\t')
}



write.table(df_DEGs, file = paste0('DEGs_vs_', ControlSample, ".txt"), sep = "\t")
write.table(df_DEGs_filtered, file = 'DEGsfiltered_vs_correctControl_clean_combinedControl_20210911.txt', sep = "\t")
write.table(df_nrDEGs, file = paste0('CountDEGsfiltered_vs_', ControlSample, ".txt"), sep = "\t")

write.table(pca_ALL_top25_log2Norm, file = 'pca_ALLgenes_log2Norm_top25.txt', sep = '\t')
write.table(pca_ALL_top25_log2Norm_Genes, file = 'pca_Geneset_log2Norm_top25.txt', sep = '\t')
write.table(pca_ALL_top25genes_log2FC, file = 'pca_ALLgenes_log2FC_top25_vs_', ControlSample, '.txt', sep = '\t')
write.table(pca_ALL_top25genes_log2FC_Genes, file = 'pca_Geneset_log2FC_top25_vs_', ControlSample, '.txt', sep = '\t')

save.image(paste0(inputDir, "0_BioSpyder analysis pipeline_General_DataAssembly_LCM_Batch1.RData"))

log2NormFull <- log2Norm
log2NormFull$gene <- NULL
meta_Filtered_full <- meta_Filtered
for(i in 1:length(unique(meta_Filtered_full$CONCENTRATION))){
  DOSE <- unique(meta_Filtered_full$CONCENTRATION)[i]
  meta_Filtered <- meta_Filtered_full[meta_Filtered_full$CONCENTRATION == DOSE,]
  sample <- meta_Filtered$SAMPLE_ID
  log2Norm <- log2NormFull[, sample]

pdf(paste0("//VUW/Personal$/Homes/W/wijayals/Desktop/Experiments/2020/20200501_LCM_RealExperiment/HeatMap_allSample_allgenes_", DOSE,".pdf"), height =75, width = 50)
pheatmap(log2Norm, cluster_cols = TRUE, cluster_rows = TRUE,
         border_color = "grey90",clustering_distance_cols = "euclidean",
         clustering_distance_rows = "euclidean",clustering_method = "ward.D2",
         fontsize = 10, width = 10, height=10, legend = TRUE,
         fontsize_col = 10, show_rownames = T, show_colnames = T )
dev.off()
}

log2Norm$gene <- rownames(log2Norm)
