
#' ### Description:  
#' Summarize differential protein expression analysis for RER Thoroughbred study.  
#' Results of permutation test performed in Scaffold Q+.  
#'  
#' ***  
#'  
#' **Code:**  
#' Parent Directory:  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics  
#'  
#' Directory/File:  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.R  
#'  
#' **Input files:**
#' Directory/File:
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata    
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/RER_Thoroughbred_Animal_Information.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/RER_Thoroughbred_Melissa_20180817_scaffold_permutation_test.txt  
#' > &nbsp;&nbsp;&nbsp;EquCab3/EquCab3_Annotation.Rdata  
#'  
#' **Output files:**  
#'  
#' Directory:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results
#'   
#' Files:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;Summary_Results_RER_Thoroughbred.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;Summary_Results_RER_Thoroughbred.Rdata  
#'  
#' Render R Script  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.qsub  
#'  
#' ***  

#' ### Code  
#' Clear Environment
rm(list=ls())

#' ### Load Data    
#' Load EquCab3 Annotation
load("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/EquCab3/EquCab3_Annotation.Rdata")

#' Load Proteins from Scaffold
dir="/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/RER_Thoroughbred_Melissa_20180817"
load(paste(dir, "AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata", sep="/"))

#' Load Animal Information
anim <- read.table(paste(dir, "RER_Thoroughbred_Animal_Information.txt", sep="/"), header=TRUE, sep="\t")
head(anim)

#' Extract summary tables from Scaffolds output files
dir <- paste(dir, "/AnnotateProteins/PermutationTest/", sep="")
system(paste("bash ", dir, "extract_permutation_summary.sh", sep=""))
system(paste("mv ", dir, "/RER_Thoroughbred_Melissa_20180817_scaffold_permutation_test.txt ", getwd(), sep=""))

# Summary Tables
PTsum <- read.table("RER_Thoroughbred_Melissa_20180817_scaffold_permutation_test.txt", header=TRUE, sep="\t")
dim(PTsum)

#' Benjamini-Hochberg cutoff for a FDR < 0.05
pval.cut <- 0.01333

#' Filter out proteins with low expression across samples
idx <- rownames(NormInt)
PTsum <- lapply(idx, function(x) data.frame(ProtEquCab3[as.character(ProtEquCab3$Genbank) %in% x,],
    PTsum[grep(x, as.character(PTsum$Protein)), 2:4]))
PTsum <- do.call(rbind, PTsum)

#' ### Differentily expressed proteins
sig.PT <-  PTsum[PTsum[,"P.Value"] < pval.cut,]
sig.PT <- data.frame(sig.PT, FoldChange=ifelse(sig.PT$Log2FoldChange > 0, 2^sig.PT$Log2FoldChange,
    (2^abs(sig.PT$Log2FoldChange) * -1)))
sig.PT <- sig.PT[order(sig.PT$FoldChange, decreasing = TRUE),]
dim(sig.PT)

#' Up-regulated proteins  
# Proteins
idx <- sig.PT$FoldChange > 0
sig.PT[idx, c(1,5,12,8)]

# Number
sum(idx)

#' Down-regulated proteins  
# Proteins
idx <- sig.PT$FoldChange < 0
sig.PT[idx, c(1,5,12,8)]

# Number
sum(idx)


#' ### Write results to file
write.table(sig.PT, 
    file=paste(getwd(), "Summary_Results_RER_Thoroughbred.txt", sep="/"), 
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")

#' Save R Data 
save(anim, sig.PT, pval.cut, file=paste(getwd(), "Summary_Results_RER_Thoroughbred.Rdata", sep="/"))


#' ### Run R Script
#+ eval = FALSE
htmlRunR
Summary_Results.R nodes=1,cpus-per-task=1,time=02:00:00,mem=5G \
+RER Thoroughbred Differentially Expressed Proteins

