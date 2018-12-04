#' ### Description:  
#' Summarize differential protein expression analysis for KER glycogen study.  
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
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.R  
#'  
#' **Input files:**
#' Directory/File:
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/Glycogen_20180622_Protein_Annotation.Rdata  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Tables/*  
#' > &nbsp;&nbsp;&nbsp;EquCab3/EquCab3_Annotation.Rdata  
#'  
#' **Output files:**  
#'  
#' Directory/File:  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.Rdata  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Tables/*.txt  
#'  
#' Render R Script  
#'  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.qsub  
#'  
#' ***  

#' ### Code  
#' Clear Environment
rm(list=ls())

#' ### Load Data    
#' Load EquCab3 Annotation
load("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/EquCab3/EquCab3_Annotation.Rdata")

#' Load Proteins from Scaffold
load("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/AnnotateProteins/Glycogen_20180622_Protein_Annotation.Rdata")

#' Directory Files
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/"

#' Extract summary tables from Scaffolds output files
system(paste("bash ", dir, "Summary_Tables/extract_permutation_summary.sh", sep=""))

#' Summary Tables
fl <- list.files(paste(dir, "Summary_Tables", sep="/"))
fl <- fl[grep("txt", fl)]
PTsum <- lapply(fl, function(x) read.table(paste(dir, "Summary_Tables/", x, sep=""), 
    header=TRUE, sep="\t"))
names(PTsum) <- unlist(lapply(strsplit(fl, "[.]"), function(x) x[[1]][1]))

#' Remove temporary files from directory
system(paste("rm ", dir, "Summary_Tables/*.txt", sep=""))

#' Load Benjamini-Hochberg cutoff for a FDR < 0.05
pval.cut <- t(read.table(paste(dir, "Summary_Tables/Benjamini-Hochberg_Cutoff", sep=""), 
    header=FALSE, row.names=1))[1,]
pval.cut

#' Filter out proteins with low expression across samples
idx <- rownames(NormInt)
PTsum <- lapply(PTsum, function(x) do.call(rbind, lapply(idx, function(y) 
    x[grep(y, as.character(x$Protein)),])))

#' Extract Protein Information
ProtInfo <- do.call(rbind, lapply(idx, function(x) ProtEquCab3[as.character(ProtEquCab3$Genbank) %in% x,]))
head(ProtInfo)


#' ### Differentily expressed proteins
sig.PT <- lapply(names(pval.cut), function(x) PTsum[[x]][PTsum[[x]]$P.Value < pval.cut[x],])
names(sig.PT) <- names(pval.cut)
unlist(lapply(sig.PT, nrow))  

#' Create summary of differentially expressed genes with gene information
sig.Info <- lapply(sig.PT, function(x) do.call(rbind, lapply(1:nrow(x), function(y) 
    data.frame(x[y,1:4], FoldChange=ifelse(x[y,"Log2FoldChange"] > 0, 2^abs(x[y,"Log2FoldChange"]),  
        2^abs(x[y,"Log2FoldChange"]) * -1), 
    ProtInfo[as.character(ProtInfo$Genbank) %in% as.character(x[y,"Protein"]),-7]))))


#' ### Write results to file
system("mkdir Tables")
x <- lapply(names(sig.Info), function(x) write.table(sig.Info[x], 
    file=paste(getwd(), "/Tables/", x, ".txt", sep=""), 
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t"))

#' Save R Data 
save(sig.Info, pval.cut, file=paste(getwd(), "Summary_Results.Rdata", sep="/"))


#' ### Run R Script
#+ eval = FALSE
htmlRunR
Summary_Results.R nodes=1,cpus-per-task=1,time=02:00:00,mem=5G \
+KER Glycogen Differentially Expressed Proteins

