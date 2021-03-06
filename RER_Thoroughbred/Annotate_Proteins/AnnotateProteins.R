#' ### Description  
#' Given the proteins normalized spectra intensities from scaffold filter out
#' those proteins not expressed among all animals in the study. Annotate proteins
#' and save for downstream analysis. This is Melissa's RER project.  
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics  
#'   
#' Directory/File:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/AnnotateProteins.R  
#'   
#' **Input files:**  
#' Directory/Files:
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/Accession_Number_Report_KER_Glycogen_Marisa_20180817_merged.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/Quantitative_Samples_View_Report_KER_Glycogen_Marisa_20180817_merged.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/annotation_functions.Rdata  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/EquCabGenes.Rdata  
#' 
#' **Output files:**  
#'  * Directory/FIle:  
#'      /RER_Thoroughbred_Melissa_20180817/AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata  
#' 
#' **Render R**  
#'  * Directory/File  
#'     /RER_Thoroughbred_Melissa_20180817/AnnotateProteins/AnnotateProteins.qsub  
#'   
#' ***  
#' ### Environment Setup and Summary  
#' **Required Packages**  
library(rentrez)
library(biomaRt)

#' **Session Information**  
sessionInfo()

#' **Clear environment**  
rm(list=ls())

#' ### Load required R objects  

#' > Protein accession numbers from Scaffold
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/RER_Thoroughbred_Melissa_20180817"
ProtAnnot <- read.table(paste(dir, 
    "Accession_Number_Report_RER_Thoroughbred_Melissa_20180817_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
str(ProtAnnot)

#' > Normalized intesnsity values for each protein and animal from Scaffold 
NormInt <- read.table(paste(dir, 
    "Quantitative_Samples_View_Report_RER_Thoroughbred_Melissa_20180817_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")

# Standard deviations
StdDev <- NormInt[,grep("StdDev", colnames(NormInt))]
dim(StdDev)

# Normalized intensities
NormInt <- NormInt[,grep("StdDev", colnames(NormInt), invert=TRUE)]
NormInt$AccessionNumber <- sapply(strsplit(NormInt$AccessionNumber, split=" "), 
    function(x) x[[1]][1])
rownames(NormInt) <- NormInt$AccessionNumber
dim(NormInt)

#' > Load Functions
Ld <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/AnnotateProteins"
load(paste(Ld, "annotation_functions.Rdata", sep="/"))

#' > EquCab3 gene information
load(paste(Ld, "EquCabGenes.Rdata", sep="/"))


#' ### Annotate Protein  
#' Filter proteins with missing values
sum(is.na(NormInt))
ND <- rowSums(data.frame(apply(NormInt, 1, function(x) length(grep("No values", x))),
    apply(NormInt, 1, function(x) length(grep("No data", x)))))
summSD(ND[ND > 0])
length(ND[ND > 0])

#' Remove proteins with missing values
NormInt <- NormInt[ND == 0,]
StdDev <- StdDev[ND == 0,]
dim(NormInt)

#' Remove decoys from data
rmD <- c("HUMAN", "PIG", "HORSE", "DECOY")
idx <- unlist(lapply(rmD, function(x) grep(x, rownames(NormInt))))
Decoy <- NormInt[idx,]

# Number of Decoy
nrow(Decoy)

# Number of proteins left for downstream analysis
NormInt <- NormInt[-idx,]
StdDev <- StdDev[-idx,]
dim(NormInt)

#' Obtain Annotation Information
ProtID <- rownames(NormInt)

#' Run Annotation of ProtID in chunks of 10
seq1 <- seq(1, length(ProtID), 10)
seq2 <- c((seq1-1)[-1], length(ProtID))

rstAnnot <- list()
for (i in 1:length(seq1)){
    IDs <- ProtID[seq1[i]:seq2[i]]
    rstAnnot[[i]] <- AnnotProt(ids=IDs, annot=EquCab)
}


#' Merge annotation history for proteins expressed in Melissa's RER project
idx <- unlist(lapply(rstAnnot, function(x) names(x$AnnotHist)))
for(i in 1:length(rstAnnot)){
    AnnotHist <-    lapply(idx, function(x) rstAnnot[[i]]$AnnotHist[[x]])
}
names(AnnotHist) <- idx
length(AnnotHist)

#' Merge gene information
GeneInfo <- do.call(rbind, lapply(rstAnnot, function(x) x$GeneInfo))
rownames(GeneInfo) <- NULL
dim(GeneInfo)



#' ### Save protein information
save(NormInt, StdDev, Decoy, GeneInfo, AnnotHist, ProtAnnot,
    file=paste(getwd(), "RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata", sep="/"))


#' ### Run R Script
#+ eval = FALSE
htmlRunR
AnnotateProteins.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+RER Thoroughbred 20180817 Protein Annotation

