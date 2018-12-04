#' # KER Glycogen 20180817 Protein Annotation  
#' *Deborah Velez-Irizarry*  
#' Date: October 8, 2018  
#'   
#' ### Description  
#' Given the proteins normalized spectra intensities from scaffold filter out
#' those proteins not expressed among all animals in the study. Annotate proteins
#' and save for downstream analysis. This is Marisa's ROS project.  
#'   
#' ***  
#' **Code:**  
#'  * Directory/File:  
#'    /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins/AnnotateProteins.R  
#'   
#' **Input files:**  
#'  * Directory:  
#'     /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817  
#'   
#'  * File:  
#'     Accession_Number_Report_KER_Glycogen_Marisa_20180817_merged.txt  
#'     Quantitative_Samples_View_Report_KER_Glycogen_Marisa_20180817_merged.txt  
#' 
#'  * Directory:  
#'     /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins  
#' 
#'  * Files:  
#'     annotation_functions.Rdata  
#'     EquCabGenes.Rdata  
#'  
#' **Output files:**  
#'  * Directory:  
#'     /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins  
#'  * File:  
#'     KER_Glycogen_20180817_Protein_Annotation.Rdata  
#' 
#' **Render R**  
#'  * Directory/File  
#'    /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins   
#'  
#'  * File:  
#'    AnnotateProteins.qsub  
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
#' Load protein accession numbers from Scaffold
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817"
ProtAnnot <- read.table(paste(dir, 
    "Accession_Number_Report_KER_Glycogen_Marisa_20180817_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
str(ProtAnnot)

#' Load normalized intesnsity values for each protein and animal from Scaffold 
NormInt <- read.table(paste(dir, 
    "Quantitative_Samples_View_Report_KER_Glycogen_Marisa_20180817_merged.txt", sep="/"),
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

#' Load Functions
Ld <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/AnnotateProteins"
load(paste(Ld, "annotation_functions.Rdata", sep="/"))

#' Load EquCab3 gene information
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


#' Merge annotation history for proteins expressed in Marisa's Glycogen project
idx <- unlist(lapply(rstAnnot, function(x) names(x$AnnotHist)))
AnnotHist <- lapply(rstAnnot, function(x) x$AnnotHist)
AnnotHist <- unlist(AnnotHist, recursive=FALSE)
length(AnnotHist)

#' Merge gene information
GeneInfo <- do.call(rbind, lapply(rstAnnot, function(x) x$GeneInfo))
rownames(GeneInfo) <- NULL
dim(GeneInfo)

#' Save protein information
save(NormInt, StdDev, Decoy, GeneInfo, AnnotHist, ProtAnnot,
    file=paste(getwd(), "KER_Glycogen_20180817_Protein_Annotation.Rdata", sep="/"))


#' ### Run R Script
#+ eval = FALSE
runR
AnnotateProteins.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G

