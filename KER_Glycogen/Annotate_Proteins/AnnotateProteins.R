#' ### Description  
#' Given the proteins normalized spectra intensities from scaffold filter out  
#' those proteins not expressed among all animals in the study. Annotate proteins  
#' and save for downstream analysis. This is Kennedy's glycogen project.  
#'  
#' ***  
#' **Code:**  
#' Parent Directory:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622  
#'  
#' Directory/File:  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/AnnotateProteins.R  
#' 
#' **Input files:**  
#' Directory/FIles:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Accession_Number_Report_Glycogen_20180622_merged.txt  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/Quantitative_Samples_View_Report_Glycogen_20180622_merged.txt  
#'  
#' **Output files:**  
#' Directory/File:  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/annotation_functions.Rdata  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/EquCabGenes.Rdata  
#' > &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/Glycogen_20180622_Protein_Annotation.Rdata  
#'  
#' **Render R**  
#'  * Directory/File  
#' 
#' > &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/AnnotateProteins.qsub  
#'   
#' ***  

#' ### Environment Setup    
#' **Required Packages**  
library(rentrez)
library(biomaRt)

#' **Session Information**
sessionInfo()

#' **Clear environment**
rm(list=ls())

#' ### Load required R objects  

#' > Load protein accession numbers from Scaffold
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622"
ProtAnnot <- read.table(paste(dir,
    "Accession_Number_Report_Glycogen_20180622_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
str(ProtAnnot)

#' > Load normalized intesnsity values for each protein and animal from Scaffold
NormInt <- read.table(paste(dir,
    "Quantitative_Samples_View_Report_Glycogen_20180622_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")

# Standard deviations
StdDev <- NormInt[,grep("StdDev", colnames(NormInt))]
dim(StdDev)

# Normalized intensities
NormInt <- NormInt[,grep("StdDev", colnames(NormInt), invert=TRUE)]
dim(NormInt)
NormInt$AccessionNumber <- sapply(strsplit(NormInt$AccessionNumber, split=" "),
    function(x) x[[1]][1])
rownames(NormInt) <- NormInt$AccessionNumber


#' ### Functions  

#' **Summary function**: includes standard deviation  

#' > `x` is a vector of numerals  
#' > `dec` is a scalar number indicating the decimal points to show in output
summSD <- function(x, dec=3) round(c(summary(x),
    Std.Dev.=sd(x)), dec)[c("Min.", "1st Qu.", "Median", "Mean",
    "Std.Dev.", "3rd Qu.", "Max.")]


#' **Annotate Proteins Function:** Obtain Nucleotide IDs from Protein IDs 

#' > `ids` is a vector of protein ids
#' > `annot` is a data frame with annotation information for EquCap3  
AnnotProt <- function(ids, annot){
    # Obtain the links available for transcripts given the protein accession number
    links <- entrez_link(dbfrom="protein", db="nuccore",
        id=ids, by_id = TRUE)
    linkids <- sapply(links, function(x) x$links$protein_nuccore_mrna)
    idxN <- sapply(links, function(x) length(x$links$protein_nuccore_mrna) > 0)

    # Obtain available transcript ids
    linkNuc <- entrez_summary(id = linkids[idxN], db = "nuccore")
    NucID <- sapply(linkNuc, function(x) strsplit(as.character(x["extra"]),
            split="[|]")[[1]][4])
    names(NucID) <- ids[idxN]

    # Obtain gene IDs from protein IDs
    links <- entrez_link(dbfrom="protein", db="gene",
        id=ids, by_id = TRUE)
    linkids <- sapply(links, function(x) x$links$protein_gene)
    idxG <- sapply(links, function(x) length(x$links$protein_gene) > 0)

    # Stop function if no information was retrieved
    if(sum(idxN, idxG) == 0){
        print("No information obtained for:", ids)
        return(NA)
    }

    # Obtain gene information
    linkGene <- entrez_summary(id = linkids[idxG], db = "gene")
    pos <- do.call(rbind, lapply(linkGene,
        function(x) ifelse(is.null(nrow(x$genomicinfo)),
            return(t(data.frame(NA, c("chrloc", "chraccver",
                "chrstart", "chrstop", "exoncount"), row.names=2))),
            return(x$genomicinfo))))
    GeneInfo <- cbind(GeneID=unlist(lapply(linkGene, function(x) x$uid)),
            pos, Organism=unlist(lapply(linkGene, function(x) x$organism$scientificname)))
    rownames(GeneInfo) <- ids[idxG]

    # Add missing proteins to GeneInfo
    nMis <- sum(!ids %in% rownames(GeneInfo))
        if(nMis > 0){
        df <- t(data.frame(NA, colnames(GeneInfo), row.names=2))
        df <- df[rep(seq_len(nrow(df)), each=nMis),]
        GeneInfo <- rbind(GeneInfo, df)
        GeneInfo <- GeneInfo[ids,]
        rownames(GeneInfo) <- ids
    }

    # Merge GeneInfo with NucID
    GeneInfo <- data.frame(ids=rownames(GeneInfo), GeneInfo, NucID=NucID[ids])

    # Save previous annotation for each gene
    AnnotHist <- lapply(linkGene, function(x) x$locationhist)
    names(linkGene) <- ids[idxG]

    # Obtain ensemble IDs
    ensem <- lapply(as.character(GeneInfo$GeneID), function(x)
        ifelse(is.na(x) | !x %in% annot$entrezgene,
            return(t(data.frame(NA, colnames(annot), row.names=2))),
            return(annot[as.character(annot$entrezgene) %in% x,])))
    names(ensem) <- ids

    # Check for entrez gene IDs with more than one ensemble ID
    ck <- unlist(lapply(ensem, nrow))
    idx <- ck > 1

    # Repeat rows in GeneInfo
    if(sum(idx) > 0){
        df <- do.call(rbind, lapply(names(ck)[idx], function(x)
            GeneInfo[x,][rep(seq_len(nrow(GeneInfo[x,])), each=ck[x]-1),]))
        GeneInfo <- rbind(GeneInfo, df)
        rownames(GeneInfo) <- NULL
        # Order
        GeneInfo <- do.call(rbind, lapply(ids, function(x)
            GeneInfo[as.character(GeneInfo$ids) %in% x,]))
    }

    # Merge GeneInfo with ensemble IDs
    GeneInfo <- cbind(GeneInfo, do.call(rbind, ensem))

    # Return gene information and annotation history in a list
    rst <- list(GeneInfo=GeneInfo, AnnotHist=AnnotHist)
    return(rst)
}

#' Save functions to Rdata file for further use
save(summSD, AnnotProt, file=paste(getwd(), "annotation_functions.Rdata", sep="/"))


#' ### Obtain Ensemble IDs  
#'  Ensemble Gene, Transcripts and Proteins IDs using entrez gene IDs  
# Genes in Equ Cap with entrez gene IDs  
ensembl <- useDataset("ecaballus_gene_ensembl",mart=useMart("ensembl"))
ensembl <- useMart("ensembl",dataset="ecaballus_gene_ensembl")
attributes <- c("entrezgene", "vgnc","ensembl_gene_id",
    "ensembl_transcript_id", "ensembl_peptide_id")
EquCab <- getBM(attributes=attributes, filters="with_entrezgene",
    values=TRUE, mart=ensembl)
dim(EquCab)

# Genes in Equ Cap with vgnc gene IDs
vgnc <- getBM(attributes=attributes, filters="with_vgnc",
    values=TRUE, mart=ensembl)
dim(vgnc)

#' Merge Gene list information
EquCab <- rbind(EquCab, vgnc)
rownames(EquCab) <- NULL
dim(EquCab)

#' Obtain Gene information
attributes <- c("chromosome_name", "start_position", "end_position", "strand", "external_gene_name",
    "entrezgene", "vgnc", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id",
    "transcript_start", "transcript_end", "transcription_start_site")
info <- getBM(attributes=attributes, values=EquCab$entrezgene, mart=ensembl)

#' Group transcripts per gene
idx <- unique(info$ensembl_gene_id)
EquCabTrans <- lapply(idx, function(x) info[info$ensembl_gene_id %in% x,])
names(EquCabTrans) <- idx

#' Number of genes in EquCap3
numb <- unlist(lapply(EquCabTrans, nrow))
table(numb)

#' RYR2 transcripts
EquCabTrans["ENSECAG00000006795"]

#' Retain only gene information for EquCabGenes
EquCabGenes <- do.call(rbind, lapply(EquCabTrans, function(x) x[1, c(1:5)]))
dim(EquCabGenes)


#' Save EquCab Gene Annotation
save(EquCab, EquCabGenes, EquCabTrans,  file=paste(getwd(), "EquCabGenes.Rdata", sep="/"))


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


#' Merge annotation history for proteins expressed in Kennedy's Glycogen project
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
    file=paste(getwd(), "Glycogen_20180622_Protein_Annotation.Rdata", sep="/"))


#' ### Run R Script
#+ eval = FALSE
htmlRunR
AnnotateProteins.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Glycogen 20180622 Protein Annotation

