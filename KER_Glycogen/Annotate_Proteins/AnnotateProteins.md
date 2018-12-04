---
title: Glycogen 20180622 Protein Annotation
author: Deborah Velez-Irizarry
date: Tue Oct 23 10:28:57 EDT 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description  
Given the proteins normalized spectra intensities from scaffold filter out  
those proteins not expressed among all animals in the study. Annotate proteins  
and save for downstream analysis. This is Kennedy's glycogen project.  
 
***  
**Code:**  
Parent Directory:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622  
 
Directory/File:  

> &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/AnnotateProteins.R  

**Input files:**  
Directory/FIles:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/Accession_Number_Report_Glycogen_20180622_merged.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/Quantitative_Samples_View_Report_Glycogen_20180622_merged.txt  
 
**Output files:**  
Directory/File:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/annotation_functions.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/EquCabGenes.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/Glycogen_20180622_Protein_Annotation.Rdata  
 
**Render R**  
 * Directory/File  

> &nbsp;&nbsp;&nbsp;&nbsp;/AnnotateProteins/AnnotateProteins.qsub  
  
***  
### Environment Setup    
**Required Packages**  


```r
library(rentrez)
library(biomaRt)
```

**Session Information**


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.2.20-GCC-6.4.0-2.28/lib/libopenblas_haswellp-r0.2.20.so
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] biomaRt_2.36.1 rentrez_1.2.1  knitr_1.20    
## 
## loaded via a namespace (and not attached):
##  [1] Rcpp_0.12.17         AnnotationDbi_1.42.1 magrittr_1.5        
##  [4] hms_0.4.2            BiocGenerics_0.26.0  progress_1.2.0      
##  [7] IRanges_2.14.12      bit_1.1-14           R6_2.2.2            
## [10] rlang_0.2.1          httr_1.3.1           stringr_1.3.1       
## [13] blob_1.1.1           tools_3.5.1          parallel_3.5.1      
## [16] Biobase_2.40.0       DBI_1.0.0            assertthat_0.2.0    
## [19] bit64_0.9-7          digest_0.6.15        crayon_1.3.4        
## [22] S4Vectors_0.18.3     bitops_1.0-6         RCurl_1.95-4.10     
## [25] memoise_1.1.0        RSQLite_2.1.1        evaluate_0.10.1     
## [28] stringi_1.2.3        compiler_3.5.1       prettyunits_1.0.2   
## [31] stats4_3.5.1         XML_3.98-1.16        jsonlite_1.5        
## [34] pkgconfig_2.0.1
```

**Clear environment**


```r
rm(list=ls())
```

### Load required R objects  
> Load protein accession numbers from Scaffold


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622"
ProtAnnot <- read.table(paste(dir,
    "Accession_Number_Report_Glycogen_20180622_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
str(ProtAnnot)
```

```
## 'data.frame':	837 obs. of  2 variables:
##  $ Accession.Number: chr  "NP_001075228.1" "NP_001075247.1" "NP_001075249.1" "NP_001075250.1" ...
##  $ Protein.Name    : chr  "myosin-1 [Equus caballus]" "heat shock cognate 71 kDa protein [Equus caballus]" "lumican precursor [Equus caballus]" "elongation factor 1-alpha 1 [Equus caballus]" ...
```

> Load normalized intesnsity values for each protein and animal from Scaffold


```r
NormInt <- read.table(paste(dir,
    "Quantitative_Samples_View_Report_Glycogen_20180622_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")

# Standard deviations
StdDev <- NormInt[,grep("StdDev", colnames(NormInt))]
dim(StdDev)
```

```
## [1] 507  44
```

```r
# Normalized intensities
NormInt <- NormInt[,grep("StdDev", colnames(NormInt), invert=TRUE)]
dim(NormInt)
```

```
## [1] 507  47
```

```r
NormInt$AccessionNumber <- sapply(strsplit(NormInt$AccessionNumber, split=" "),
    function(x) x[[1]][1])
rownames(NormInt) <- NormInt$AccessionNumber
```

### Functions  
**Summary function**: includes standard deviation  
> `x` is a vector of numerals  
> `dec` is a scalar number indicating the decimal points to show in output


```r
summSD <- function(x, dec=3) round(c(summary(x),
    Std.Dev.=sd(x)), dec)[c("Min.", "1st Qu.", "Median", "Mean",
    "Std.Dev.", "3rd Qu.", "Max.")]
```

**Annotate Proteins Function:** Obtain Nucleotide IDs from Protein IDs 
> `ids` is a vector of protein ids
> `annot` is a data frame with annotation information for EquCap3  


```r
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
```

Save functions to Rdata file for further use


```r
save(summSD, AnnotProt, file=paste(getwd(), "annotation_functions.Rdata", sep="/"))
```

### Obtain Ensemble IDs  
 Ensemble Gene, Transcripts and Proteins IDs using entrez gene IDs  


```r
# Genes in Equ Cap with entrez gene IDs  
ensembl <- useDataset("ecaballus_gene_ensembl",mart=useMart("ensembl"))
ensembl <- useMart("ensembl",dataset="ecaballus_gene_ensembl")
attributes <- c("entrezgene", "vgnc","ensembl_gene_id",
    "ensembl_transcript_id", "ensembl_peptide_id")
EquCab <- getBM(attributes=attributes, filters="with_entrezgene",
    values=TRUE, mart=ensembl)
dim(EquCab)
```

```
## [1] 7869    5
```

```r
# Genes in Equ Cap with vgnc gene IDs
vgnc <- getBM(attributes=attributes, filters="with_vgnc",
    values=TRUE, mart=ensembl)
dim(vgnc)
```

```
## [1] 15255     5
```

Merge Gene list information


```r
EquCab <- rbind(EquCab, vgnc)
rownames(EquCab) <- NULL
dim(EquCab)
```

```
## [1] 23124     5
```

Obtain Gene information


```r
attributes <- c("chromosome_name", "start_position", "end_position", "strand", "external_gene_name",
    "entrezgene", "vgnc", "ensembl_gene_id", "ensembl_transcript_id", "ensembl_peptide_id",
    "transcript_start", "transcript_end", "transcription_start_site")
info <- getBM(attributes=attributes, values=EquCab$entrezgene, mart=ensembl)
```

Group transcripts per gene


```r
idx <- unique(info$ensembl_gene_id)
EquCabTrans <- lapply(idx, function(x) info[info$ensembl_gene_id %in% x,])
names(EquCabTrans) <- idx
```

Number of genes in EquCap3


```r
numb <- unlist(lapply(EquCabTrans, nrow))
table(numb)
```

```
## numb
##     1     2     3     4     5     6     7     8     9    10 
## 25198  1329   321    91    15    10     4    16     2     5
```

RYR2 transcripts


```r
EquCabTrans["ENSECAG00000006795"]
```

```
## $ENSECAG00000006795
##       chromosome_name start_position end_position strand
## 23731               1       73898255     74403258     -1
## 23732               1       73898255     74403258     -1
## 23733               1       73898255     74403258     -1
## 23734               1       73898255     74403258     -1
## 23735               1       73898255     74403258     -1
## 23736               1       73898255     74403258     -1
##       external_gene_name entrezgene       vgnc    ensembl_gene_id
## 23731               RYR2         NA VGNC:22643 ENSECAG00000006795
## 23732               RYR2         NA VGNC:22643 ENSECAG00000006795
## 23733               RYR2         NA VGNC:22643 ENSECAG00000006795
## 23734               RYR2         NA VGNC:22643 ENSECAG00000006795
## 23735               RYR2         NA VGNC:22643 ENSECAG00000006795
## 23736               RYR2         NA VGNC:22643 ENSECAG00000006795
##       ensembl_transcript_id ensembl_peptide_id transcript_start
## 23731    ENSECAT00000009674 ENSECAP00000007395         73936064
## 23732    ENSECAT00000009665 ENSECAP00000007386         73936064
## 23733    ENSECAT00000009649 ENSECAP00000007371         73936064
## 23734    ENSECAT00000009645 ENSECAP00000007367         73936064
## 23735    ENSECAT00000000088 ENSECAP00000000065         73898255
## 23736    ENSECAT00000029145 ENSECAP00000022768         73898258
##       transcript_end transcription_start_site
## 23731       74299190                 74299190
## 23732       74299190                 74299190
## 23733       74299190                 74299190
## 23734       74299190                 74299190
## 23735       73927986                 73927986
## 23736       74403258                 74403258
```

Retain only gene information for EquCabGenes


```r
EquCabGenes <- do.call(rbind, lapply(EquCabTrans, function(x) x[1, c(1:5)]))
dim(EquCabGenes)
```

```
## [1] 26991     5
```

Save EquCab Gene Annotation


```r
save(EquCab, EquCabGenes, EquCabTrans,  file=paste(getwd(), "EquCabGenes.Rdata", sep="/"))
```

### Annotate Protein  
Filter proteins with missing values


```r
sum(is.na(NormInt))
```

```
## [1] 0
```

```r
ND <- rowSums(data.frame(apply(NormInt, 1, function(x) length(grep("No values", x))),
    apply(NormInt, 1, function(x) length(grep("No data", x)))))
summSD(ND[ND > 0])
```

```
##     Min.  1st Qu.   Median     Mean Std.Dev.  3rd Qu.     Max. 
##     4.00    11.00    22.00    21.36    10.49    33.00    44.00
```

```r
length(ND[ND > 0])
```

```
## [1] 114
```

Remove proteins with missing values


```r
NormInt <- NormInt[ND == 0,]
StdDev <- StdDev[ND == 0,]
dim(NormInt)
```

```
## [1] 393  47
```

Remove decoys from data


```r
rmD <- c("HUMAN", "PIG", "HORSE", "DECOY")
idx <- unlist(lapply(rmD, function(x) grep(x, rownames(NormInt))))
Decoy <- NormInt[idx,]

# Number of Decoy
nrow(Decoy)
```

```
## [1] 6
```

```r
# Number of proteins left for downstream analysis
NormInt <- NormInt[-idx,]
StdDev <- StdDev[-idx,]
dim(NormInt)
```

```
## [1] 387  47
```

Obtain Annotation Information


```r
ProtID <- rownames(NormInt)
```

Run Annotation of ProtID in chunks of 10


```r
seq1 <- seq(1, length(ProtID), 10)
seq2 <- c((seq1-1)[-1], length(ProtID))

rstAnnot <- list()
for (i in 1:length(seq1)){
    IDs <- ProtID[seq1[i]:seq2[i]]
    rstAnnot[[i]] <- AnnotProt(ids=IDs, annot=EquCab)
}
```

Merge annotation history for proteins expressed in Kennedy's Glycogen project


```r
idx <- unlist(lapply(rstAnnot, function(x) names(x$AnnotHist)))
AnnotHist <- lapply(rstAnnot, function(x) x$AnnotHist)
AnnotHist <- unlist(AnnotHist, recursive=FALSE)
length(AnnotHist)
```

```
## [1] 384
```

Merge gene information


```r
GeneInfo <- do.call(rbind, lapply(rstAnnot, function(x) x$GeneInfo))
rownames(GeneInfo) <- NULL
dim(GeneInfo)
```

```
## [1] 581  14
```

Save protein information


```r
save(NormInt, StdDev, Decoy, GeneInfo, AnnotHist, ProtAnnot,
    file=paste(getwd(), "Glycogen_20180622_Protein_Annotation.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
AnnotateProteins.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G \
+Glycogen 20180622 Protein Annotation
```

