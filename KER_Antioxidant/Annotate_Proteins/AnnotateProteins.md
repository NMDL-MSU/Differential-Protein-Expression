# KER Glycogen 20180817 Protein Annotation  
*Deborah Velez-Irizarry*  
Date: October 8, 2018  
  
### Description  
Given the proteins normalized spectra intensities from scaffold filter out
those proteins not expressed among all animals in the study. Annotate proteins
and save for downstream analysis. This is Marisa's ROS project.  
  
***  
**Code:**  
 * Directory/File:  
   /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins/AnnotateProteins.R  
  
**Input files:**  
 * Directory:  
    /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817  
  
 * File:  
    Accession_Number_Report_KER_Glycogen_Marisa_20180817_merged.txt  
    Quantitative_Samples_View_Report_KER_Glycogen_Marisa_20180817_merged.txt  

 * Directory:  
    /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins  

 * Files:  
    annotation_functions.Rdata  
    EquCabGenes.Rdata  
 
**Output files:**  
 * Directory:  
    /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins  
 * File:  
    KER_Glycogen_20180817_Protein_Annotation.Rdata  

**Render R**  
 * Directory/File  
   /mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817/AnnotateProteins   
 
 * File:  
   AnnotateProteins.qsub  
***  
### Environment Setup and Summary  
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
Load protein accession numbers from Scaffold


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817"
ProtAnnot <- read.table(paste(dir, 
    "Accession_Number_Report_KER_Glycogen_Marisa_20180817_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")
str(ProtAnnot)
```

```
## 'data.frame':	711 obs. of  2 variables:
##  $ AccessionNumber: chr  "NP_001075247.1" "NP_001075249.1" "NP_001075254.1" "NP_001075257.1" ...
##  $ ProteinName    : chr  "heat shock cognate 71 kDa protein [Equus caballus]" "lumican precursor [Equus caballus]" "elongation factor 1-gamma [Equus caballus]" "band 3 anion transport protein [Equus caballus]" ...
```

Load normalized intesnsity values for each protein and animal from Scaffold 


```r
NormInt <- read.table(paste(dir, 
    "Quantitative_Samples_View_Report_KER_Glycogen_Marisa_20180817_merged.txt", sep="/"),
    header=TRUE, sep="\t", stringsAsFactors = FALSE, quote="")

# Standard deviations
StdDev <- NormInt[,grep("StdDev", colnames(NormInt))]
dim(StdDev)
```

```
## [1] 387  22
```

```r
# Normalized intensities
NormInt <- NormInt[,grep("StdDev", colnames(NormInt), invert=TRUE)]
NormInt$AccessionNumber <- sapply(strsplit(NormInt$AccessionNumber, split=" "), 
    function(x) x[[1]][1])
rownames(NormInt) <- NormInt$AccessionNumber
dim(NormInt)
```

```
## [1] 387  25
```

Load Functions


```r
Ld <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/AnnotateProteins"
load(paste(Ld, "annotation_functions.Rdata", sep="/"))
```

Load EquCab3 gene information


```r
load(paste(Ld, "EquCabGenes.Rdata", sep="/"))
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
##   11.000   11.000   11.000   11.379    2.043   11.000   22.000
```

```r
length(ND[ND > 0])
```

```
## [1] 29
```

Remove proteins with missing values


```r
NormInt <- NormInt[ND == 0,]
StdDev <- StdDev[ND == 0,]
dim(NormInt)
```

```
## [1] 358  25
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
## [1] 352  25
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

Merge annotation history for proteins expressed in Marisa's Glycogen project


```r
idx <- unlist(lapply(rstAnnot, function(x) names(x$AnnotHist)))
AnnotHist <- lapply(rstAnnot, function(x) x$AnnotHist)
AnnotHist <- unlist(AnnotHist, recursive=FALSE)
length(AnnotHist)
```

```
## [1] 349
```

Merge gene information


```r
GeneInfo <- do.call(rbind, lapply(rstAnnot, function(x) x$GeneInfo))
rownames(GeneInfo) <- NULL
dim(GeneInfo)
```

```
## [1] 526  14
```

Save protein information


```r
save(NormInt, StdDev, Decoy, GeneInfo, AnnotHist, ProtAnnot,
    file=paste(getwd(), "KER_Glycogen_20180817_Protein_Annotation.Rdata", sep="/"))
```

### Run R Script


```r
runR
AnnotateProteins.R nodes=1,cpus-per-task=1,time=03:00:00,mem=10G
```

