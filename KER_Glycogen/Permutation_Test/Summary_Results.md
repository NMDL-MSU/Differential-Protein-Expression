---
title: KER Glycogen Differentially Expressed Proteins
author: Deborah Velez-Irizarry
date: Tue Dec 4 11:16:28 EST 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Summarize differential protein expression analysis for KER glycogen study.  
Results of permutation test performed in Scaffold Q+.  
 
***  
 
**Code:**  
Parent Directory:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics  
 
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.R  
 
**Input files:**
Directory/File:
 
> &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/Glycogen_20180622_Protein_Annotation.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Tables/*  
> &nbsp;&nbsp;&nbsp;EquCab3/EquCab3_Annotation.Rdata  
 
**Output files:**  
 
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Tables/*.txt  
 
Render R Script  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.qsub  
 
***  
### Code  
Clear Environment


```r
rm(list=ls())
```

### Load Data    
Load EquCab3 Annotation


```r
load("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/EquCab3/EquCab3_Annotation.Rdata")
```

Load Proteins from Scaffold


```r
load("/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/AnnotateProteins/Glycogen_20180622_Protein_Annotation.Rdata")
```

Directory Files


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/Glycogen_Kennedy_20180622/AnnotateProteins/PermutationTest/"
```

Extract summary tables from Scaffolds output files


```r
system(paste("bash ", dir, "Summary_Tables/extract_permutation_summary.sh", sep=""))
```

Summary Tables


```r
fl <- list.files(paste(dir, "Summary_Tables", sep="/"))
fl <- fl[grep("txt", fl)]
PTsum <- lapply(fl, function(x) read.table(paste(dir, "Summary_Tables/", x, sep=""), 
    header=TRUE, sep="\t"))
names(PTsum) <- unlist(lapply(strsplit(fl, "[.]"), function(x) x[[1]][1]))
```

Remove temporary files from directory


```r
system(paste("rm ", dir, "Summary_Tables/*.txt", sep=""))
```

Load Benjamini-Hochberg cutoff for a FDR < 0.05


```r
pval.cut <- t(read.table(paste(dir, "Summary_Tables/Benjamini-Hochberg_Cutoff", sep=""), 
    header=FALSE, row.names=1))[1,]
pval.cut
```

```
##             Diet_Depl              Diet_Pre            Diet_Rep24 
##               0.00243               0.00092               0.00256 
##            Diet_Rep72          Fat_Depl-Pre          Fat_OverTime 
##               0.00046               0.00434               0.01150 
##        Fat_Rep24-Depl         Fat_Rep24-Pre        Fat_Rep72-Depl 
##               0.00099               0.00296               0.00420 
##         Fat_Rep72-Pre    SweetFeed_Depl-Pre    SweetFeed_OverTime 
##               0.00874               0.00013               0.00972 
## SweetFeed_Repl24-Depl  SweetFeed_Repl24-Pre SweetFeed_Repl72-Depl 
##               0.00079               0.00197               0.00499 
##  SweetFeed_Repl72-Pre 
##               0.00854
```

Filter out proteins with low expression across samples


```r
idx <- rownames(NormInt)
PTsum <- lapply(PTsum, function(x) do.call(rbind, lapply(idx, function(y) 
    x[grep(y, as.character(x$Protein)),])))
```

Extract Protein Information


```r
ProtInfo <- do.call(rbind, lapply(idx, function(x) ProtEquCab3[as.character(ProtEquCab3$Genbank) %in% x,]))
head(ProtInfo)
```

```
##       seqid    start      end strand  gene    GeneID        Genbank
## 31376 chr11 53325909 53350918      -  MYH1    791235 NP_001075228.1
## 43130 chr18 33649007 33862519      -   NEB 100057773 XP_023478641.1
## 32690 chr12 28604054 28614904      -  PYGM 100050345 NP_001138725.1
## 31068 chr11 49908755 49914315      -  ENO3 100061187 XP_005597080.1
## 26380 chr10 16076279 16085270      -   CKM 100065641 XP_001502572.1
## 34376 chr13 20695261 20701160      + ALDOA 100066121 XP_003362760.1
##                                           product
## 31376                                    myosin-1
## 43130                          nebulin isoform X1
## 32690       glycogen phosphorylase%2C muscle form
## 31068                     beta-enolase isoform X1
## 26380                      creatine kinase M-type
## 34376 fructose-bisphosphate aldolase A isoform X2
```

### Differentily expressed proteins


```r
sig.PT <- lapply(names(pval.cut), function(x) PTsum[[x]][PTsum[[x]]$P.Value < pval.cut[x],])
names(sig.PT) <- names(pval.cut)
unlist(lapply(sig.PT, nrow))  
```

```
##             Diet_Depl              Diet_Pre            Diet_Rep24 
##                    28                    10                    32 
##            Diet_Rep72          Fat_Depl-Pre          Fat_OverTime 
##                     6                    53                   139 
##        Fat_Rep24-Depl         Fat_Rep24-Pre        Fat_Rep72-Depl 
##                    10                    40                    55 
##         Fat_Rep72-Pre    SweetFeed_Depl-Pre    SweetFeed_OverTime 
##                   108                     1                   119 
## SweetFeed_Repl24-Depl  SweetFeed_Repl24-Pre SweetFeed_Repl72-Depl 
##                    10                    23                    62 
##  SweetFeed_Repl72-Pre 
##                   112
```

Create summary of differentially expressed genes with gene information


```r
sig.Info <- lapply(sig.PT, function(x) do.call(rbind, lapply(1:nrow(x), function(y) 
    data.frame(x[y,1:4], FoldChange=ifelse(x[y,"Log2FoldChange"] > 0, 2^abs(x[y,"Log2FoldChange"]),  
        2^abs(x[y,"Log2FoldChange"]) * -1), 
    ProtInfo[as.character(ProtInfo$Genbank) %in% as.character(x[y,"Protein"]),-7]))))
```

### Write results to file


```r
system("mkdir Tables")
x <- lapply(names(sig.Info), function(x) write.table(sig.Info[x], 
    file=paste(getwd(), "/Tables/", x, ".txt", sep=""), 
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t"))
```

Save R Data 


```r
save(sig.Info, pval.cut, file=paste(getwd(), "Summary_Results.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Summary_Results.R nodes=1,cpus-per-task=1,time=02:00:00,mem=5G \
+KER Glycogen Differentially Expressed Proteins
```

