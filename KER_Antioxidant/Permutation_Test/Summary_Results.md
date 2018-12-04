---
title: KER Anntioxidant Pre-Excercise Differentially Expressed Proteins
author: Deborah Velez-Irizarry
date: Tue Dec 4 15:56:57 EST 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Summarize differential protein expression analysis for KER antioxidant study.  
Results of permutation test performed in Scaffold Q+.  
 
***  
 
**Code:**  
Parent Directory:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics  
 
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/KER_Marisa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.R  
 
**Input files:**
Directory/File:
 
> &nbsp;&nbsp;&nbsp;&nbsp;/KER_Marisa_20180817/AnnotateProteins/KER_Glycogen_20180817_Protein_Annotation.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/KER_Marisa_20180817/SportvsRaceHorsesMatrix.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/KER_Marisa_20180817/AnnotateProteins/PermutationTest/Marisa_20180817_scaffold_permutation_test.txt  
> &nbsp;&nbsp;&nbsp;EquCab3/EquCab3_Annotation.Rdata  
 
**Output files:**  
 
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/KER_Marisa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results_KER_antioxidant.Rdata  
> &nbsp;&nbsp;&nbsp;&nbsp;/KER_Marisa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results_KER_antioxidant.txt  
 
Render R Script  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/KER_Marisa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.qsub  
 
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
dir="/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/KER_Marisa_20180817"
load(paste(dir, "AnnotateProteins/KER_Glycogen_20180817_Protein_Annotation.Rdata", sep="/"))
```

Load Animal Information


```r
anim <- read.table(paste(dir, "SportvsRaceHorsesMatrix.txt", sep="/"), header=TRUE, sep="\t")
head(anim)
```

```
##        Animal  Type Treatment    Assay    Pre X10m.Post X1hr.Post
## 1      Bullet Sport   Control Cysteine 41.856    12.863    78.523
## 2   Easy_Lion Sport   Control Cysteine 15.037    12.796    10.710
## 3 Fair_Dinkum Sport   Control Cysteine 32.644    50.917    16.811
## 4       Jamie Sport   Control Cysteine 18.935    19.919    12.842
## 5      Stoney Sport   Control Cysteine 97.527    12.074    12.250
## 6      Stella Sport   Control Cysteine 23.053    70.292    67.699
##   X4hr.Post
## 1    73.601
## 2    10.655
## 3    40.173
## 4    44.411
## 5    15.191
## 6    62.462
```

Extract summary tables from Scaffolds output files


```r
dir <- paste(dir, "/AnnotateProteins/PermutationTest/", sep="")
system(paste("bash ", dir, "extract_permutation_summary.sh", sep=""))
system(paste("mv ", dir, "/Marisa_20180817_scaffold_permutation_test.txt ", getwd(), sep=""))

# Summary Tables
PTsum <- read.table("Marisa_20180817_scaffold_permutation_test.txt", header=TRUE, sep="\t")
dim(PTsum)
```

```
## [1] 387  48
```

Benjamini-Hochberg cutoff for a FDR < 0.05


```r
pval.cut <- 0.00455
```

Filter out proteins with low expression across samples


```r
idx <- rownames(NormInt)
PTsum <- lapply(idx, function(x) data.frame(ProtEquCab3[as.character(ProtEquCab3$Genbank) %in% x,],
    PTsum[grep(x, as.character(PTsum$Protein)), 2:4]))
PTsum <- do.call(rbind, PTsum)
```

### Differentily expressed proteins


```r
sig.PT <-  PTsum[PTsum[,"P.Value"] < pval.cut,]
sig.PT <- data.frame(sig.PT, FoldChange=ifelse(sig.PT$Log2FoldChange > 0, 2^sig.PT$Log2FoldChange,
    (2^abs(sig.PT$Log2FoldChange) * -1)))
sig.PT <- sig.PT[order(sig.PT$FoldChange, decreasing = TRUE),]
dim(sig.PT)
```

```
## [1] 40 12
```

Up-regulated proteins  


```r
# Proteins
idx <- sig.PT$FoldChange > 0
sig.PT[idx, c(1,5,12,8)]
```

```
##       seqid         gene FoldChange
## 9666   chr3          ALB   1.301342
## 19429  chr7        APOA1   1.265757
## 19310  chr7        CRYAB   1.214195
## 5763   chr2        FABP3   1.205808
## 55592 chr28           MB   1.164734
## 16785  chr6 LOC100061692   1.164734
## 2850   chr1        COX5A   1.132884
## 1873   chr1        ACTN2   1.117287
## 1642   chr1        VDAC2   1.101905
## 8965   chr3       COX4I1   1.094294
## 11875  chr4         FLNC   1.086735
## 2329   chr1         IDH2   1.086735
## 9288   chr3       PDLIM5   1.086735
## 55445 chr28       MYBPC1   1.079228
## 55895 chr28         ACO2   1.079228
## 33984 chr13         MDH2   1.079228
## 56979 chr30           FH   1.079228
## 31165 chr11       ACADVL   1.079228
## 18092  chr6      ATP5F1B   1.071773
## 54858 chr27        ACSL1   1.071773
## 10924  chr4         OGDH   1.064370
## 47977 chr21          NNT   1.049717
##                                                                           product
## 9666                                                      serum albumin precursor
## 19429                                                          apolipoprotein A-I
## 19310                                                    alpha-crystallin B chain
## 5763                                          fatty acid-binding protein%2C heart
## 55592                                                                   myoglobin
## 16785                                                       alpha-2-macroglobulin
## 2850                             cytochrome c oxidase subunit 5A%2C mitochondrial
## 1873                                                   alpha-actinin-2 isoform X1
## 1642                          voltage-dependent anion-selective channel protein 2
## 8965                    cytochrome c oxidase subunit 4 isoform 1%2C mitochondrial
## 11875                                                        filamin-C isoform X1
## 2329                             isocitrate dehydrogenase [NADP]%2C mitochondrial
## 9288                                      PDZ and LIM domain protein 5 isoform X9
## 55445                            myosin-binding protein C%2C slow-type isoform X2
## 55895                                        aconitate hydratase%2C mitochondrial
## 33984                                       malate dehydrogenase%2C mitochondrial
## 56979                                         fumarate hydratase%2C mitochondrial
## 31165 very long-chain specific acyl-CoA dehydrogenase%2C mitochondrial isoform X6
## 18092                                  ATP synthase subunit beta%2C mitochondrial
## 54858                              long-chain-fatty-acid--CoA ligase 1 isoform X2
## 10924                    2-oxoglutarate dehydrogenase%2C mitochondrial isoform X3
## 47977                         NAD(P) transhydrogenase%2C mitochondrial isoform X1
```

```r
# Number
sum(idx)
```

```
## [1] 22
```

Down-regulated proteins  


```r
# Proteins
idx <- sig.PT$FoldChange < 0
sig.PT[idx, c(1,5,12,8)]
```

```
##       seqid   gene FoldChange
## 17391  chr6   PFKM  -1.049717
## 23382  chr8  MYOM1  -1.057018
## 34975 chr13    SRL  -1.064370
## 3016   chr1    PKM  -1.071773
## 21697  chr7   LDHA  -1.071773
## 31377 chr11   MYH1  -1.079228
## 49924 chr23   FBP2  -1.079228
## 55014 chr27  MYOM2  -1.086735
## 35998 chr14  ANXA6  -1.094294
## 25560 chr10    GPI  -1.101905
## 58900  chrX  PHKA1  -1.101905
## 8088   chr3   PHKB  -1.101905
## 31069 chr11   ENO3  -1.109569
## 15482  chr5   PGM1  -1.109569
## 26829 chr10 MYBPC2  -1.125058
## 34431 chr13 ATP2A1  -1.125058
## 48868 chr22  MYLK2  -1.148698
## 57731 chr31  ARMT1  -1.319508
##                                                                                     product
## 17391                         ATP-dependent 6-phosphofructokinase%2C muscle type isoform X3
## 23382                                                                 myomesin-1 isoform X4
## 34975                                                               sarcalumenin isoform X1
## 3016                                                         pyruvate kinase PKM isoform M1
## 21697                                                       L-lactate dehydrogenase A chain
## 31377                                                                   myosin-1 isoform X1
## 49924                                               fructose-1%2C6-bisphosphatase isozyme 2
## 55014                                                                 myomesin-2 isoform X3
## 35998                                                                 annexin A6 isoform X1
## 25560                                                         glucose-6-phosphate isomerase
## 58900 phosphorylase b kinase regulatory subunit alpha%2C skeletal muscle isoform isoform X1
## 8088                                         phosphorylase b kinase regulatory subunit beta
## 31069                                                                          beta-enolase
## 15482                                                       phosphoglucomutase-1 isoform X2
## 26829                                                 myosin-binding protein C%2C fast-type
## 34431                        sarcoplasmic/endoplasmic reticulum calcium ATPase 1 isoform X1
## 48868                                myosin light chain kinase 2%2C skeletal/cardiac muscle
## 57731                                                 protein-glutamate O-methyltransferase
```

```r
# Number
sum(idx)
```

```
## [1] 18
```

### Write results to file


```r
write.table(sig.PT, 
    file=paste(getwd(), "Summary_Results_KER_antioxidant.txt", sep="/"), 
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
```

Save R Data 


```r
save(anim, sig.PT, pval.cut, file=paste(getwd(), "Summary_Results_KER_antioxidant.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Summary_Results.R nodes=1,cpus-per-task=1,time=02:00:00,mem=5G \
+KER Anntioxidant Pre-Excercise Differentially Expressed Proteins
```

