---
title: RER Thoroughbred Differentially Expressed Proteins
author: Deborah Velez-Irizarry
date: Tue Dec 4 16:34:54 EST 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---
### Description:  
Summarize differential protein expression analysis for RER Thoroughbred study.  
Results of permutation test performed in Scaffold Q+.  
 
***  
 
**Code:**  
Parent Directory:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics  
 
Directory/File:  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.R  
 
**Input files:**
Directory/File:
 
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata    
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/RER_Thoroughbred_Animal_Information.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/RER_Thoroughbred_Melissa_20180817_scaffold_permutation_test.txt  
> &nbsp;&nbsp;&nbsp;EquCab3/EquCab3_Annotation.Rdata  
 
**Output files:**  
 
Directory:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results
  
Files:  
  
> &nbsp;&nbsp;&nbsp;&nbsp;Summary_Results_RER_Thoroughbred.txt  
> &nbsp;&nbsp;&nbsp;&nbsp;Summary_Results_RER_Thoroughbred.Rdata  
 
Render R Script  
 
> &nbsp;&nbsp;&nbsp;&nbsp;/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results.qsub  
 
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
dir="/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred/Proteomics/RER_Thoroughbred_Melissa_20180817"
load(paste(dir, "AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata", sep="/"))
```

Load Animal Information


```r
anim <- read.table(paste(dir, "RER_Thoroughbred_Animal_Information.txt", sep="/"), header=TRUE, sep="\t")
head(anim)
```

```
##   X    ID Age   Sex Last_Work AST   CK Biopsy_Changes  Dx Owner_Vet
## 1 1 12401   7 FALSE       2.0 921  663              2 RER Slaughter
## 2 2 12402   3 FALSE       3.0 448  521              2 RER     Tores
## 3 3 12403   3 FALSE       3.0 280  181              3 RER Slaughter
## 4 4 12610   3 FALSE       1.5 381  262              0 RER Slaughter
## 5 5 12611   3 FALSE       2.0 654 1086              1 RER Slaughter
## 6 6 12612   2 FALSE       2.3 503  543              0 RER Slaughter
##              Tx               Name
## 1 no_dantrolene   DeterminedYankee
## 2 no_dantrolene ChoppyChoppyChoppy
## 3 no_dantrolene            Basheba
## 4          <NA>     PromiseOfPeace
## 5    dantrolene            LaDonia
## 6 no_dantrolene             LaDama
```

Extract summary tables from Scaffolds output files


```r
dir <- paste(dir, "/AnnotateProteins/PermutationTest/", sep="")
system(paste("bash ", dir, "extract_permutation_summary.sh", sep=""))
system(paste("mv ", dir, "/RER_Thoroughbred_Melissa_20180817_scaffold_permutation_test.txt ", getwd(), sep=""))

# Summary Tables
PTsum <- read.table("RER_Thoroughbred_Melissa_20180817_scaffold_permutation_test.txt", header=TRUE, sep="\t")
dim(PTsum)
```

```
## [1] 401  48
```

Benjamini-Hochberg cutoff for a FDR < 0.05


```r
pval.cut <- 0.01333
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
## [1] 125  12
```

Up-regulated proteins  


```r
# Proteins
idx <- sig.PT$FoldChange > 0
sig.PT[idx, c(1,5,12,8)]
```

```
##       seqid         gene FoldChange
## 35432 chr13          HBA   3.031433
## 24122  chr9          CA1   2.928171
## 21169  chr7          HBB   2.566852
## 8651   chr3 LOC100067869   2.234574
## 19429  chr7        APOA1   2.158456
## 16785  chr6 LOC100061692   2.143547
## 24120  chr9          CA2   2.143547
## 41555 chr16           CP   2.128740
## 20125  chr7        PRDX2   2.114036
## 41345 chr16           TF   2.099433
## 9666   chr3          ALB   1.958841
## 29566 chr11       SLC4A1   1.958841
## 22173  chr8 LOC106783470   1.741101
## 51643 chr24 LOC100065068   1.693491
## 44665 chr19         AHSG   1.658639
## 8198   chr3 LOC100050685   1.602140
## 9682   chr3           GC   1.569168
## 51204 chr24 LOC100147142   1.558329
## 13675  chr5        SPTA1   1.547565
## 26041 chr10        BLVRB   1.536875
## 50970 chr24         SPTB   1.505247
## 27581 chr10         A1BG   1.474269
## 44669 chr19         KNG1   1.356604
## 46343 chr20 LOC100059239   1.248331
## 31472 chr11          UBC   1.197479
## 13612  chr5       ATP1A2   1.197479
## 24506  chr9       UBE2V2   1.197479
## 44673 chr19       ADIPOQ   1.189207
## 54421 chr27         ANK1   1.180993
## 38581 chr15        CALM2   1.148698
## 2725   chr1        PSMA4   1.148698
## 14240  chr5     SELENBP1   1.140764
## 32772 chr12        CAPN1   1.140764
## 26665 chr10          HRC   1.140764
## 10248  chr3         QDPR   1.132884
## 58305  chrX          DMD   1.117287
## 5144   chr2        PRDX1   1.109569
## 52618 chr25          GSN   1.109569
## 51959 chr25          VCP   1.094294
## 8207   chr3 LOC100051065   1.086735
## 44670 chr19       EIF4A2   1.086735
## 24933  chr9        SNTB1   1.079228
## 34975 chr13          SRL   1.071773
## 46706 chr20         GLO1   1.071773
## 57770 chr31        PCMT1   1.071773
## 16635  chr6        GAPDH   1.064370
## 13604  chr5        CASQ1   1.064370
## 55857 chr28         ST13   1.057018
## 14938  chr5          AGL   1.049717
## 25913 chr10         RYR1   1.049717
## 35998 chr14        ANXA6   1.035265
## 58900  chrX        PHKA1   1.035265
##                                                                                     product
## 35432                                                              hemoglobin subunit alpha
## 24122                                                                  carbonic anhydrase 1
## 21169                                                               hemoglobin subunit beta
## 8651                                                                            haptoglobin
## 19429                                                                    apolipoprotein A-I
## 16785                                                                 alpha-2-macroglobulin
## 24120                                                                  carbonic anhydrase 2
## 41555                                                              ceruloplasmin isoform X1
## 20125                                                            peroxiredoxin-2 isoform X1
## 41345                                                             serotransferrin precursor
## 9666                                                                serum albumin precursor
## 29566                                                        band 3 anion transport protein
## 22173                                                                   LOW QUALITY PROTEIN
## 51643                                               alpha-1-antiproteinase 2-like precursor
## 44665                                                               alpha-2-HS-glycoprotein
## 8198                                                                 liver carboxylesterase
## 9682                                                              vitamin D-binding protein
## 51204                                                                   LOW QUALITY PROTEIN
## 13675                                                spectrin alpha chain%2C erythrocytic 1
## 26041                                                              flavin reductase (NADPH)
## 50970                                        spectrin beta chain%2C erythrocytic isoform X1
## 27581                                                                 alpha-1B-glycoprotein
## 44669                                                                kininogen-1 isoform X2
## 46343                                                                       complement C4-A
## 31472                                                                           ubiquitin C
## 13612                                  sodium/potassium-transporting ATPase subunit alpha-2
## 24506                                             ubiquitin-conjugating enzyme E2 variant 2
## 44673                                                                           adiponectin
## 54421                                                                  ankyrin-1 isoform X4
## 38581                                                                            calmodulin
## 2725                                             proteasome subunit alpha type-4 isoform X1
## 14240                                                            selenium-binding protein 1
## 32772                                                           calpain-1 catalytic subunit
## 26665                         sarcoplasmic reticulum histidine-rich calcium-binding protein
## 10248                                                 dihydropteridine reductase isoform X1
## 58305                                                                dystrophin isoform X11
## 5144                                                                        peroxiredoxin-1
## 52618                                                                              gelsolin
## 51959                                             transitional endoplasmic reticulum ATPase
## 8207                                                      liver carboxylesterase isoform X1
## 44670                                                    eukaryotic initiation factor 4A-II
## 24933                                                          beta-1-syntrophin isoform X2
## 34975                                                               sarcalumenin isoform X1
## 46706                                                              lactoylglutathione lyase
## 57770                    protein-L-isoaspartate(D-aspartate) O-methyltransferase isoform X2
## 16635                                              glyceraldehyde-3-phosphate dehydrogenase
## 13604                                                                       calsequestrin-1
## 55857                                                             hsc70-interacting protein
## 14938                                                glycogen debranching enzyme isoform X1
## 25913                                                       ryanodine receptor 1 isoform X1
## 35998                                                                 annexin A6 isoform X1
## 58900 phosphorylase b kinase regulatory subunit alpha%2C skeletal muscle isoform isoform X1
```

```r
# Number
sum(idx)
```

```
## [1] 52
```

Down-regulated proteins  


```r
# Proteins
idx <- sig.PT$FoldChange < 0
sig.PT[idx, c(1,5,12,8)]
```

```
##             seqid         gene FoldChange
## 11875        chr4         FLNC  -1.035265
## 16706        chr6         PHB2  -1.035265
## 3016         chr1          PKM  -1.042466
## 19319        chr7         DLAT  -1.042466
## 58186        chrX        PDHA1  -1.042466
## 31069       chr11         ENO3  -1.049717
## 55895       chr28         ACO2  -1.049717
## 18092        chr6      ATP5F1B  -1.057018
## 16692        chr6         TPI1  -1.057018
## 50105       chr23      ALDH1A1  -1.057018
## 16596        chr6       NDUFA9  -1.057018
## 36397       chr14        HSPA9  -1.057018
## 34419       chr13         TUFM  -1.057018
## 47977       chr21          NNT  -1.064370
## 8965         chr3       COX4I1  -1.064370
## 31165       chr11       ACADVL  -1.064370
## 36458       chr14         MYOT  -1.064370
## 32902       chr12        ACTN3  -1.071773
## 44131       chr18       NDUFS1  -1.071773
## 36548       chr14        VDAC1  -1.071773
## 16228        chr6      NDUFA10  -1.071773
## 26380       chr10          CKM  -1.079228
## 43492       chr18     SLC25A12  -1.079228
## 49924       chr23         FBP2  -1.079228
## 47694       chr21      NDUFA13  -1.079228
## 56767       chr29      ATP5F1C  -1.086735
## 1642         chr1        VDAC2  -1.086735
## 33984       chr13         MDH2  -1.094294
## 39029       chr15        HADHA  -1.094294
## 18049        chr6           CS  -1.094294
## 42150       chr17       SUCLA2  -1.094294
## 38024       chr15       SUCLG1  -1.094294
## 23673        chr8      ATP5F1A  -1.101905
## 1562         chr1        MYOZ1  -1.101905
## 39960       chr16         PDHB  -1.101905
## 54395       chr27        VDAC3  -1.101905
## 4528         chr1         CFL2  -1.101905
## 55445       chr28       MYBPC1  -1.109569
## 31953       chr12       NDUFS3  -1.109569
## 10585        chr3       ATP5ME  -1.109569
## 11968        chr4       CHCHD3  -1.109569
## 14715        chr5       ATP5PB  -1.117287
## 2850         chr1        COX5A  -1.117287
## 36562       chr14 LOC111767815  -1.117287
## 38347       chr15         MDH1  -1.125058
## 37713       chr15        COX5B  -1.125058
## 60627 Contig62106       NDUFA2  -1.125058
## 44459       chr19       NDUFB5  -1.125058
## 54854       chr27      SLC25A4  -1.132884
## 19237        chr7        ACAT1  -1.132884
## 33982       chr13        HSPB1  -1.140764
## 19310        chr7        CRYAB  -1.140764
## 33014       chr12       NDUFS8  -1.140764
## 27080       chr10       NDUFA3  -1.140764
## 25661       chr10       COX6B1  -1.164734
## 37125       chr14        CKMT2  -1.172835
## 29015       chr11       ATP5PD  -1.180993
## 3319         chr1        ANXA2  -1.180993
## 24653        chr9        UQCRB  -1.189207
## 43885       chr18 LOC100070523  -1.189207
## 27303       chr10        TNNT1  -1.197479
## 22797        chr8       ATP2A2  -1.205808
## 60892          MT         ATP8  -1.205808
## 53737       chr26       ATP5PF  -1.205808
## 42086       chr17        KPNA3  -1.214195
## 31377       chr11         MYH1  -1.222640
## 55592       chr28           MB  -1.239708
## 58444        chrX      NDUFB11  -1.265757
## 5858         chr2      ATP5IF1  -1.301342
## 57354       chr30        TNNI1  -1.310393
## 40160       chr16        TNNC1  -1.310393
## 22770        chr8         MYL2  -1.394744
## 40620       chr16         MYL3  -1.404445
##                                                                                                          product
## 11875                                                                                       filamin-C isoform X1
## 16706                                                                                               prohibitin-2
## 3016                                                                              pyruvate kinase PKM isoform M1
## 19319 dihydrolipoyllysine-residue acetyltransferase component of pyruvate dehydrogenase complex%2C mitochondrial
## 58186                         pyruvate dehydrogenase E1 component subunit alpha%2C somatic form%2C mitochondrial
## 31069                                                                                               beta-enolase
## 55895                                                                       aconitate hydratase%2C mitochondrial
## 18092                                                                 ATP synthase subunit beta%2C mitochondrial
## 16692                                                                                  triosephosphate isomerase
## 50105                                                                                    retinal dehydrogenase 1
## 16596                              NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 9%2C mitochondrial
## 36397                                                                         stress-70 protein%2C mitochondrial
## 34419                                                                      elongation factor Tu%2C mitochondrial
## 47977                                                        NAD(P) transhydrogenase%2C mitochondrial isoform X1
## 8965                                                   cytochrome c oxidase subunit 4 isoform 1%2C mitochondrial
## 31165                                very long-chain specific acyl-CoA dehydrogenase%2C mitochondrial isoform X6
## 36458                                                                                        myotilin isoform X1
## 32902                                                                                            alpha-actinin-3
## 44131                                  NADH-ubiquinone oxidoreductase 75 kDa subunit%2C mitochondrial isoform X1
## 36548                                                        voltage-dependent anion-selective channel protein 1
## 16228                  NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 10%2C mitochondrial isoform X1
## 26380                                                                                     creatine kinase M-type
## 43492                                           calcium-binding mitochondrial carrier protein Aralar1 isoform X1
## 49924                                                                    fructose-1%2C6-bisphosphatase isozyme 2
## 47694                                              NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 13
## 56767                                                     ATP synthase subunit gamma%2C mitochondrial isoform X3
## 1642                                                         voltage-dependent anion-selective channel protein 2
## 33984                                                                      malate dehydrogenase%2C mitochondrial
## 39029                                                        trifunctional enzyme subunit alpha%2C mitochondrial
## 18049                                                                          citrate synthase%2C mitochondrial
## 42150                                          succinate--CoA ligase [ADP-forming] subunit beta%2C mitochondrial
## 38024                          succinate--CoA ligase [ADP/GDP-forming] subunit alpha%2C mitochondrial isoform X1
## 23673                                                                ATP synthase subunit alpha%2C mitochondrial
## 1562                                                                                                  myozenin-1
## 39960                                          pyruvate dehydrogenase E1 component subunit beta%2C mitochondrial
## 54395                                             voltage-dependent anion-selective channel protein 3 isoform X1
## 4528                                                                                        cofilin-2 isoform X1
## 55445                                                           myosin-binding protein C%2C slow-type isoform X2
## 31953                          NADH dehydrogenase [ubiquinone] iron-sulfur protein 3%2C mitochondrial isoform X1
## 10585                                                                    ATP synthase subunit e%2C mitochondrial
## 11968                                                                     MICOS complex subunit MIC19 isoform X2
## 14715                                                      ATP synthase F(0) complex subunit B1%2C mitochondrial
## 2850                                                            cytochrome c oxidase subunit 5A%2C mitochondrial
## 36562                                                                          cytochrome b-c1 complex subunit 8
## 38347                                                                        malate dehydrogenase%2C cytoplasmic
## 37713                                                           cytochrome c oxidase subunit 5B%2C mitochondrial
## 60627                                               NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 2
## 44459                               NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 5%2C mitochondrial
## 54854                                                                                      ADP/ATP translocase 1
## 19237                                                              acetyl-CoA acetyltransferase%2C mitochondrial
## 33982                                                                                  heat shock protein beta-1
## 19310                                                                                   alpha-crystallin B chain
## 33014                                     NADH dehydrogenase [ubiquinone] iron-sulfur protein 8%2C mitochondrial
## 27080                                               NADH dehydrogenase [ubiquinone] 1 alpha subcomplex subunit 3
## 25661                                                                           cytochrome c oxidase subunit 6B1
## 37125                                                                    creatine kinase S-type%2C mitochondrial
## 29015                                                         ATP synthase subunit d%2C mitochondrial isoform X1
## 3319                                                                                       annexin A2 isoform X1
## 24653                                                                          cytochrome b-c1 complex subunit 7
## 43885                                                      10 kDa heat shock protein%2C mitochondrial isoform X2
## 27303                                                              troponin T%2C slow skeletal muscle isoform X1
## 22797                                             sarcoplasmic/endoplasmic reticulum calcium ATPase 2 isoform X1
## 60892                                                                                  ATP synthase F0 subunit 8
## 53737                                                            ATP synthase-coupling factor 6%2C mitochondrial
## 42086                                                                        importin subunit alpha-4 isoform X1
## 31377                                                                                        myosin-1 isoform X1
## 55592                                                                                                  myoglobin
## 58444                              NADH dehydrogenase [ubiquinone] 1 beta subcomplex subunit 11%2C mitochondrial
## 5858                                                                ATPase inhibitor%2C mitochondrial isoform X1
## 57354                                                                         troponin I%2C slow skeletal muscle
## 40160                                                            troponin C%2C slow skeletal and cardiac muscles
## 22770                                      myosin regulatory light chain 2%2C ventricular/cardiac muscle isoform
## 40620                                                                                       myosin light chain 3
```

```r
# Number
sum(idx)
```

```
## [1] 73
```

### Write results to file


```r
write.table(sig.PT, 
    file=paste(getwd(), "Summary_Results_RER_Thoroughbred.txt", sep="/"), 
    col.names=TRUE, row.names=FALSE, quote=FALSE, sep="\t")
```

Save R Data 


```r
save(anim, sig.PT, pval.cut, file=paste(getwd(), "Summary_Results_RER_Thoroughbred.Rdata", sep="/"))
```

### Run R Script


```r
htmlRunR
Summary_Results.R nodes=1,cpus-per-task=1,time=02:00:00,mem=5G \
+RER Thoroughbred Differentially Expressed Proteins
```

