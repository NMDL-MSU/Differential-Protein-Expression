---
title: Protein Set Enrichment RER Thoroughbred Project
author: Deborah Velez-Irizarry
date: Tue Dec 11 10:24:46 EST 2018
output:
  prettydoc::html_pretty:
    theme: tactile
    highlight: github
    toc: true
---


```r
### Description:  
```

Performe protein set enrichment analysis for RER Thoroughbred differential protein expression. 
Look at enrichment of differentially expressed proteins in horses presenting recurrent exertional  
rhabdomyolysis (RER) compared to controls.  
  
***  
**Code:**  
Parent Directory:  

>&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred    
  
Directory/File:  
 
&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/Enrichment/RER_Thoroughbred/RER_Thoroughbred.R  
 
**Input files:**  
Directory/File:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results_RER_Thoroughbred.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/EquCab3/EquCab3_Annotation.Rdata  
>&nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata  
  
**Output files:**  
  
Directory:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/Enrichment/RER_Thoroughbred/  
  
Files:  
  
>&nbsp;&nbsp;&nbsp;&nbsp;/Merged/*.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;/UpRegulated/*.txt  
>&nbsp;&nbsp;&nbsp;&nbsp;/DownRegulated/*.txt  
  
Render R Script  
  
> &nbsp;&nbsp;&nbsp;&nbsp;RER_Thoroughbred.qsub  
 
***  
### R Environment    
Load required libraries


```r
library(edgeR)
```

```
## Loading required package: limma
```

```r
library(limma)
library(DOSE)
```

```
## 
```

```
## DOSE v3.8.0  For help: https://guangchuangyu.github.io/DOSE
## 
## If you use DOSE in published research, please cite:
## Guangchuang Yu, Li-Gen Wang, Guang-Rong Yan, Qing-Yu He. DOSE: an R/Bioconductor package for Disease Ontology Semantic and Enrichment analysis. Bioinformatics 2015, 31(4):608-609
```

```r
library(GO.db)
```

```
## Loading required package: AnnotationDbi
```

```
## Loading required package: stats4
```

```
## Loading required package: BiocGenerics
```

```
## Loading required package: parallel
```

```
## 
## Attaching package: 'BiocGenerics'
```

```
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
```

```
## The following object is masked from 'package:limma':
## 
##     plotMA
```

```
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
```

```
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, basename, cbind,
##     colMeans, colnames, colSums, dirname, do.call, duplicated,
##     eval, evalq, Filter, Find, get, grep, grepl, intersect,
##     is.unsorted, lapply, lengths, Map, mapply, match, mget, order,
##     paste, pmax, pmax.int, pmin, pmin.int, Position, rank, rbind,
##     Reduce, rowMeans, rownames, rowSums, sapply, setdiff, sort,
##     table, tapply, union, unique, unsplit, which, which.max,
##     which.min
```

```
## Loading required package: Biobase
```

```
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
```

```
## Loading required package: IRanges
```

```
## Loading required package: S4Vectors
```

```
## 
## Attaching package: 'S4Vectors'
```

```
## The following object is masked from 'package:base':
## 
##     expand.grid
```

```r
library(GSEABase)
```

```
## Loading required package: annotate
```

```
## Loading required package: XML
```

```
## Loading required package: graph
```

```
## 
## Attaching package: 'graph'
```

```
## The following object is masked from 'package:XML':
## 
##     addNode
```

```r
library(clusterProfiler)
```

```
## clusterProfiler v3.10.0  For help: https://guangchuangyu.github.io/software/clusterProfiler
## 
## If you use clusterProfiler in published research, please cite:
## Guangchuang Yu, Li-Gen Wang, Yanyan Han, Qing-Yu He. clusterProfiler: an R package for comparing biological themes among gene clusters. OMICS: A Journal of Integrative Biology. 2012, 16(5):284-287.
```

```r
library(org.Hs.eg.db)
```

```
## 
```

```r
library(GenomicFeatures)
```

```
## Loading required package: GenomeInfoDb
```

```
## Loading required package: GenomicRanges
```

Clear Environment


```r
rm(list=ls())
```

Session Information


```r
sessionInfo()
```

```
## R version 3.5.1 (2018-07-02)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: CentOS Linux 7 (Core)
## 
## Matrix products: default
## BLAS/LAPACK: /opt/software/OpenBLAS/0.3.1-GCC-7.3.0-2.30/lib/libopenblas_sandybridgep-r0.3.1.so
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
## [1] parallel  stats4    stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] GenomicFeatures_1.34.1 GenomicRanges_1.34.0   GenomeInfoDb_1.18.1   
##  [4] org.Hs.eg.db_3.7.0     clusterProfiler_3.10.0 GSEABase_1.44.0       
##  [7] graph_1.60.0           annotate_1.60.0        XML_3.98-1.16         
## [10] GO.db_3.7.0            AnnotationDbi_1.44.0   IRanges_2.16.0        
## [13] S4Vectors_0.20.1       Biobase_2.42.0         BiocGenerics_0.28.0   
## [16] DOSE_3.8.0             edgeR_3.24.1           limma_3.38.3          
## [19] knitr_1.20            
## 
## loaded via a namespace (and not attached):
##  [1] fgsea_1.8.0                 colorspace_1.3-2           
##  [3] ggridges_0.5.1              qvalue_2.14.0              
##  [5] XVector_0.22.0              farver_1.1.0               
##  [7] urltools_1.7.0              ggrepel_0.8.0              
##  [9] bit64_0.9-7                 xml2_1.2.0                 
## [11] splines_3.5.1               GOSemSim_2.8.0             
## [13] jsonlite_1.5                Rsamtools_1.34.0           
## [15] ggforce_0.1.3               compiler_3.5.1             
## [17] httr_1.3.1                  rvcheck_0.1.3              
## [19] assertthat_0.2.0            Matrix_1.2-14              
## [21] lazyeval_0.2.1              tweenr_1.0.0               
## [23] prettyunits_1.0.2           tools_3.5.1                
## [25] bindrcpp_0.2.2              igraph_1.2.1               
## [27] gtable_0.2.0                glue_1.3.0                 
## [29] GenomeInfoDbData_1.2.0      reshape2_1.4.3             
## [31] DO.db_2.9                   dplyr_0.7.8                
## [33] fastmatch_1.1-0             Rcpp_1.0.0                 
## [35] enrichplot_1.2.0            Biostrings_2.50.1          
## [37] rtracklayer_1.42.1          ggraph_1.0.2               
## [39] stringr_1.3.1               europepmc_0.3              
## [41] MASS_7.3-51.1               zlibbioc_1.28.0            
## [43] scales_1.0.0                hms_0.4.2                  
## [45] SummarizedExperiment_1.12.0 RColorBrewer_1.1-2         
## [47] memoise_1.1.0               gridExtra_2.3              
## [49] ggplot2_3.1.0.9000          UpSetR_1.3.3               
## [51] biomaRt_2.38.0              triebeard_0.3.0            
## [53] stringi_1.2.3               RSQLite_2.1.1              
## [55] BiocParallel_1.16.2         rlang_0.3.0.1              
## [57] pkgconfig_2.0.2             bitops_1.0-6               
## [59] matrixStats_0.54.0          evaluate_0.12              
## [61] lattice_0.20-38             purrr_0.2.5                
## [63] bindr_0.1.1                 GenomicAlignments_1.18.0   
## [65] cowplot_0.9.2               bit_1.1-14                 
## [67] tidyselect_0.2.5            plyr_1.8.4                 
## [69] magrittr_1.5                R6_2.3.0                   
## [71] DelayedArray_0.8.0          DBI_1.0.0                  
## [73] pillar_1.2.3                units_0.6-2                
## [75] RCurl_1.95-4.11             tibble_1.4.2               
## [77] crayon_1.3.4                viridis_0.5.1              
## [79] progress_1.2.0              locfit_1.5-9.1             
## [81] grid_3.5.1                  data.table_1.11.8          
## [83] blob_1.1.1                  digest_0.6.18              
## [85] xtable_1.8-2                tidyr_0.8.1                
## [87] gridGraphics_0.3-0          munsell_0.5.0              
## [89] viridisLite_0.3.0           ggplotify_0.0.3
```

### Load Data for Pathway Analysis
Load Proteomic Data


```r
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred"
dirP <- "/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata"
load(paste(dir, dirP, sep=""))
```

Load DE Results


```r
dirR <- "/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results_RER_Thoroughbred.Rdata"
load(paste(dir, dirR, sep=""))
```

Load EquCab3 Annotation 


```r
dirA <- "/Proteomics/EquCab3/EquCab3_Annotation.Rdata"
load(paste(dir, dirA, sep=""))
```

Load enrichment function


```r
Fun <- "/RNA_Seq/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata"
load(paste(dir, Fun, sep=""))
```

### Animal Groups


```r
# Animals with proteomics data
idx <- c("12613", "12910", "12915", "12916", "12918", "12401", "12402", "12403", "12620", "12913")
id <- anim$ID
anim <- anim[,-2]
rownames(anim) <- id
anim <- anim[idx,]
group <- anim$Dx
table(group)
```

```
## group
## Control     RER 
##       5       5
```

### Differential Expression Analysis Results
Significant gene names


```r
sigG <- as.character(sig.PT$gene)
names(sigG) <- sig.PT$Genbank
```

Seperate genes without a name


```r
# Number of unknown gene transcripts
na.sigG <- sigG[is.na(sigG)]
length(na.sigG)
```

```
## [1] 0
```

```r
# Number of known genes
sigG <- sigG[!is.na(sigG)]
length(sigG)
```

```
## [1] 125
```

### Background for Enrichment Analysis
Background Gene List


```r
annot <- do.call(rbind, lapply(NormInt$AccessionNumber, function(x) ProtEquCab3[as.character(ProtEquCab3$Genbank) %in% x,]))
```

Seperate genes without a name


```r
na.Bkg <- annot[is.na(annot$gene),]
nrow(na.Bkg)
```

```
## [1] 0
```

```r
Bkg <- annot[!is.na(annot$gene),]
nrow(Bkg)
```

```
## [1] 375
```

Obtain human EntrezIDs for gene enrichment analysis


```r
# DE gene list
sigG.entrez <- bitr(unique(sigG), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Warning in bitr(unique(sigG), fromType = "SYMBOL", toType =
## c("ENTREZID"), : 8.8% of input gene IDs are fail to map...
```

```r
nrow(sigG.entrez)
```

```
## [1] 114
```

```r
# Background gene list
bkg.entrez <- bitr(as.character(unique(Bkg$gene)), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
```

```
## 'select()' returned 1:1 mapping between keys and columns
```

```
## Warning in bitr(as.character(unique(Bkg$gene)), fromType = "SYMBOL", toType
## = c("ENTREZID"), : 8.27% of input gene IDs are fail to map...
```

```r
nrow(bkg.entrez)
```

```
## [1] 344
```

### Enrichment Analysis  
Global gene enrichment analysis: Merge all DE genes


```r
enrich.rst <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##   114     0     0     8     0
```

Summary of significant GO terms


```r
lapply(enrich.rst, function(x) head(x[,2]))
```

```
## $Genes
## [1] "759"  "3043" "335"  "760"  "1356" "7001"
## 
## $GO.BP
## character(0)
## 
## $GO.CC
## character(0)
## 
## $GO.MF
## [1] "transporter activity"                                         
## [2] "transmembrane transporter activity"                           
## [3] "ion transmembrane transporter activity"                       
## [4] "inorganic molecular entity transmembrane transporter activity"
## [5] "inorganic cation transmembrane transporter activity"          
## [6] "cation transmembrane transporter activity"                    
## 
## $Kegg
## character(0)
```

Save enrichment analysis results to file


```r
# Enrichemnt analysis results
data_merged <- lapply(enrich.rst, function(x) data.frame(x))
data_merged <- data_merged[unlist(lapply(data_merged, nrow)) > 0]

# Save results to file
system("mkdir Merged")
lapply(names(data_merged), function(x) 
    write.table(data_merged[[x]], file=paste(getwd(), "/Merged/", x, "_merged.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
```

### Enrichment of Up-regulated Genes  
Gene set enrichment for upregulated genes


```r
up <- sig.PT[sig.PT$FoldChange > 0,]
up <- lapply(as.character(up$gene), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
up <- do.call(rbind, up[unlist(lapply(up, nrow)) > 0])
nrow(up)
```

```
## [1] 43
```

Gene enrichment analysis for upregulated genes:


```r
enrich.upreg <- enrich(lst=up, bkg=bkg.entrez)
unlist(lapply(enrich.upreg, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##    43    58    27     0     0
```

Summary of significant GO terms


```r
lapply(enrich.upreg, function(x) head(x[,2]))
```

```
## $Genes
## [1] "759"  "3043" "335"  "760"  "1356" "7001"
## 
## $GO.BP
## [1] "cellular response to stimulus"          
## [2] "maintenance of location"                
## [3] "cellular protein modification process"  
## [4] "protein modification process"           
## [5] "macromolecule modification"             
## [6] "post-translational protein modification"
## 
## $GO.CC
## [1] "blood microparticle"         "endoplasmic reticulum lumen"
## [3] "extracellular region"        "extracellular region part"  
## [5] "extracellular organelle"     "extracellular vesicle"      
## 
## $GO.MF
## character(0)
## 
## $Kegg
## character(0)
```

Save results


```r
data_upreg <- lapply(enrich.upreg, function(x)  
        data.frame(x))
data_upreg <- data_upreg[unlist(lapply(data_upreg, nrow)) > 0]

# Save results to file
system("mkdir UpRegulated")
lapply(names(data_upreg), function(x) 
    write.table(data_upreg[[x]], file=paste(getwd(), "/UpRegulated/", x, "_upregulated.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
```

### Enrichment of Down-regulated Genes  
Gene set enrichment for down-regulated genes


```r
down <- sig.PT[sig.PT$FoldChange < 0,]
down <- lapply(as.character(down$gene), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
down <- do.call(rbind, down[unlist(lapply(down, nrow)) > 0])
nrow(down)
```

```
## [1] 71
```

Gene enrichment analysis for down-regulated genes:


```r
enrich.downreg <- enrich(lst=down, bkg=bkg.entrez)
unlist(lapply(enrich.downreg, nrow))
```

```
## Genes GO.BP GO.CC GO.MF  Kegg 
##    71    30    19     8     8
```

Summary of significant GO terms


```r
lapply(enrich.downreg, function(x) head(x[,2]))
```

```
## $Genes
## [1] "2318"  "11331" "5315"  "1737"  "5160"  "2027" 
## 
## $GO.BP
## [1] "generation of precursor metabolites and energy"
## [2] "organophosphate metabolic process"             
## [3] "drug metabolic process"                        
## [4] "nucleoside phosphate metabolic process"        
## [5] "nucleotide metabolic process"                  
## [6] "aerobic respiration"                           
## 
## $GO.CC
## [1] "mitochondrial part"            "mitochondrion"                
## [3] "mitochondrial protein complex" "mitochondrial inner membrane" 
## [5] "organelle inner membrane"      "mitochondrial membrane"       
## 
## $GO.MF
## [1] "ion transmembrane transporter activity"                        
## [2] "inorganic molecular entity transmembrane transporter activity" 
## [3] "proton transmembrane transporter activity"                     
## [4] "transmembrane transporter activity"                            
## [5] "monovalent inorganic cation transmembrane transporter activity"
## [6] "transporter activity"                                          
## 
## $Kegg
## [1] "Huntington disease"        "Parkinson disease"        
## [3] "Metabolic pathways"        "Oxidative phosphorylation"
## [5] "Thermogenesis"             "Alzheimer disease"
```

Save results


```r
data_downreg <- lapply(enrich.downreg, function(x)  
        data.frame(x))
names(data_downreg) <- names(enrich.downreg)
data_downreg <- data_downreg[unlist(lapply(data_downreg, nrow)) > 0]

# Add gene names to Kegg results (currently has the geneIDs)
kegg.gene <- sapply(data_downreg$Kegg$geneID, function(x) strsplit(x, "/")[[1]])
names(kegg.gene) <- NULL
data_downreg$Kegg$geneID <- sapply(kegg.gene, function(x) 
    paste(data_downreg$Gene$SYMBOL[data_downreg$Gene$ENTREZID %in% x], collapse="/"))

# Save results to file
system("mkdir DownRegulated")
lapply(names(data_downreg), function(x) 
    write.table(data_downreg[[x]], file=paste(getwd(), "/DownRegulated/", x, "_downregulated.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))
```

```
## [[1]]
## NULL
## 
## [[2]]
## NULL
## 
## [[3]]
## NULL
## 
## [[4]]
## NULL
## 
## [[5]]
## NULL
```

### Run R Script


```r
htmlRunR
RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Protein Set Enrichment RER Thoroughbred Project
```

