
### Description:  
#' Performe protein set enrichment analysis for RER Thoroughbred differential protein expression. 
#' Look at enrichment of differentially expressed proteins in horses presenting recurrent exertional  
#' rhabdomyolysis (RER) compared to controls.  
#'   
#' ***  
#' **Code:**  
#' Parent Directory:  
#' 
#' >&nbsp;&nbsp;&nbsp;&nbsp;/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred    
#'   
#' Directory/File:  
#'  
#' &nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/Enrichment/RER_Thoroughbred/RER_Thoroughbred.R  
#'  
#' **Input files:**  
#' Directory/File:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results_RER_Thoroughbred.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/EquCab3/EquCab3_Annotation.Rdata  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/RNA_Seq/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata  
#'   
#' **Output files:**  
#'   
#' Directory:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Proteomics/Enrichment/RER_Thoroughbred/  
#'   
#' Files:  
#'   
#' >&nbsp;&nbsp;&nbsp;&nbsp;/Merged/*.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/UpRegulated/*.txt  
#' >&nbsp;&nbsp;&nbsp;&nbsp;/DownRegulated/*.txt  
#'   
#' Render R Script  
#'   
#' > &nbsp;&nbsp;&nbsp;&nbsp;RER_Thoroughbred.qsub  
#'  
#' ***  

#' ### R Environment    

#' Load required libraries
library(edgeR)
library(limma)
library(DOSE)
library(GO.db)
library(GSEABase)
library(clusterProfiler)
library(org.Hs.eg.db)
library(GenomicFeatures)

#' Clear Environment
rm(list=ls())

#' Session Information
sessionInfo()


#' ### Load Data for Pathway Analysis

#' Load Proteomic Data
dir <- "/mnt/research/NMDL/KER_Glycogen_and_RER_Thoroughbred"
dirP <- "/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/RER_Thoroughbred_Melissa_20180817_Protein_Annotation.Rdata"
load(paste(dir, dirP, sep=""))

#' Load DE Results
dirR <- "/Proteomics/RER_Thoroughbred_Melissa_20180817/AnnotateProteins/PermutationTest/Summary_Results/Summary_Results_RER_Thoroughbred.Rdata"
load(paste(dir, dirR, sep=""))

#' Load EquCab3 Annotation 
dirA <- "/Proteomics/EquCab3/EquCab3_Annotation.Rdata"
load(paste(dir, dirA, sep=""))

#' Load enrichment function
Fun <- "/RNA_Seq/Enrichment/Enrichment_KER_Glycogen_Depl/enrichment_function.Rdata"
load(paste(dir, Fun, sep=""))


#' ### Animal Groups
# Animals with proteomics data
idx <- c("12613", "12910", "12915", "12916", "12918", "12401", "12402", "12403", "12620", "12913")
id <- anim$ID
anim <- anim[,-2]
rownames(anim) <- id
anim <- anim[idx,]
group <- anim$Dx
table(group)


#' ### Differential Expression Analysis Results

#' Significant gene names
sigG <- as.character(sig.PT$gene)
names(sigG) <- sig.PT$Genbank

#' Seperate genes without a name
# Number of unknown gene transcripts
na.sigG <- sigG[is.na(sigG)]
length(na.sigG)

# Number of known genes
sigG <- sigG[!is.na(sigG)]
length(sigG)


#' ### Background for Enrichment Analysis

#' Background Gene List
annot <- do.call(rbind, lapply(NormInt$AccessionNumber, function(x) ProtEquCab3[as.character(ProtEquCab3$Genbank) %in% x,]))

#' Seperate genes without a name
na.Bkg <- annot[is.na(annot$gene),]
nrow(na.Bkg)

Bkg <- annot[!is.na(annot$gene),]
nrow(Bkg)

#' Obtain human EntrezIDs for gene enrichment analysis
# DE gene list
sigG.entrez <- bitr(unique(sigG), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
nrow(sigG.entrez)

# Background gene list
bkg.entrez <- bitr(as.character(unique(Bkg$gene)), fromType="SYMBOL", 
    toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")
nrow(bkg.entrez)


#' ### Enrichment Analysis  

#' Global gene enrichment analysis: Merge all DE genes
enrich.rst <- enrich(lst=sigG.entrez, bkg=bkg.entrez)
unlist(lapply(enrich.rst, nrow))

#' Summary of significant GO terms
lapply(enrich.rst, function(x) head(x[,2]))

#' Save enrichment analysis results to file
# Enrichemnt analysis results
data_merged <- lapply(enrich.rst, function(x) data.frame(x))
data_merged <- data_merged[unlist(lapply(data_merged, nrow)) > 0]

# Save results to file
system("mkdir Merged")
lapply(names(data_merged), function(x) 
    write.table(data_merged[[x]], file=paste(getwd(), "/Merged/", x, "_merged.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))


#' ### Enrichment of Up-regulated Genes  

#' Gene set enrichment for upregulated genes
up <- sig.PT[sig.PT$FoldChange > 0,]
up <- lapply(as.character(up$gene), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
up <- do.call(rbind, up[unlist(lapply(up, nrow)) > 0])
nrow(up)


#' Gene enrichment analysis for upregulated genes:
enrich.upreg <- enrich(lst=up, bkg=bkg.entrez)
unlist(lapply(enrich.upreg, nrow))

#' Summary of significant GO terms
lapply(enrich.upreg, function(x) head(x[,2]))

#' Save results
data_upreg <- lapply(enrich.upreg, function(x)  
        data.frame(x))
data_upreg <- data_upreg[unlist(lapply(data_upreg, nrow)) > 0]

# Save results to file
system("mkdir UpRegulated")
lapply(names(data_upreg), function(x) 
    write.table(data_upreg[[x]], file=paste(getwd(), "/UpRegulated/", x, "_upregulated.txt", sep=""),
    col.names=TRUE, row.names=TRUE, quote=FALSE, sep="\t"))



#' ### Enrichment of Down-regulated Genes  

#' Gene set enrichment for down-regulated genes
down <- sig.PT[sig.PT$FoldChange < 0,]
down <- lapply(as.character(down$gene), 
    function(x) sigG.entrez[sigG.entrez$SYMBOL %in% x,])
down <- do.call(rbind, down[unlist(lapply(down, nrow)) > 0])
nrow(down)


#' Gene enrichment analysis for down-regulated genes:
enrich.downreg <- enrich(lst=down, bkg=bkg.entrez)
unlist(lapply(enrich.downreg, nrow))

#' Summary of significant GO terms
lapply(enrich.downreg, function(x) head(x[,2]))

#' Save results
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



#' ### Run R Script
#+ eval = FALSE
htmlRunR
RER_Thoroughbred.R nodes=1,cpus-per-task=1,time=03:00:00,mem=50G \
+Protein Set Enrichment RER Thoroughbred Project

