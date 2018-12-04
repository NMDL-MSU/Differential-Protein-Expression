# Differential-Protein-Expression
## Project KER Glycogen

**Critical pathways limiting glycogen repletion and performance in the horse defined by
proteomic and transcriptomic analyses**

*Objective:*  
To determine the time course of differential expression specific genes and biological pathways during
glycogen depletion and repletion in fit thoroughbred horses fed a high starch (Sweet Feed) compared to an
isocaloric low starch, high fat (Fat) diet. 

### Study design
Five Thoroughbred geldings were given two diets (low and high starch) with appropriate washout periods between diet. 
Muscle biopsies were taken from the gluteal muscle at four-time points: pre-depletion of glycogen, depletion, repletion at 24h and repletion at 72h. 

### Filtering and Annotating Proteins
Normalized protein expressions were obtained from all horses using Scaffold. Proteins not expressed in all of the samples including any technical replicate were discarded from further analysis. The script used for this step can be viewed [here](https://htmlpreview.github.io/?https://github.com/NMDL-MSU/Differential-Protein-Expression/blob/master/KER_Glycogen/Permutation_Test/Summary_Results.html).

### Differential Protein Expression Analysis
Determined differentially expressed protein was identified for the following comparisons in Scaffold Q+ using a permutation test and FDR 0.05:  
* Between diets per time point  
* Between diets over time  
* Between time points per diet  
* Between time points per diet over time  

#### &nbsp;&nbsp;&nbsp;&nbsp;DE Analysis:
The scripts used to generate the summary tables of the permutation test can be viewed in detail [here](https://htmlpreview.github.io/?https://github.com/NMDL-MSU/Differential-Protein-Expression/blob/master/KER_Glycogen/Permutation_Test/Summary_Results.html).

#### &nbsp;&nbsp;&nbsp;&nbsp;DE Results:
The DE results for all comparisons can be downloaded [here](https://github.com/NMDL-MSU/Differential-Protein-Expression/blob/master/KER_Glycogen/Permutation_Test/Results_DE_Proteins.xlsx?raw=true).

