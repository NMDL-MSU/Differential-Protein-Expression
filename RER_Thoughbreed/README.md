# Differential-Protein-Expression
## Project RER Thoroughbred

**Underlying Cause for Recurrent Exertional Rhabdomyolysis (RER)**

*Hypothesis:*  
Fit Thoroughbred fillies susceptible to RER in race training have differential expression of biological pathways and genes involved in intramuscular Ca2+ regulation compared to horses that are not susceptible to RER. 

### Filtering and Annotating Proteins
Normalized protein expressions were obtained from all ten horses (five RER and five controls) using Scaffold. Proteins not expressed in all of the samples including any technical replicate were discarded from further analysis. The script used for this step can be viewed [here](https://htmlpreview.github.io/?https://github.com/NMDL-MSU/Differential-Protein-Expression/blob/master/RER_Thoughbreed/Annotate_Proteins/AnnotateProteins.html).

### Differential Protein Expression Analysis
Differentially expressed protein between horses with RER and controls were identified in Scaffold Q+ using a permutation test and FDR 0.05.

#### &nbsp;&nbsp;&nbsp;&nbsp;DE Analysis:
The scripts used to generate the summary tables of the permutation test can be viewed in detail [here](https://htmlpreview.github.io/?https://github.com/NMDL-MSU/Differential-Protein-Expression/blob/master/RER_Thoughbreed/Permutation_Test/Summary_Results.html).

#### &nbsp;&nbsp;&nbsp;&nbsp;DE Results:
The DE results for all comparisons can be downloaded [here](https://github.com/NMDL-MSU/Differential-Protein-Expression/blob/master/RER_Thoughbreed/Permutation_Test/Summary_Results_RER_Thoroughbred.xlsx?raw=true).
