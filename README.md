# Survival Analysis: Stomach Cancer
This project involves conducting a survival analysis of stomach cancer, using gene expression data from the TCGA database. Various survival analysis techniques were employed, including:

1) Kaplan-Meier (KM) Estimation
2) Cox Proportional Hazards (CPH)
3) Random Survival Forest (RSF): identifying top predictors of stomach cancer, out of 20,000 genes.

All code was done in R and is available in the file, "tcga_survival_analysis.R"


### Data
Data were derived from the the stomach cancer cohort in the TCGA database. The following three data sets were used:
1) IlluminaHiSeq pancan normalized: gene expression data
2) Phenotypes
3) Curated survival data

After joining all three datasets, there were a total of 450 samples from 418 individual patients. Samples includes both primary tumor and normal solid tissue samples. 

The **gene expression data** consisted of 20,5311 genes. Gene expression values quantified the abundance of RNA transcribed from a given genes.

The **phenotype data** consisted of tissue sample type, cancer histology, cancer stage, patient characteristics, etc.

The **survival data** consisted of four metrics measuring patient survival and recurrence, including overall survival (OS), disease specific survival (DSS), disease free interval (DFI), and progression free interval (PFI), each measured in days.

Detailed exploratory analysis is included in the code file, "tcga_prelim_analysis.R"

Data are publicly available through the UCSC Xena Browser.
https://xenabrowser.net/datapages/?cohort=TCGA%20Stomach%20Cancer%20(STAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443 

### Survival Analysis

#### Kaplan-Meier Estimation
Prior literature shows that high expression of the PRTG gene is associated with poor gastric cancer prognosis (Xiang, et al. 2021). The figure below shows the KM curve, stratified by PRTG gene expression levels. Gene expression levels were categorized into two groups: low gene expression and high gene expression, using the median gene expression level as a threshold.

Code is included in "tcga_prelim_analysis.R"



##### Disease Specific Survival of Stomach Cancer by PRTG Gene Expression
<img width="1279" height="814" alt="prtg" src="https://github.com/user-attachments/assets/74ae393e-abf4-4c5a-8acd-518b11e42291" />

There appeared to be a statistically significant association between PRTG gene expression and disease specific survival, with high PRTG gene expression associated with poorer survival. Similar curves were created to analyze the relationship between different genes and stomach cancer prognosis. 

#### Cox Proportional Hazards: PRTG Gene Expression
<img width="633" height="318" alt="prtg_cox" src="https://github.com/user-attachments/assets/815fcc5d-c649-4340-8d6b-abe13a6cbc99" />

There appears to be a statistically significant, negative correlation between low PRTG gene expression and disease specific survival of stomach cancer.

Code is included in "tcga_prelim_analysis.R"


#### Cox Proportional Hazards: All Genes, Multiple Testing Correction
Additionally, we performed Cox PH iteratively across each individual gene, to understand which genes were most signficantly associated with disease specific survival. We used the Benjamini-Hochberg multiple testing correction method to obtain adjusted p-values.

##### Top five genes, along with hazard ratios and p-values.
<img width="342" height="177" alt="top_cox" src="https://github.com/user-attachments/assets/3fa42d72-ba83-4a35-a24a-0c82627ca5f7" />

Code is included in "tcga_survival_prediction.Rmd"


#### Random Survival Forest
We next aimed to understand how gene expression could predict disease specific survival of stomach cancer, using the Random Survival Forest technique. After training and testing the model and tuning its hyperparameters (number of trees & number of variables randomly selected at each node split). 

The top five genes obtained from the iterative CPH analysis were used to visualize a decision tree in the RSF model.
<img width="508" height="420" alt="decisiontree" src="https://github.com/user-attachments/assets/6758bc6c-6c5b-4165-91f7-6bf6729509ef" />

Code is included in "tcga_survival_prediction.Rmd" Used package randomforestSRC for RSF analyses

### Sources
Xiang T, Yuan C, Guo X, et al. The novel ZEB1-upregulated protein PRTG induced by Helicobacter pylori infection promotes gastric carcinogenesis through the cGMP/PKG signaling pathway. Cell Death Dis. 2021;12(2):150. doi:10.1038/s41419-021-03440-1

