# DrGA: driver gene analysis in an automatic manner
#### I. Introduction
---
DrGA is a novel R package that has been developed based on the idea of our most recent driver gene analysis scheme [here](https://www.nature.com/articles/s41598-020-77318-1). Its aim is to wholy automate the analysis process of driver genes at a low investment of time for this process by merging state-of-the-art statistical tools into one.

#### II. Understanding the tool and Data Structure
---
The following are parameters included in DrGA and their role:
- organism: organism name. Organism names are constructed by concatenating the first letter of the name and the family name. Example: human - hsapiens, mouse - mmusculus. Default is organism = "hsapiens".

- sources: possible biological mechanisms allowed (e.g., Gene Ontology - GO:BP, GO:MF,GO:CC; KEGG; REAC; TF; MIRNA; CORUM; HP; HPA; WP;… Please see the g:GOSt web tool for the comprehensive list and details on incorporated data sources). Default is sources = c("GO:BP", "KEGG").

- methodCC: Correlation method type. Allowed values are spearman (default), pearson, kendall.

- exp: a data frame or matrix. exp has its rows are samples and its columns are genes. It is input data to serve to run the second and third modules.

- clinicalEXP: a data frame or matrix. It includes its rows are samples, and its columns are clinical features of your choice. Notice that the clinical data must have two columns as overall survival time (continuous variable) and overall survival status (binary variable; usually coded 1 as death and 0 as live) of all the subjects.

- timeEXP: a vector of overall survival time. It is a column vector of clinicalEXP.

- statusEXP: a vector of overall survival time. It is a column vector of clinicalEXP.

- datMODULE4: a data frame or matrix. datMODULE4 has its rows are samples and its columns are genes. It is input data to serve to run the forth module.

- cliMODULE4: a data frame or matrix. It includes its rows are samples, and its columns are clinical features of your choice. Notice that the clinical data must have two columns as overall survival time (continuous variable) and overall survival status (binary variable; usually coded 1 as death and 0 as live) of all the subjects.

- timeMODULE4: a vector of overall survival time. It is a column vector of clinicalMODULE4.

- statusMODULE4: a vector of overall survival time. It is a column vector of clinicalMODULE4.

- minClusterSize: Minimum cluster size. minClusterSize = 10 is default.

- NetworkType: network type. Allowed values are (unique abbreviations of) "unsigned", "signed", "signed hybrid". Default value is signed.

- hm_row_names: logical. If hm_row_names = TRUE (default value), gene names appear in rows of the heatmap.  If due to the large number of driver genes leading to impossibly showing gene names in rows of the heatmap, users can turn them off by hm_row_names = FALSE.

Please download datasets [data_n_code](https://github.com/huynguyen250896/DrGA/tree/master/data_n_code) as examples to well grasp DrGA's requirement on data structure.

#### III. Pipeline
---
![Figure](https://imgur.com/mquJy2O.png)
**Figure:** Pipeline of the package DrGA.

#### IV. Implementation
---
Use the following command to install directly from GitHub;
```sh
devtools::install_github("huynguyen250896/DrGA")
```
Call the library;
```sh
x = c("DrGA", "dplyr", "survival", "tibble", "tidyr", "ComplexHeatmap", 
     'cluster', 'mclust', 'clValid', 'Biobase', 'annotate', 'GO.db', 
     'mygene', "dynamicTreeCut", "flashClust", "Hmisc", "WGCNA","purrr",
     "gprofiler2", "table1", "compareGroups")
lapply(x, require, character.only = TRUE)
```
running example:
```sh
DriverGeneAnalysis(exp = exp, clinicalEXP = clinicalEXP, timeEXP = clinicalEXP$time, statusEXP = clinicalEXP$status, 
                   datMODULE4 = cna,  cliMODULE4 = clinicalCNA, timeMODULE4 = clinicalCNA$time, statusMODULE4 = clinicalCNA$status)
```

Feel free to contact [Quang-Huy Nguyen](https://github.com/huynguyen250896) <huynguyen96.dnu AT gmail DOT com> for any questions about the code and results.
