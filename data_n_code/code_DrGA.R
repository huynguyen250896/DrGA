#library
devtools::install_github("huynguyen250896/DrGA") #NOTE: Download all dependencies of the tool!
libraries=c("DrGA", "dplyr", "survival", "tibble", "tidyr", "ComplexHeatmap", 
     'cluster', 'mclust', 'clValid', 'Biobase', 'annotate', 'GO.db', 
     'mygene', "dynamicTreeCut", "flashClust", "Hmisc", "WGCNA","purrr",
     "gprofiler2", "table1", "compareGroups")
lapply(libraries, require, character.only = TRUE)

#-------------------------------------o0o-------------------------------------#
# PREPROCESSING PROCEDURES
#-------------------------------------o0o-------------------------------------#

#load raw data
exp = read.table('data_mRNA_median_Zscores.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = NULL)
cna = read.table('data_CNA.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = 1)
clinical = read.table('data_clinical_patient.txt', sep = '\t', check.names = FALSE, header = TRUE, row.names = 1,fill=TRUE)

#identified driver genes
driver=c("MAP2K4", "ARID1A", "PIK3CA", "TBX3", "MAP3K1", "TP53", "AKT1", "GATA3", "CDH1", "RB1", "CDKN1B", "NCOR1", "CDKN2A", "ERBB2", "KRAS", "BRCA2", "BAP1", "PTEN", "CBFB", "KMT2C", "RUNX1", "NF1", "PIK3R1", "ERBB3", "FOXO3", "SMAD4", "GPS2", "AGTR2", "ZFP36L1", "MEN1","SF3B1")
length(driver) #31 driver genes

#only keep the 31 driver genes in exp and cna
exp=exp %>%
  dplyr::filter(.$Hugo_Symbol %in% driver) %>%
  tibble::column_to_rownames('Hugo_Symbol') %>%
  dplyr::select(-Entrez_Gene_Id)

cna=cna[driver, ]

#check dimension
dim(exp) # 31 1904
dim(cna) # 31 2173
dim(clinical) # 2509   21

#match patients sharing between exp versus clinical, and cna versus clinical
#exp and cna are two matrices whose rows are samples and columns are genes
exp = exp[,intersect(colnames(exp), rownames(clinical))] %>% t()
clinicalEXP = clinical[intersect(rownames(exp), rownames(clinical)), ]

cna = cna[,intersect(colnames(cna), rownames(clinical))] %>% t()
clinicalCNA = clinical[intersect(rownames(cna), rownames(clinical)), ]

#preprocess clinicalEXP
clinicalEXP = clinicalEXP %>%
  tibble::rownames_to_column('sample') %>%
  dplyr::select(c(sample, LYMPH_NODES_EXAMINED_POSITIVE, NPI, stage, OS_MONTHS, OS_STATUS)) %>%
  dplyr::mutate(status = ifelse(clinicalEXP$OS_STATUS == "DECEASED",1,0)) %>%
  tibble::column_to_rownames('sample') %>%
  dplyr::select(-OS_STATUS)
colnames(clinicalEXP)[1:4] =  c("lymph", "npi", "stage", "time")
str(clinicalEXP)
# 'data.frame':	1904 obs. of  5 variables:
# $ lymph : int  1 5 8 1 0 1 0 2 0 6 ...
# $ npi   : num  4.04 6.03 6.03 5.04 3.05 ...
# $ stage : int  2 2 3 2 2 2 1 2 2 4 ...
# $ time  : num  47 20.4 138.1 119.8 101.2 ...
# $ status: num  1 1 0 0 0 0 0 0 0 1 ...

#preprocess clinicalCNA
clinicalCNA = clinicalCNA %>%
  tibble::rownames_to_column('sample') %>%
  dplyr::select(c(sample, LYMPH_NODES_EXAMINED_POSITIVE, NPI, stage, OS_MONTHS, OS_STATUS)) %>%
  dplyr::mutate(status = ifelse(clinicalCNA$OS_STATUS == "DECEASED",1,0)) %>%
  tibble::column_to_rownames('sample') %>%
  dplyr::select(-OS_STATUS)
colnames(clinicalCNA)[1:4] =  c("lymph", "npi", "stage", "time")
clinicalCNA$stage = as.character(clinicalCNA$stage)
str(clinicalCNA)
# 'data.frame':	2173 obs. of  5 variables:
#   $ lymph : int  10 0 3 24 1 3 0 0 0 1 ...
# $ npi   : num  6.04 2.04 5.04 6.07 4.05 ...
# $ stage : chr  "2" "1" "2" "2" ...
# $ time  : num  140.5 163.5 164.9 14.1 103.8 ...
# $ status: num  0 0 0 1 0 0 0 0 0 0 ...

#-------------------------------------o0o-------------------------------------#
# RUN DrGA
#-------------------------------------o0o-------------------------------------#

#RUN!!!!
DriverGeneAnalysis(exp = exp, clinicalEXP = clinicalEXP, timeEXP = clinicalEXP$time, statusEXP = clinicalEXP$status, 
                   datMODULE4 = cna,  cliMODULE4 = clinicalCNA, timeMODULE4 = clinicalCNA$time, statusMODULE4 = clinicalCNA$status)
