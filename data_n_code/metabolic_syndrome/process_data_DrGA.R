#library
devtools::install_github("huynguyen250896/DrGA", force = T) #NOTE: Download all dependencies of the tool!
libraries=c("DrGA", "dplyr", "survival", "tibble", "tidyr", "ComplexHeatmap", 
    'cluster', 'mclust', 'clValid', 'Biobase', 'annotate', 'GO.db', 
    'mygene', "dynamicTreeCut", "flashClust", "Hmisc", "WGCNA","purrr",
    "gprofiler2", "table1", "compareGroups")
lapply(libraries, require, character.only = TRUE)

#load raw file
exp = read.table("LiverFemale3600.csv", header = T, check.names = F, sep=",")
cli = read.table("ClinicalTraits.csv", header = T, check.names = F, sep=",", row.names = 1)

#remove missing gene names and duplicated genes in exp
exp = exp[which(exp$gene_symbol != "0"),]

dup = duplicated(exp$gene_symbol)
exp = exp[which(dup == FALSE),]

#turn exp1 into satisfactory format of DrGA
#DrGA requires data whose rows are samples and columns are genes.
exp = exp %>%
  dplyr::select(-c(substanceBXH, LocusLinkID, ProteomeID, cytogeneticLoc,
            CHROMOSOME, StartPosition, EndPosition)) %>%
  tibble::remove_rownames() %>%
  tibble::column_to_rownames('gene_symbol') %>%
  drop_na() %>% t()

#detect outliers
sampleTree = hclust(dist(exp), method = "average")
par(cex = 0.6);par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 12.2, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 12.2, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
exp = exp[keepSamples, ]

#match mouses that share between cli versus exp
cli = cli[cli$Mice %in% rownames(exp),]
cli = cli %>%
  remove_rownames() %>%
  tibble::column_to_rownames('Mice') %>% 
  dplyr::select(-c(Number, sex, Mouse_ID, Strain, DOB, parents, Western_Diet,
            Sac_Date, comments, Note))

#how the clinical traitsrelate to the sample dendrogram.
# Re-cluster samples
sampleTree2 = hclust(dist(exp), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(cli, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(cli),
                    main = "Sample dendrogram and trait heatmap")
#white means a low value, red a high value, and grey a missing entry.
#Only keep several clinical features with red color
cli = cli %>%
  dplyr::select(weight_g, length_cm, ab_fat,
                total_fat, UC, FFA, Glucose, 
                LDL_plus_VLDL)

#make sure that mice that share between exp and cli are included at their rows and
#in exactly the same order 
all(rownames(exp) == rownames(cli))
#[1] FALSE
exp = exp[rownames(cli), ]

#check dimension
dim(exp)
#[1] 134 2281
dim(cli)
#[1] 134  8

#-------------------------------------o0o-------------------------------------#
# RUN DrGA
#-------------------------------------o0o-------------------------------------#

#RUN!!!!
drga = DriverGeneAnalysis(exp = exp, clinicalEXP = cli,
                   datMODULE4 = exp,  cliMODULE4 = cli, organism = 'mmusculus',
                   hm_row_names = F)

