#' @title DrGA: driver gene analysis in an automatic manner
#'
#' @description DrGA is a novel R package that has been developed based on the idea of our recent driver gene analysis scheme. Its aim is to wholy automate the analysis process of driver genes at a low investment of time for this process by merging state-of-the-art statistical tools into one.
#'
#' @docType package
#'
#' @author Quang-Huy Nguyen
#'
#' @usage DriverGeneAnalysis(organism, sources, methodCC,
#'        exp, clinicalEXP, timeEXP, statusEXP,
#'        datMODULE4, cliMODULE4, timeMODULE4, statusMODULE4,
#'        minClusterSize, verbose,
#'        NetworkType, hm_row_names)
#'
#' @param organism organism name. Organism names are constructed by concatenating the first letter of
#' the name and the family name. Example: human - \code{hsapiens}, mouse - \code{mmusculus}. Default is
#' \code{organism = "hsapiens"}
#'
#' @param sources possible biological mechanisms allowed (e.g., Gene Ontology - \code{GO:BP}, \code{GO:MF},
#' \code{GO:CC}; \code{KEGG}; \code{REAC}; \code{TF}; \code{MIRNA}; \code{CORUM}; \code{HP}; \code{HPA};
#' \code{WP};… Please see the g:GOSt web tool for the comprehensive list and details on incorporated data
#' sources). Default is \code{sources = c("GO:BP", "KEGG")}
#'
#' @param methodCC Correlation method type. Allowed values are \code{spearman} (default), \code{pearson},
#' \code{kendall}, \code{kendall})
#'
#' @param exp a data frame or matrix. \code{exp} has its rows are samples and its columns are genes.
#' It is input data to serve to run the second and third modules.
#'
#' @param clinicalEXP  a data frame or matrix.
#' It includes its rows are samples, and its columns are clinical features of your choice.
#' Notice that the clinical data must have two columns as overall survival time (continuous
#' variable) and overall survival status (binary variable; usually coded 1 as death and 0 as
#' live) of all the subjects.
#'
#' @param timeEXP a vector of overall survival time. It is a column vector of \code{clinicalEXP}.
#'
#' @param statusEXP a vector of overall survival time. It is a column vector of \code{clinicalEXP}.
#'
#' @param datMODULE4 a data frame or matrix. \code{datMODULE4} has its rows are samples and its columns are genes.
#' It is input data to serve to run the forth module.
#'
#' @param cliMODULE4 a data frame or matrix.
#' It includes its rows are samples, and its columns are clinical features of your choice.
#' Notice that the clinical data must have two columns as overall survival time (continuous
#' variable) and overall survival status (binary variable; usually coded 1 as death and 0 as
#' live) of all the subjects.
#'
#' @param timeMODULE4 a vector of overall survival time. It is a column vector of \code{clinicalMODULE4}
#'
#' @param statusMODULE4 a vector of overall survival time. It is a column vector of \code{clinicalMODULE4}
#'
#' @param minClusterSize Minimum cluster size. \code{minClusterSize = 10} is default.
#'
#' @param verbose Defaul value is \code{TRUE}
#'
#' @param NetworkType network type. Allowed values are (unique abbreviations of) "unsigned", "signed",
#' "signed hybrid". Defaul value is \code{signed}
#'
#' @param hm_row_names logical. If \code{hm_row_names = TRUE} (default value), gene names appear in rows of
#' the heatmap.  If due to the large number of driver genes leading to impossibly showing gene names in
#' rows of the heatmap, users can turn them off by \code{hm_row_names = FALSE}.
#'
#' @return NULL
#'
#' @examples DriverGeneAnalysis(exp = exp, clinicalEXP = clinicalEXP, timeEXP = clinicalEXP$time, statusEXP = clinicalEXP$status, datMODULE4 = cna,  cliMODULE4 = clinicalCNA, timeMODULE4 = clinicalCNA$time, statusMODULE4 = clinicalCNA$status)
#'
#' @export


DriverGeneAnalysis = function(organism = "hsapiens", sources = c("GO:BP", "KEGG"), methodCC="spearman",
                              exp=NULL, clinicalEXP=NULL, timeEXP=NULL,
                              statusEXP=NULL,
                              datMODULE4=NULL, cliMODULE4 = NULL, timeMODULE4 = NULL, statusMODULE4 = NULL,
                              minClusterSize = 10, verbose = T,
                              NetworkType = "signed", hm_row_names = T){
  now0 = Sys.time()
  seed = round(rnorm(1)*10^6)

  #Errors
  if(missing(exp)){
    stop("Error: exp is missing. \n")
  }
  if(missing(datMODULE4)){
    stop("Error: datMODULE4 is missing. \n")
  }

  if(missing(clinicalEXP)){
    stop("Error: clinical data dedicated to exp is missing. \n")
  }
  if(missing(cliMODULE4)){
    stop("Error: clinical data dedicated to datMODULE4 is missing. \n")
  }

  if(missing(timeEXP)){
    stop("Error: overall survival time of patients dedicated to exp is missing. \n")
  }
  if(missing(timeMODULE4)){
    stop("Error: overall survival time of patients dedicated to datMODULE4 is missing. \n")
  }

  if(missing(statusEXP)){
    stop("Error: overall survival status of patients dedicated to exp is missing. \n")
  }
  if(missing(statusMODULE4)){
    stop("Error: overall survival status of patients dedicated to datMODULE4 is missing. \n")
  }

  if(all(!(colnames(exp) %in% colnames(datMODULE4)))){
    stop("Error: make sure you put all identified driver genes into columns of exp and datMODULE4. \n")
  }

  if(all(rownames(exp) != rownames(clinicalEXP))){
    stop("Error: make sure patients sharing between exp and clinicalEXP are included and in exactly the same order. \n")
  }

  if(all(rownames(datMODULE4) != rownames(cliMODULE4))){
    stop("Error: make sure patients sharing between datMODULE4 and cliMODULE4 are included and in exactly the same order. \n")
  }

  # defined log function
  mlog <- if(!verbose) function(...){} else function(...){
    message(...)
    flush.console()
  }

  # defined computeQ function to automatically compute Q-value following the Benjamini-Hochberg procedure
  computeQ <- function(x){
    (x$P.value*nrow(x))/(x$rank)
  }

  # defined geneSA function
  geneSA = function(genename=NULL, event=NULL){
    #run SA
    df1=lapply(genename,

               function(x) {

                 formula <- as.formula(paste('Surv(time,event)~',as.factor(x)))
                 coxFit <- coxph(formula, data = event)
                 summary(coxFit)
               })

    cc = data.frame(My_name_is = paste("Huy ",1:length(df1)), HR=NA, confidence_intervals=NA, P.value=NA)

    for (i in c(1:length(df1))) {
      cc$HR[i] = round(df1[[i]][["coefficients"]][2],3) #hazard ratio
      cc$confidence_intervals[i] = paste(round(df1[[i]][["conf.int"]][[3]],3), "-", round(df1[[i]][["conf.int"]][[4]],3)) #95% CI
      cc$P.value[i] = df1[[i]][["logtest"]][3] #P-value
      rownames(cc)[i] =rownames(df1[[i]][["conf.int"]])
      order.pvalue = order(cc$P.value)
      cc = cc[order.pvalue,] #re-order rows following p-value
      cc$rank = c(1:length(df1)) #rank of P.value
      cc$Q.value = computeQ(cc) #compute Q-value
      rownames(cc) <- gsub("up","",rownames(cc)) #remove the word "up" in row names
    }
    cc = dplyr::select(cc, -c(rank, My_name_is)) #remove the 'rank' and 'No.' columns
    cc = cc %>% subset(P.value <= 0.05) #only retain Genes with P <=0.05
    cc = cc %>% subset(Q.value <= 0.05) #only retain Genes with Q <=0.05
    write.table(cc,"gene_SA.txt",sep = "\t", quote = FALSE)

    #warning
    writeLines("\nNOTE: \n*gene_SA.txt placed in your current working directory.\n*Please check to identify which genes significantly associated with prognostic value (survival rates of patients).\n*In any case, the numerator is up-expression level and the denominator is down-expression level. In other words, the down-expression level of genes considered is the reference.")
  }

  # defined computeC function
  computeC = function(data,var,x)
  {
    #implementation
    cc1 <- data.frame(My_name_is=paste("Huy", 1:ncol(data)),CC=NA ,P.value=NA)
    estimates = numeric(ncol(data))
    pvalues = numeric(ncol(data))
    for (i in c(1:ncol(data))) {
      cc=cor.test(data[,i],var[,x],
                  method = methodCC)
      cc1$CC[i]=cc$estimate
      cc1$P.value[i]=cc$p.value
      rownames(cc1) = colnames(data)[1:ncol(data)]
    }
    order.pvalue = order(cc1$P.value)
    cc1 = cc1[order.pvalue,] #order rows following p-value
    cc1$rank = rank(cc1$P.value) #re-order
    cc1$Q.value = computeQ(cc1) #compute Q-value
    cc1 = cc1 %>% subset(P.value <= 0.05) #only retain Genes with P <=0.05
    cc1 = cc1 %>% subset(Q.value <= 0.05) #only retain Genes with Q <=0.05
    cc1 = dplyr::select(cc1, -c(rank, My_name_is))
    cc1 = list(cc1 %>% subset(CC > 0),cc1 %>% subset(CC < 0)) # [1] cor coefficient > 0 - [2] cor coefficient <0
    return(cc1)}

  #-------------------------------------o0o-------------------------------------#
  #### MODULE 1: DrGA-EA: Enrichment Analysis
  #-------------------------------------o0o-------------------------------------#

  set.seed(seed)
  mlog("MODULE 1. DrGA-EA: Enrichment Analysis")
  cat("- Starting to perform enrichment analysis of individual driver genes using g:Profiler...", "\n")
  gostres <- gost(query = colnames(exp),
                  organism = organism, ordered_query = FALSE,
                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE,
                  measure_underrepresentation = FALSE, evcodes = FALSE,
                  user_threshold = 0.05, correction_method = "g_SCS",
                  domain_scope = "annotated", custom_bg = NULL,
                  numeric_ns = "", sources = sources, as_short_link = FALSE)
  #result
  writeLines("NOTE:\n*EnrichmentAnalysis.txt placed in your current working directory.\n*Please check to specifically see the full result of this process.")
  gostres$result <- apply(gostres$result,2,as.character); gostres$result= as.matrix(gostres$result)
  write.table(gostres$result,"EnrichmentAnalysis.txt", sep="\t", quote = F)

  #time difference
  timediff0 = Sys.time() - now0;
  mlog("Done in ", timediff0, " ", units(timediff0), ".\n")

  #-------------------------------------o0o-------------------------------------#
  #### MODULE 2: VinCor
  #-------------------------------------o0o-------------------------------------#

  mlog("MODULE 2. DrGA-Cor: Correlation")
  now = Sys.time()
  #### MODULE 2.1: Association of driver genes with Overall survival of patients
  # create event vector for RNA expression data
  #> median is up-regulated genes and < median is down regulated genes
  event_rna <- apply(t(exp),1, function(x) ifelse(x > median(x),"up","down"))
  event_rna <- as.data.frame(event_rna) #should be as data frame || rows: patients, columns: driver genes
  event <- statusEXP #numeric;  death = 1, survival = 0
  time = timeEXP #numeric;
  #add time and event columns of clinical_exp to event_rna
  event_rna <- cbind(event_rna,time) #time
  event_rna <- cbind(event_rna,event) #status

  #type of data
  event_rna <- as.data.frame(event_rna)
  event_rna$time = as.double(as.character(event_rna$time))
  event_rna$event = as.double(as.character(event_rna$event))

  # Identify which driver genes significantly associated with prognostic value
  set.seed(seed)
  cat("\n", "- Starting to perform the association analysis of individual driver genes with survival rates of patients...", "\n")
  geneSA(genename = colnames(exp), event=event_rna)

  #### MODULE 2.2: Association of driver genes with other clinical features
  #create the necessary df
  cor= exp %>% as.data.frame() #should be as data frame

  #create clinical data with only clinical features of our choice
  remove_OSstatus = statusEXP
  remove_OStime = timeEXP

  index_status = c() #empty vector
  index_time = c() #empty vector
  for (i in 1:ncol(clinicalEXP)){
    index_status[i] = identical(remove_OSstatus,clinicalEXP[,i])
    index_time[i] = identical(remove_OStime,clinicalEXP[,i])
  }

  featureEXP = clinicalEXP[, -c(which(index_status == TRUE), which(index_time == TRUE))]  %>%
    as.data.frame() #should be as data frame

  #####correlation between driver genes and lymph
  set.seed(seed)
  cat("\n", "- Starting to perform association analysis of individual driver genes with the remaining clinical features of your choice...", "\n")
  listCC=list()
  for (i in 1:length(names(featureEXP)) ) {
    listCC[i] = lapply(names(featureEXP)[i], computeC, data = cor, var = featureEXP)
    names(listCC)[i] <- names(featureEXP)[i]
  }; print(listCC)

  #time difference
  timediff = Sys.time() - now;
  mlog("Done in ", timediff, " ", units(timediff), ".\n")

  #-------------------------------------o0o-------------------------------------#
  #### MODULE 3: DrGA-WGCNA
  #-------------------------------------o0o-------------------------------------#

  mlog("MODULE 3. DrGA-WGCNA: Weighted Gene Correlation Network Analysis")
  now1 = Sys.time()

  #### MODULE 3.1: "Weighted Gene Correlation Network Analysis" construction
  # some important settings
  options(stringsAsFactors = FALSE);
  invisible(capture.output(allowWGCNAThreads())) ### Allowing parallel execution

  # Choose a set of soft-thresholding powers
  cat("- Starting to choose the soft-thresholding power...", "\n")
  # Choose a set of soft-thresholding powers
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  # Call the network topology analysis function
  invisible(capture.output(sft <- WGCNA::pickSoftThreshold(exp, powerVector = powers,
                                                           verbose = 5, networkType = NetworkType)))
  # Plot the results:
  sizeGrWindow(9, 5)
  cex1 = 0.9;
  # Scale-free topology fit index as a function of the soft-thresholding power
  print(plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
             xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
             main = paste("Scale independence")));
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red");
  sft$Power = as.numeric(readline('What value of soft-thresholding power β will you select based on the Figure? \n'))


  adjacency = WGCNA::adjacency(exp, power = sft$Power,
                               type = NetworkType);

  # Turn adjacency into topological overlap
  invisible(capture.output(TOM <- TOMsimilarity(adjacency, TOMType = NetworkType)));
  dissTOM = 1-TOM

  #hierichical clustering
  cat("\n","- Starting to seek the optimal agglomeration method for grouping driver genes into functional modules...", "\n")
  # methods to assess
  m <- c( "average", "single", "complete", "ward")
  names(m) <- c( "average", "single", "complete", "ward")
  # function to compute agglomerative coefficient
  set.seed(seed)
  ac <- function(x) {
    agnes(exp, method = x)$ac
  }

  #automatically choose the best agglomeration method
  agg_method = purrr::map_dbl(m, ac) # Agglomerative coefficient of each agglomeration method
  agg_method = as.data.frame(agg_method)
  rownames(agg_method)[4]  = "ward.D2"
  agg_method$method = rownames(agg_method)
  agg_method=agg_method[which.max(agg_method$agg_method),]
  cat(">>>>> The best agglomeration method identified in this step is:", agg_method$method, "\n")

  # Call the hierarchical clustering function
  geneTree = hclust(as.dist(dissTOM), method = agg_method$method);

  # We set the minimum module size at 10:
  # Module identification using dynamic tree cut:
  invisible(capture.output(dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                                                        minClusterSize = minClusterSize)));

  # Convert numeric lables into colors
  moduleColors = labels2colors(dynamicMods)
  cat("\n",">>>>> The number of driver genes assigned to each of colored modules is: ", "\n"); print(table(moduleColors))

  # Plot the dendrogram and colors underneath
  # Open a pdf file
  pdf("Dendro-MolduColor.pdf", width=5.5, height=4)
  # Create the plot
  plotDendroAndColors(geneTree, moduleColors, "Module Colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05,
                      main = "Driver Gene dendrogram and module colors")
  # Close the file
  dev.off()

  #### MODULE 3.2: Clinical feature-gene modules Assocation
  #Define numbers of genes and samples
  nGenes = ncol(exp);
  nSamples = nrow(exp);
  # Recalculate MEs with color labels
  MEs0 = moduleEigengenes(exp, moduleColors)$eigengenes
  MEs = orderMEs(MEs0)
  moduleTraitCor = cor(MEs, featureEXP, use = "p");
  moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

  sizeGrWindow(20,6)
  # Will display correlations and their p-values
  textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitPvalue, 1), ")", sep = "");
  dim(textMatrix) = dim(moduleTraitCor)
  par(mar = c(6, 8.5, 3, 3));
  # Display the correlation values within a heatmap plot
  # Open a pdf file
  pdf("Assoc-CliModul.pdf", width = 6,height = 5)
  # Plot
  labeledHeatmap(Matrix = moduleTraitCor,
                 xLabels = names(featureEXP),
                 yLabels = names(table(moduleColors)),
                 ySymbols = names(MEs),
                 colorLabels = FALSE,
                 colors = blueWhiteRed(50),
                 textMatrix = textMatrix,
                 setStdMargins = FALSE,
                 cex.text = 0.8,
                 zlim = c(-1,1),
                 main = paste("Module-clinical feature relationships"))
  # Close the file
  dev.off()

  #Identify hub-genes in each module: high intramodular connectivity ~ high kwithin
  cat("\n","- Starting to detect top 5 hub-genes in each discovered module...")
  connectivity = intramodularConnectivity(adjacency, moduleColors)
  con=list()
  for (i in 1:length(unique(moduleColors)) ){
    con[[i]]=connectivity[colnames(exp)[moduleColors==unique(moduleColors)[i]],]
    order.kWithin = order(con[[i]]$kWithin, decreasing = TRUE)
    con[[i]] = con[[i]][order.kWithin,] #order rows following kWithin
    con[[i]] = con[[i]][1:5,] #top 5 genes that have a high connectivity to other genes in each module
  }

  for (j in 1:length(con)){
    cat("\n",">>>> Top 5 hub genes identified in the",unique(moduleColors)[[j]],"module is:",rownames(con[[j]]), "\n")
  }

  #time difference
  timediff1 = Sys.time() - now1;
  mlog("Done in ", timediff1, " ", units(timediff1), "\n")

  #-------------------------------------o0o-------------------------------------#
  #### MODULE 4: DrGA-PS
  #-------------------------------------o0o-------------------------------------#

  #MODULE 4.1. Classification
  mlog("MODULE 4. DrGA-PS: Patient Stratification")
  now2 = Sys.time()

  # function to compute agglomerative coefficient
  cat("- Starting to re-seek the optimal agglomeration method for grouping individuals into distinct subgroups...", "\n")

  set.seed(seed)
  ac <- function(x) {
    agnes(datMODULE4, method = x)$ac
  }

  agg_method1 = purrr::map_dbl(m, ac) # Agglomerative coefficient of each agglomeration method
  agg_method1 = as.data.frame(agg_method1)
  agg_method1$method = rownames(agg_method1)
  agg_method1=agg_method1[which.max(agg_method1$agg_method),]
  cat(">>>>> The best agglomeration method identified in this step is:", agg_method1$method, "\n", "\n")

  #find the number of cluster
  cat("- Starting to seek the optimal number of patient subgroups...", "\n")

  #Dunn's index
  set.seed(seed)
  transposed_cna=as.data.frame(datMODULE4)
  v <- clValid::clValid(transposed_cna, 2:15, clMethods="hierarchical",
                        validation="internal", metric = "euclidean", method = agg_method1$method,
                        maxitems = rownames(datMODULE4))
  #Plot result
  sizeGrWindow(4, 4) #call WGCNA
  pdf("optimal-group-number.pdf", width = 4, height = 4)
  plot(v)
  dev.off()

  #result
  optimalnumber=optimalScores(v)
  if((identical(optimalnumber$Clusters[1],optimalnumber$Clusters[2])==TRUE) | (identical(optimalnumber$Clusters[2],optimalnumber$Clusters[3]) == TRUE) | (identical(optimalnumber$Clusters[1],optimalnumber$Clusters[3])==TRUE) | (identical(optimalnumber$Clusters[1],optimalnumber$Clusters[3])==TRUE & identical(optimalnumber$Clusters[2],optimalnumber$Clusters[3])==TRUE & identical(optimalnumber$Clusters[1],optimalnumber$Clusters[2])==TRUE)){
    optimalnumber=names(sort(summary(as.factor(optimalnumber$Clusters)), decreasing=T)[1])
    cat(">>>> the optimal number of patient subgroups identified in this step is:", optimalnumber, "subgroups of patients","\n")}
  else{
    stop(">>>> No optimal subgroup number can be found in this step \n")}

  # Cut tree into the identified optimal subgroup numbers
  set.seed(seed)
  hc_a <- agnes(datMODULE4, method = agg_method1$method)
  sub_grp= cutree(as.hclust(hc_a), k = optimalnumber)

  # Number of members in each cluster
  cat(">>>> The number of patients is distributed to each of the", optimalnumber,"identified subgroups is:", "\n"); print(table(sub_grp))

  ## make a named vector from the vector
  info =as.data.frame(sub_grp)
  info$patient = rownames(info)
  info <-info[order(info$sub_grp),]
  info = dplyr::select(info,-patient)
  colnames(info) <- c('groups')
  info$groups = as.character(info$groups)
  datMODULE4 = datMODULE4[rownames(info),] #change the order of column/patients of cna data following the variable 'info'

  ## Heatmap annotation
  ha <- ComplexHeatmap::columnAnnotation(df = info)

  #Plot heatmap
  pdf("heatmap.pdf", width=6, height=6)
  hm<-ComplexHeatmap::Heatmap(t(datMODULE4), name = "CNA scale",
                              show_row_names = hm_row_names, show_column_names = FALSE,
                              cluster_columns = FALSE,show_column_dend = FALSE,
                              show_row_dend = FALSE,top_annotation = ha, column_order = rownames(info),
                              row_names_gp = gpar(fontsize = 6.5), use_raster = TRUE)
  draw(hm)
  dev.off()

  #MODULE 4.2. comparision between the identified subgroups
  #survival rate
  #timeMODULE4:  num;
  #statusMODULE4:  num; death = 1, survival = 0
  survData = data.frame(timeMODULE4, statusMODULE4)
  colnames(survData)[1:2] = c("time", "status")
  survData$time = as.double(as.character(survData$time))
  rownames(survData)<- rownames(cliMODULE4)

  #run SA
  set.seed(seed)
  coxFit <- survival::coxph(
    Surv(time, status) ~ as.factor(sub_grp),
    data = survData,
    ties = "exact"
  )

  #message
  cat("\n","- Starting to perform a comparison between the identified", optimalnumber, "patient subgroups in term of survival rates...", "\n")
  pcox= summary(coxFit)$logtest[3]#Cox p-value
  cat(">>>> The Cox P-value gained from comparing patient outcomes between the identified", optimalnumber, "patient subgroups is: ", pcox[[1]])
  cat("\n", ">>>> And the Hazard ratio between the identified", optimalnumber, "patient subgroups is: "); print(exp(coxFit[["coefficients"]]))
  cat("With its 95% Confidence Interval is: ", paste(round(summary(coxFit)[["conf.int"]][[3]],3), "-", round(summary(coxFit)[["conf.int"]][[4]],3)))

  #remaining clinical feature
  #message
  cat("\n", "\n", "- Starting to perform comparisons between the identified", optimalnumber, "patient subgroups in terms of remaining clinical features...", "\n")

  set.seed(seed)
  #create clinical data with only clinical features of our choice
  remove_OSstatus1 = statusMODULE4
  remove_OStime1 = timeMODULE4

  index_status1 = c() #empty vector
  index_time1 = c() #empty vector
  for (i in 1:ncol(cliMODULE4)){
    index_status1[i] = identical(remove_OSstatus1,cliMODULE4[,i])
    index_time1[i] = identical(remove_OStime1,cliMODULE4[,i])
  }

  featureCNA = cliMODULE4[, -c(which(index_status1 == TRUE), which(index_time1 == TRUE))]  %>%
    as.data.frame() #should be as data frame

  featureCNA = featureCNA[rownames(datMODULE4),] #change the order of column/patients of featureCNA data following the variable 'info'
  featureCNA$groups = info$groups
  des=compareGroups::createTable(compareGroups::compareGroups(groups ~ ., data = featureCNA, method = NA))
  #word
  compareGroups::export2xls(des, file = "tableSTAT.xlsx", header.labels = c(p.overall = "p-value"))
  #message
  cat(">>>> The following are the remaining clinical features used and their own statistical description", "\n"); print(des$avail[,4])
  writeLines("NOTE:\n*tableSTAT.txt placed in your current working directory\n*Please check to observe the statistical differences in remaining clinical features between identified subgroups.")

  #time difference
  timediff2 = Sys.time() - now2;
  mlog("Done in ", timediff2, " ", units(timediff2), ".\n")
}
