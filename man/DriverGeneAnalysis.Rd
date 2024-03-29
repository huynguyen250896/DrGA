% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DrGA.R
\docType{package}
\name{DriverGeneAnalysis}
\alias{DriverGeneAnalysis}
\title{DrGA: driver gene analysis in an automatic manner}
\usage{
DriverGeneAnalysis(organism, sources, methodCC,
       exp, clinicalEXP, timeEXP, statusEXP,
       datMODULE4, cliMODULE4, timeMODULE4, statusMODULE4,
       minClusterSize, verbose,
       NetworkType, hm_row_names)
}
\arguments{
\item{organism}{organism name. Organism names are constructed by concatenating the first letter of
the name and the family name. Example: human - \code{hsapiens}, mouse - \code{mmusculus}. Default is
\code{organism = "hsapiens"}}

\item{sources}{possible biological mechanisms allowed (e.g., Gene Ontology - \code{GO:BP}, \code{GO:MF},
\code{GO:CC}; \code{KEGG}; \code{REAC}; \code{TF}; \code{MIRNA}; \code{CORUM}; \code{HP}; \code{HPA};
\code{WP};… Please see the g:GOSt web tool for the comprehensive list and details on incorporated data
sources). Default is \code{sources = c("GO:BP", "KEGG")}}

\item{methodCC}{Correlation method type. Allowed values are \code{spearman} (default), \code{pearson},
\code{kendall}}

\item{exp}{a data frame or matrix. \code{exp} has its rows are samples and its columns are genes.
It is input data to serve to run the second and third modules.}

\item{clinicalEXP}{a data frame or matrix.
It includes its rows are samples, and its columns are clinical features of your choice.
Note that if users want to perform survival analysis, \code{clinicalEXP} must have two 
columns as overall survival time (continuous variable) and overall survival status (binary 
variable; usually coded 1 as death and 0 as live) of all the subjects.}

\item{timeEXP}{a vector of overall survival time. It is a column vector of \code{clinicalEXP}.}

\item{statusEXP}{a vector of overall survival time. It is a column vector of \code{clinicalEXP}.}

\item{datMODULE4}{a data frame or matrix. \code{datMODULE4} has its rows are samples and its columns are genes.
It is input data to serve to run the forth module.}

\item{cliMODULE4}{a data frame or matrix.
It includes its rows are samples, and its columns are clinical features of your choice.
Note that if users want to perform survival analysis, \code{cliMODULE4} must have two 
columns as overall survival time (continuous variable) and overall survival status (binary 
variable; usually coded 1 as death and 0 as live) of all the subjects.}

\item{timeMODULE4}{a vector of overall survival time. It is a column vector of \code{cliMODULE4}}

\item{statusMODULE4}{a vector of overall survival time. It is a column vector of \code{cliMODULE4}}

\item{minClusterSize}{Minimum cluster size. \code{minClusterSize = 10} is default.}

\item{verbose}{Default value is \code{TRUE}. A logical specifying whether to print details of analysis processes.}

\item{NetworkType}{network type. Allowed values are (unique abbreviations of) "unsigned", "signed",
"signed hybrid". Default value is \code{signed}}

\item{hm_row_names}{logical. If \code{hm_row_names = TRUE} (default value), gene names appear in rows of
the heatmap.  If due to the large number of driver genes leading to impossibly showing gene names in
rows of the heatmap, users can turn them off by \code{hm_row_names = FALSE}.}
}
\description{
DrGA is a novel R package that has been developed based on the idea of our recent driver gene analysis scheme. Its aim is to wholy automate the analysis process of driver genes at a low investment of time for this process by merging state-of-the-art statistical tools into one.
}
\examples{
DriverGeneAnalysis(exp = exp, clinicalEXP = clinicalEXP, timeEXP = clinicalEXP$time, statusEXP = clinicalEXP$status, datMODULE4 = cna,  cliMODULE4 = clinicalCNA, timeMODULE4 = clinicalCNA$time, statusMODULE4 = clinicalCNA$status)
}

\references{
Quang-Huy Nguyen, Duc-Hau Le. (2022). 
DrGA: cancer driver gene analysis in a simpler manner.
BMC Genomics, 23(1):86.
}  

\author{
Quang-Huy Nguyen
}
