\name{MM2S.mouse}
\alias{MM2S.mouse}
\title{
MM2S Prediction of Mouse Medulloblastoma Samples
}
\description{
This function generates MM2S subtype predictions for Mouse samples of interest. 
Users are provided the option to save this the predictions as a XLS file. 
}
\usage{
MM2S.mouse(InputMatrix,xls_output,parallelize)
}
\arguments{
  \item{InputMatrix}{Matrix containing normalized gene expression data, with genes in rows and samples in columns. Gene IDs MUST be represented as Entres Gene IDs. }
  \item{xls_output}{Option to save MM2S subtype predictions as a XLS file.}
  \item{parallelize}{Option to set number of cores to run ssGSEA calculations in parallel. Default is 1 (no parallelization)}
}
\value{
\item{Predictions}{MM2S Percent Confidence Predictions of Human Subtypes (Group3, Group4, Normal, SHH, WNT) for a given Mouse sample.}
\item{MM2S_Subtype}{List of Sample names the designated Human MB subtype from MM2S classification.}
\item{RankMatrixTesting}{ssGSEA rank matrix of the test data, using selected genesets common between the test data and training set. These genesets are based on a filtered list using the top24 common genesets, for each subtype}
\item{RankMatrixTraining}{ssGSEA rank matrix of the trained MM2S human data, using selected genesets common between the test data and training set. These genesets are based on a filtered list using the top24 common genesets, for each subtype}
}
\references{
Gendoo, D. M., Smirnov, P., Lupien, M. & Haibe-Kains, B. Personalized diagnosis of medulloblastoma subtypes across patients and model systems. 
Genomics, doi:10.1016/j.ygeno.2015.05.002 (2015)

Manuscript URL: http://www.sciencedirect.com/science/article/pii/S0888754315000774
}
\author{
Deena M.A. Gendoo
}
\seealso{
\code{\link[MM2S]{MM2S.human}},\code{\link[MM2S]{PredictionsHeatmap}},\code{\link[MM2S]{PCARender}}
}
\examples{
# Load Mouse gene expression data 
# Data contains selected samples from a heterogenous WNT Mouse model
# data(WNT_Mouse_Expr)
# Generate Subtype Predictions
# SubtypePreds<-MM2S.mouse(InputMatrix=WNT_Mouse_Expr[2:3],xls_output=TRUE,parallelize=2)

## Not Run
## Generate Heatmap of Predictions
# PredictionsHeatmap(InputMatrix=SubtypePreds$Predictions, pdf_output=TRUE,pdfheight=5,pdfwidth=5)
## Generate projections of the selected genesets from Mouse model onto the training set, using PCA
# PCARender(GSVAmatrixTesting=SubtypePreds$RankMatrixTesting,
# GSVAmatrixTraining=SubtypePreds$RankMatrixTraining)
}
\keyword{ classification }
