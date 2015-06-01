\name{MM2S_Human}
\alias{MM2S_Human}
\title{
MM2S Prediction of Human Medulloblastoma Samples
}
\description{
This function generates MM2S subtype predictions for human samples of interest. 
Users are provided the option to save this the predictions as a XLS file. 
}
\usage{
MM2S_Human(InputMatrix,xls_output,parallelize)
}
\arguments{
  \item{InputMatrix}{Matrix containing normalized gene expression data, with genes in rows and samples in columns. Gene IDs MUST be represented as Entres Gene IDs. }
  \item{xls_output}{Option to save MM2S subtype predictions as a XLS file.}
  \item{parallelize}{Option to set number of cores to run ssGSEA calculations in parallel. Default is 1 (no parallelization)}
}
\value{
\item{Predictions}{MM2S Percent Confidence Predictions of Human Subtypes (Group3, Group4, Normal, SHH, WNT) for a given sample.}
\item{MM2S_Subtype}{List of Sample names the designated Human MB subtype from MM2S classification.}
\item{RankMatrixTesting}{ssGSEA rank matrix of the test data, using selected genesets common between the test data and training set. These genesets are based on a filtered list using the top24 common genesets, for each subtype}
\item{RankMatrixTraining}{ssGSEA rank matrix of the trained MM2S human data, using selected genesets common between the test data and training set. These genesets are based on a filtered list using the top24 common genesets, for each subtype}
}
\references{
Gendoo DMA, Smirnov P, Lupien M, Haibe-Kains B, (2014) "Personalized Diagnosis of Medulloblastoma Across Patients and Model Systems", \emph{Genomics}, \bold{14}(16):5158--5165.
}
\author{
Deena M.A. Gendoo
}
\seealso{
\code{\link[MM2S]{MM2S_Mouse}},\code{\link[MM2S]{PredictionsHeatmap}},\code{\link[MM2S]{PCARender}}
}
\examples{
## load Human gene expression data containing selected samples from GSE37418
data(Human_Cho_Expr)
SubtypePreds<-MM2S_Human(InputMatrix=Human_Cho_Expr,xls_output=TRUE,parallelize=1)
## Generate Heatmap of Predictions
PredictionsHeatmap(InputMatrix=SubtypePreds$Predictions, pdf_output=TRUE,pdfheight=5,pdfwidth=5)
}
\keyword{ classification }
