\name{PredictionsBarplot}
\alias{PredictionsBarplot}
\title{
Stacked Barplot of MM2S Subtype Predictions for Given Samples
}
\description{
This function generates a stacked barplot of MM2S subtype predictions for samples of interest. 
Users are provided the option to save this heatmap as a PDF file. 
}
\usage{
PredictionsBarplot(InputMatrix,pdf_output,pdfheight,pdfwidth)
}
\arguments{
  \item{InputMatrix}{Matrix with samples in rows, and columns with MM2S percentage predictions for each subtype (Gr4,Gr3,SHH,WNT, and Normal)}
  \item{pdf_output}{Option to save the heatmap as a PDF file}
  \item{pdfheight}{User-defined specification for PDF height size}
  \item{pdfwidth}{User-defined specification for PDF width size}
}
\value{
  Generated Stacked Barplot of MM2S subtype predictions. Samples are in columns. Stacks are reflective of prediction percentages across MB subtypes for a given sample.  
}
\references{
Gendoo, D. M., Smirnov, P., Lupien, M. & Haibe-Kains, B. Personalized diagnosis of medulloblastoma subtypes across patients and model systems. 
Genomics, doi:10.1016/j.ygeno.2015.05.002 (2015)

Manuscript URL: http://www.sciencedirect.com/science/article/pii/S0888754315000774
}
\author{Deena M.A. Gendoo}
\examples{
# Generate heatmap from already-computed predictions for the GTML Mouse Model
## load computed MM2S predictions for GTML mouse model
data(GTML_Mouse_Preds)
## Generate Barplot
PredictionsBarplot(InputMatrix=GTML_Mouse_Preds, pdf_output=TRUE,pdfheight=5,pdfwidth=5)

## Not Run
# Generate heatmap after running raw expression data through MM2S
# load Mouse gene expression data for the potential WNT mouse model
# data(WNT_Mouse_Expr)
# SubtypePreds<-MM2S.mouse(InputMatrix=WNT_Mouse_Expr[2:3],xls_output=TRUE,parallelize=1)
# Generate Heatmap
# PredictionsBarplot(InputMatrix=SubtypePreds$Predictions,pdf_output=TRUE,pdfheight=5,pdfwidth=5)
}
\keyword{ heatmap }
