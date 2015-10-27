\name{PredictionsDistributionBoxplot}
\alias{PredictionsDistributionBoxplot}
\title{
Boxplot of MM2S Subtype Predictions for Given Samples
}
\description{
This function generates a boxplot of MM2S subtype predictions for samples of interest. 
Users are provided the option to save this plot as a PDF file. 
}
\usage{
PredictionsDistributionBoxplot(InputMatrix,pdf_output,pdfheight,pdfwidth)
}
\arguments{
  \item{InputMatrix}{Output object from calling MM2S.Human or MM2S.Mouse predictions}
  \item{pdf_output}{Option to save the boxplot as a PDF file}
  \item{pdfheight}{User-defined specification for PDF height size}
  \item{pdfwidth}{User-defined specification for PDF width size}
}
\value{
  Generated Pie of all the MM2S subtype predictions in a testing set. 
}
\references{
Gendoo, D. M., Smirnov, P., Lupien, M. & Haibe-Kains, B. Personalized diagnosis of medulloblastoma subtypes across patients and model systems. 
Genomics, doi:10.1016/j.ygeno.2015.05.002 (2015)

Manuscript URL: http://www.sciencedirect.com/science/article/pii/S0888754315000774
}
\author{Deena M.A. Gendoo}
\examples{
# First run raw expression data through MM2S
# load Mouse gene expression data for the potential WNT mouse model
# data(WNT_Mouse_Expr)
# SubtypePreds<-MM2S.mouse(InputMatrix=WNT_Mouse_Expr,xls_output=TRUE,parallelize=2)
# Generate Heatmap
# PredictionsDistributionBoxplot(InputMatrix=SubtypePreds,pdf_output=TRUE,pdfheight=5,pdfwidth=5)
}
\keyword{ heatmap }
