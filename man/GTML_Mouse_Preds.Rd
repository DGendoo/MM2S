\name{GTML_Mouse_Preds}
\alias{GTML_Mouse_Preds}
\docType{data}
\title{
Table of Pre-computed MM2S Predictions for 32 replicates of the GTML mouse model
}
\description{
Table containing percent confidence scores for 32 mouse model replicates, across 5 Medulloblastoma subtypes (WNT, SHH, Group3, Group4, and Normal)
}
\usage{data(GTML_Mouse_Preds)}
\format{
  A data frame with 32 observations on the following 5 variables.
  \describe{
    \item{\code{Group3}}{a numeric vector}
    \item{\code{Group4}}{a numeric vector}
    \item{\code{SHH}}{a numeric vector}
    \item{\code{WNT}}{a numeric vector}
    \item{\code{NORMAL}}{a numeric vector}
  }
}
\details{
Gene Expression data of the GTML model was obtained from GSE36594 and normalized as described in (references). 
The table is an example of MM2S output for that  set of samples. 
}
\references{
Gendoo, D. M., Smirnov, P., Lupien, M. & Haibe-Kains, B. Personalized diagnosis of medulloblastoma subtypes across patients and model systems. 
Genomics, doi:10.1016/j.ygeno.2015.05.002 (2015)

Manuscript URL: http://www.sciencedirect.com/science/article/pii/S0888754315000774
}
\seealso{
\code{\link[MM2S]{PredictionsHeatmap}}
}
\examples{
data(GTML_Mouse_Preds)
PredictionsHeatmap(InputMatrix=GTML_Mouse_Preds, pdf_output=TRUE,pdfheight=20,pdfwidth=5)
}
\keyword{datasets}
