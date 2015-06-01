\name{WNT_Mouse_Expr}
\alias{WNT_Mouse_Expr}
\docType{data}
\title{
Gene Expression data pertaining to the GMTL mouse mode
}
\description{
Gene Expression data (formatted and normalized using BrainArray CDFs) for 32 samples of the GTML mouse model (see references for details). 
These samples were obtained from GSE33199 and GSE24628. 
}
\usage{data(WNT_Mouse_Expr)}
\format{
  A data frame with 17607 observations on the following 6 variables.
  \describe{
    \item{\code{GSE33199_GSM821968}}{a numeric vector}
    \item{\code{GSE33199_GSM821969}}{a numeric vector}
    \item{\code{GSE33199_GSM821970}}{a numeric vector}
    \item{\code{GSE24628_GSM607429}}{a numeric vector}
    \item{\code{GSE24628_GSM607436}}{a numeric vector}
    \item{\code{GSE24628_GSM607444}}{a numeric vector}
  }
}
\source{
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33199
http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE24628
}
\references{
Gendoo DMA, Smirnov P, Lupien M, Haibe-Kains B, (2014) "Personalized Diagnosis of Medulloblastoma Across Patients and Model Systems", \emph{Genomics}, \bold{14}(16):5158--5165.
}
\examples{
data(WNT_Mouse_Expr)
}
\keyword{datasets}
