\name{Human_Cho_Expr}
\alias{Human_Cho_Expr}
\docType{data}
\title{
Gene Expression data from selected samples of the Human CHO Medulloblastoma dataset
}
\description{
Gene Expression data (formatted and normalized using BrainArray CDFs) for 10 selected samples of the CHO dataset (see references for details). 
}
\usage{data(Human_Cho_Expr)}
\format{
  A data frame with 12129 observations on the following 10 variables.
  \describe{
    \item{\code{MD_269}}{WNT sample}
    \item{\code{MD_172}}{WNT sample}
    \item{\code{NC_1}}{Normal Cerebellum}
    \item{\code{NC_3}}{Normal Cerebellum}
    \item{\code{MD_151}}{SHH sample}
    \item{\code{MD_153}}{SHH sample}
    \item{\code{MD_268}}{Group4 sample}
    \item{\code{MD_161}}{Group4 sample}
    \item{\code{MD_242}}{Group3 sample}
    \item{\code{MD_75}}{Group3 sample}
  }
}
\source{
CHO dataset: Y.J. Cho, A. Tsherniak, P. Tamayo, S. Santagata, A. Ligon, H. Greulich, R. Berhoukim, V. Amani, L. Goumnerova, C.G. Eberhart, C.C. Lau, J.M. Olson, R.J. Gilbertson, A. Gajjar, O. Delattre, M. Kool, K. Ligon, M. Meyerson, J.P. Mesirov, S.L. Pomeroy, Integrative genomic analysis of medulloblastoma identifies a molecular subgroup that drives poor clinical outcome, J Clin Oncol, 29 (2011) 1424-1430.
}
\references{
Gendoo DMA, Smirnov P, Lupien M, Haibe-Kains B, (2014) "Personalized Diagnosis of Medulloblastoma Across Patients and Model Systems", \emph{Genomics}, \bold{14}(16):5158--5165.
}
\examples{
data(Human_Cho_Expr)
}
\keyword{datasets}
