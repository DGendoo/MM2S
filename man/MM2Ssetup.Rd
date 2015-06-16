\name{MM2Ssetup}
\alias{MM2Ssetup}
\alias{HumanGMT}
\alias{MouseGMT}
\alias{MB_SampleInfo}
\alias{genesetHuman}
\alias{GenesetStatNormal}
\alias{GenesetStatWNT}
\alias{GenesetStatSHH}
\alias{GenesetStatGroup3}
\alias{GenesetStatGroup4}
\alias{Frozen_ES_Rank_Matrix}
\docType{data}

\title{
MM2S ssGSEA Training Matrix, MB Sample Group information for the training matrix, GMT files for Human and Mouse, and list of geneset names
}
\description{
This dataset contains the GMT files required to generate ssGSEA as part of MM2S along with a list of 674 geneset names that constitute the files, as well as the precomputed ssGSEA matrix on Human training data (347 human samples) as well as a list of the original MB subtype classifications for those samples. 
}
\usage{data(MM2Ssetup)}
\format{
\code{\link[MM2S]{MM2Ssetup}} is a data file containing four objects:
\describe{
\item{HumanGMT}{GMT file 694 Human Genesets, based on Gene Ontology Biological Processes (GO-BP) for MM2S subtype prediction}
\item{MouseGMT}{GMT file 694 Mouse Genesets, based on Gene Ontology Biological Processes (GO-BP) for MM2S subtype prediction}
\item{MB_SampleInfo}{Sample Names and Subtype for 347 Human Samples of the MM2S training set}
\item{genesetHuman}{List of 694 genesets that are common to mouse and human, and which are part of the GMT files}
\item{GenesetStatNormal}{List of 674 genesets and p-values (from the Wilcoxon test) for the Normal Subtype. Of these, the top 24 genesets specific to this group, and which overlap with the testing datasets, are selected for generating the ssGSEA rank matrix}
\item{GenesetStatWNT}{List of 674 genesets and p-values (from the Wilcoxon test) for the WNT Subtype. Of these, the top 24 genesets specific to this group, and which overlap with the testing datasets, are selected for generating the ssGSEA rank matrix}
\item{GenesetStatSHH}{List of 674 genesets and p-values (from the Wilcoxon test) for the SHH Subtype. Of these, the top 24 genesets specific to this group, and which overlap with the testing datasets, are selected for generating the ssGSEA rank matrix}
\item{GenesetStatGroup3}{List of 674 genesets and p-values (from the Wilcoxon test) for the Group3 Subtype. Of these, the top 24 genesets specific to this group, and which overlap with the testing datasets, are selected for generating the ssGSEA rank matrix}
\item{GenesetStatGroup4}{List of 674 genesets and p-values (from the Wilcoxon test) for the Group4 Subtype. Of these, the top 24 genesets specific to this group, and which overlap with the testing datasets, are selected for generating the ssGSEA rank matrix}
\item{Frozen_ES_Rank_Matrix}{Ranked list of 674 genesets per sample of the human training set, sorted in descending order of enrichment scores.}
}
}
\references{
Gendoo, D. M., Smirnov, P., Lupien, M. & Haibe-Kains, B. Personalized diagnosis of medulloblastoma subtypes across patients and model systems. 
Genomics, doi:10.1016/j.ygeno.2015.05.002 (2015)

Manuscript URL: http://www.sciencedirect.com/science/article/pii/S0888754315000774
}
\author{
Deena M.A. Gendoo
}

\examples{
data(MM2Ssetup)
}
\keyword{data}
