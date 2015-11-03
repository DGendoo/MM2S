\name{PCARender}
\alias{PCARender}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PCA Projections of ssGSEA Rank Matrix for Test Samples Onto the MM2S Training Set
}
\description{
This function projects the ssGSEA ranked matrix for the Test Data onto the ssGSEA ranked matrix of the MM2S training set. 
The projections are done using ssGSEA ranked matrix containing the selected genesets common to both the Training set and test data. 
}
\usage{
PCARender(GSVAmatrixTesting,GSVAmatrixTraining)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GSVAmatrixTesting}{ Matrix of ranked enrichment scores for the tested datasets or data sample, containing the sample(s) in rows and genesets in columns }
  \item{GSVAmatrixTraining}{ Matrix of ranked enrichment scores for the training datasets, containing the samples in rows and genesets in columns }
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
  3 PDF files of projected test data onto the MM2S training set, using PCA (Principal Component Analysis) based on the selected genesets:
  PC1-PC2: Projection of the ssGSEA rank matrix from the testing set onto the training data, first and second principal component
  PC2-PC3: Projection of the ssGSEA rank matrix from the testing set onto the training data, second and third principal component
  Lattice: Lattice matrix with the projections onto PCA1-PC3
}
\references{
Gendoo, D. M., Smirnov, P., Lupien, M. & Haibe-Kains, B. Personalized diagnosis of medulloblastoma subtypes across patients and model systems. 
Genomics, doi:10.1016/j.ygeno.2015.05.002 (2015)

Manuscript URL: http://www.sciencedirect.com/science/article/pii/S0888754315000774
}
\author{
Deena M.A. Gendoo
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
## Running raw expression data through MM2S
# load Mouse gene expression data for the potential WNT mouse model
# data(WNT_Mouse_Expr)
# SubtypePreds<-MM2S.mouse(InputMatrix=WNT_Mouse_Expr[2:3],xls_output=TRUE,parallelize=1)
## Generate Heatmap
#   PCARender(GSVAmatrixTesting=SubtypePreds$RankMatrixTesting,
#         GSVAmatrixTraining=SubtypePreds$RankMatrixTraining)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ PCA }
