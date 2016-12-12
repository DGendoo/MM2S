# Deena M.A. Gendoo
# Perform ssGSEA analysis using the GSVA function
# Input: GMT files (Human and Mouse), Expression Matrix of Mouse, Expression Matrix of Human
# October 16, 2015

#GOAL:
# Perform a 10-fold cross validation on the training set to get estimate of
# - best K parameter in KKNN
# - best number of top genesets to use in the algorithm.

##########################################################
# LIBRARIES
##########################################################
library(WriteXLS)
library(gplots)
library(VennDiagram)
library(GSVA)
library(gdata)
library(genefu)
library(parallel)
library(survcomp)

##########################################################
# FUNCTION TO CALCULATE MCC FOR THE KKNN CLASSIFIER, courtesy of BHK
# This cofficient is calculated by taking the confusion matrix from the KNN classification
##########################################################
`mcc` <- 
  function(ct, nbcat=nrow(ct)) {
    if(nrow(ct) != ncol(ct)) { stop("the confusion table should be square!") }
    if(!(sum(ct)==sum(diag(ct))) &&  (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 2, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1)))) || (length(which(apply(ct, 2, sum) == 0)) == (nbcat-1) & ((length(which(apply(ct, 1, sum) == 0)) != (nbcat-1)) | (length(which(apply(ct, 1, sum) == 0)) == (nbcat-1)) & sum(diag(ct)) == 0))) { ct <- ct + matrix(1, ncol=nbcat, nrow=nbcat) } ### add element to categories if nbcat-1 predictive categories do not contain elements. Not in case where all are correct!
    
    if(sum(ct, na.rm=TRUE) <= 0) { return(NA) }
    
    myx <- matrix(TRUE, nrow=nrow(ct), ncol=ncol(ct))
    diag(myx) <- FALSE
    if(sum(ct[myx]) == 0) { return(1) }
    myperf <- 0
    for(k in 1:nbcat) {
      for(m in 1:nbcat) {
        for(l in 1:nbcat) {
          myperf <- myperf + ((ct[k, k] * ct[m, l]) - (ct[l, k] * ct[k, m]))
        }
      }
    }
    aa <- 0
    for(k in 1:nbcat) {
      cc <- 0
      for(l in 1:nbcat) { cc <- cc + ct[l, k] }
      dd <- 0
      for(f in 1:nbcat) {
        for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[g, f] } }
      }
      aa <- aa + (cc * dd)
    }
    bb <- 0
    for(k in 1:nbcat) {
      cc <- 0
      for(l in 1:nbcat) { cc <- cc + ct[k, l] }
      dd <- 0
      for(f in 1:nbcat) {
        for(g in 1:nbcat) { if(f != k) { dd <- dd + ct[f, g] } }
      }
      bb <- bb + (cc * dd)
    }
    
    myperf <- myperf / (sqrt(aa) * sqrt(bb))
    return(myperf)
  }

###################################################################
# Load Input, Annotate the data and prepare
###################################################################
HumanGMT<-HumanGMT_20_100 ## 694 genesets
MB_SampleGroups<-data.frame(Sample=rownames(HumanGSVAFull),GROUP=HumanGSVAFull$Group)
MB_SampleInfo<-MB_SampleGroups
colnames(MB_SampleInfo)<-c("Sample_ID","subtype")
colnames(MB_SampleGroups)<-c("Sample_ID","subtype")
HumanGSVA<-HumanGSVAFull 
table(HumanGSVA$Group) 
genesetHuman<-colnames(HumanGSVA[-ncol(HumanGSVA)])
Northcott<-HumanGSVA

############################################################################
############################################################################

HumanGSVAOriginal<-HumanGSVA

Iteration <- sample(1:55555, 100,replace=F)

AvgTraining<-NULL
AvgTesting<-NULL
SumConfusion<-NULL

for (Iteration in 1:100)
{
  set.seed(Iteration)
  
  message("Iteration Number: ",Iteration)
  
  library(kknn)
  
  #randomize
  ListofSamples<-sample(MB_SampleGroups$Sample_ID,replace=F)
  # fold select
  TestingSubset<-NULL
  TestingSubset[[1]]<-ListofSamples[1:35]
  TestingSubset[[2]]<-ListofSamples[36:70]
  TestingSubset[[3]]<-ListofSamples[71:105]
  TestingSubset[[4]]<-ListofSamples[106:140]
  TestingSubset[[5]]<-ListofSamples[141:175]
  TestingSubset[[6]]<-ListofSamples[176:210]
  TestingSubset[[7]]<-ListofSamples[211:245]
  TestingSubset[[8]]<-ListofSamples[246:280]
  TestingSubset[[9]]<-ListofSamples[281:315]
  TestingSubset[[10]]<-ListofSamples[316:347]
  
  MCC_LIST_TRAINING<-NULL
  MCC_LIST_TESTING<-NULL
  MCC_TRAINING_ByFOLD<-NULL
  MCC_TESTING_ByFOLD<-NULL
  
  ConfusionMatrix<-NULL
  
  for (FOLD in 1:10) #10-Fold cross validation
  { 
    SubsetTest<-as.character(TestingSubset[[FOLD]])    
    SubsetTrain<-as.character(ListofSamples[which(SubsetTest!=ListofSamples)])
    
    Northcott<-HumanGSVAOriginal[SubsetTrain,] 

    # groups 
    Normal<-(Northcott[Northcott$Group =="NORMAL",])
    Group3<-(Northcott[Northcott$Group =="Group3",])
    Group4<-(Northcott[Northcott$Group =="Group4",])
    SHH<-(Northcott[Northcott$Group =="SHH",])
    WNT<-(Northcott[Northcott$Group =="WNT",])
    Normal<-Normal[,-ncol(Normal)]
    Group3<-Group3[,-ncol(Group3)]
    Group4<-Group4[,-ncol(Group4)]
    SHH<-SHH[,-ncol(SHH)]
    WNT<-WNT[,-ncol(WNT)]
    NonNormal<-(Northcott[!Northcott$Group =="NORMAL",])
    NonGroup3<-(Northcott[!Northcott$Group =="Group3",])
    NonGroup4<-(Northcott[!Northcott$Group =="Group4",])
    NonSHH<-(Northcott[!Northcott$Group =="SHH",])
    NonWNT<-(Northcott[!Northcott$Group =="WNT",])
    NonNormal<-NonNormal[,-ncol(NonNormal)]
    NonGroup3<-NonGroup3[,-ncol(NonGroup3)]
    NonGroup4<-NonGroup4[,-ncol(NonGroup4)]
    NonSHH<-NonSHH[,-ncol(NonSHH)]
    NonWNT<-NonWNT[,-ncol(NonWNT)]
    
    #Get genesets 
    geneset<-colnames(Northcott[,-ncol(Northcott)])
    GenesetStatSHH<-NULL
    GenesetStatNormal<-NULL
    GenesetStatGroup3<-NULL
    GenesetStatGroup4<-NULL
    GenesetStatWNT<-NULL
    for (count in 1:length(geneset))
    {
      GenesetStatNormal[geneset[count]]=(wilcox.test(Normal[,count],NonNormal[,count],paired=FALSE,conf.level=0.95))$p.value
      GenesetStatGroup3[geneset[count]]=(wilcox.test(Group3[,count],NonGroup3[,count],paired=FALSE,conf.level=0.95))$p.value
      GenesetStatGroup4[geneset[count]]=(wilcox.test(Group4[,count],NonGroup4[,count],paired=FALSE,conf.level=0.95))$p.value
      GenesetStatSHH[geneset[count]]=(wilcox.test(SHH[,count],NonSHH[,count],paired=FALSE,conf.level=0.95))$p.value
      GenesetStatWNT[geneset[count]]=(wilcox.test(WNT[,count],NonWNT[,count],paired=FALSE,conf.level=0.95))$p.value
    }
  
    FeatureSelection<-as.list(NULL)
    for(NumFeature in 1:25)
    {
      FeatureSelection[[NumFeature]]<-unique(
        c(names(sort(GenesetStatSHH,decreasing=FALSE))[1:NumFeature],
          names(sort(GenesetStatNormal,decreasing=FALSE))[1:NumFeature],
          names(sort(GenesetStatGroup4,decreasing=FALSE))[1:NumFeature],
          names(sort(GenesetStatGroup3,decreasing=FALSE))[1:NumFeature]))
    }  
    
    Human_GSVA_Matrix<-NULL
    
    for(NumFeature in 1:25)
    {
      HumanGSVA<-HumanGSVAOriginal 
      HumanGSVA<-HumanGSVA[,FeatureSelection[[NumFeature]]] 
      HumanGSVA<-t(HumanGSVA) 
      genesetHuman<-rownames(HumanGSVA) 
      
      Matrix_RANK_Human<-data.frame(genesetHuman) 
      for(sample in 1:ncol(HumanGSVA)) 
      {
        TempRankHuman<-sort(HumanGSVA[,sample],decreasing=TRUE)
        Matrix_RANK_Human[,(colnames(HumanGSVA)[sample])]<-match(Matrix_RANK_Human$geneset,names(TempRankHuman))
      }
      rownames(Matrix_RANK_Human)<-Matrix_RANK_Human$geneset 
      Matrix_RANK_Human<-Matrix_RANK_Human[,-1]
      Human_GSVA_Matrix[[NumFeature]]<-Matrix_RANK_Human
    }
    
    for(NumFeature in 1:25)
    {
      TrainSet<-t(Human_GSVA_Matrix[[NumFeature]][,SubsetTrain]) 
      TestSet<-t(Human_GSVA_Matrix[[NumFeature]][,SubsetTest]) 
      
      geneset<-FeatureSelection[[NumFeature]] 
             
      TrainSet<-as.data.frame(TrainSet)
      TestSet<-as.data.frame(TestSet)
      TrainSet[,"Group"]<- MB_SampleInfo$subtype[match((rownames(TrainSet)),as.character(MB_SampleInfo$Sample_ID))]
      TestSet[,"Group"]<- MB_SampleInfo$subtype[match((rownames(TestSet)),as.character(MB_SampleInfo$Sample_ID))]  
      
      TrainKKNN<-train.kknn(formula=TrainSet$Group ~ .,data=TrainSet[,-ncol(TrainSet)],kmax=25,distance=1,kernel="rectangular") 
      
      MCC<-NULL
      ConfusionMatrix<-NULL
      for(CountOfK in 1:25)
      {
        ConfusionMatrix<-table(TrainSet$Group,TrainKKNN$fitted.values[[CountOfK]]) #Compare for every K generated from the training range above
        MCC[CountOfK]<-mcc(ct=ConfusionMatrix, nbcat=5)
      }
      MCC_LIST_TRAINING[[NumFeature]]<-MCC
      
      MCC<-NULL
      ConfusionMatrix<-NULL
      for(CountOfK in 1:25)
      {
        TestKKNN<-kknn(formula = TrainSet$Group ~ ., TrainSet[,-ncol(TrainSet)], TestSet[,-ncol(TestSet)], na.action = na.omit(),k = CountOfK, distance = 1, kernel = "rectangular", scale=TRUE)
        ConfusionMatrix<-table(TestSet$Group,TestKKNN$fitted.values)
        MCC[CountOfK]<-mcc(ct=ConfusionMatrix, nbcat=5)
      }
      MCC_LIST_TESTING[[NumFeature]]<-MCC 
      
    }
  
    MCC_TRAINING_MATRIX <- data.frame(matrix(unlist(MCC_LIST_TRAINING), nrow=25, byrow=T))
    MCC_TESTING_MATRIX <- data.frame(matrix(unlist(MCC_LIST_TESTING), nrow=25, byrow=T))
    
    MCC_TRAINING_ByFOLD[[FOLD]]<-MCC_TRAINING_MATRIX
    MCC_TESTING_ByFOLD[[FOLD]]<-MCC_TESTING_MATRIX
  }  

  AvgTraining[[Iteration]]<-((MCC_TRAINING_ByFOLD[[1]]+MCC_TRAINING_ByFOLD[[2]]+MCC_TRAINING_ByFOLD[[3]]+MCC_TRAINING_ByFOLD[[4]]+MCC_TRAINING_ByFOLD[[5]]
                +MCC_TRAINING_ByFOLD[[6]]+MCC_TRAINING_ByFOLD[[7]]+MCC_TRAINING_ByFOLD[[8]]+MCC_TRAINING_ByFOLD[[9]]+MCC_TRAINING_ByFOLD[[10]])/10)
  
  AvgTesting[[Iteration]]<-((MCC_TESTING_ByFOLD[[1]]+MCC_TESTING_ByFOLD[[2]]+MCC_TESTING_ByFOLD[[3]]+MCC_TESTING_ByFOLD[[4]]+MCC_TESTING_ByFOLD[[5]]
               +MCC_TESTING_ByFOLD[[6]]+MCC_TESTING_ByFOLD[[7]]+MCC_TESTING_ByFOLD[[8]]+MCC_TESTING_ByFOLD[[9]]+MCC_TESTING_ByFOLD[[10]])/10)
  
}  

AvgTestingSum_100Iterations<-(Reduce('+', AvgTesting))/100

write.table(AvgTestingSum_100Iterations,file="NewMCCMatrix_100Iterations.txt")

