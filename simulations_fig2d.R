use_condaenv("r-reticulate", required = TRUE)
library(MOFA)
library(MOFAdata)
library(MultiAssayExperiment)
library(CIMLR)
library(SIMLR)
library(mogsa)
library(clue)
library(BiDAG)
library(gRbase)
library(pcalg)
library(mclust)
library(iClusterPlus)
library(bnclustOmics)
#
source("bnclustSimCode/R/generate.R")
source("bnclustSimCode/R/otheralgos.R")
source("bnclustSimCode/R/simclust.R")
source("bnclustSimCode/R/helpfns.R")
source("bnclustSimCode/R/clustfns.R")
source("bnclustSimCode/R/comparemodels.R")

glob<-list()
glob$res<-NULL
glob$nint<-NULL
glob$vnarMOFA<-vector()
#glob$nvarmo<-vector() for testing moCluster
ss<-20
k<-3
n<-1000
nbin<-100
args<-c(150,20,1,70,4) # c(150,30,1,70,6), c(150,10,1,70,2) 
path<-""
base<-"featsel1000"
#run 50 runs for each value of args, 5 for reasonable time
for(i in 1:5) {
  
  mofares<-FALSE
  sseed<-100*i
  set.seed(sseed)
  mixttest<-genMixture(k=k,type="mixed",centersignal="medium",sigma0=0.3, eqval=0,
                       ssvec=rep(ss,k),n=n,avpar=1,deltamu=as.numeric(args[5]),lB=0.5, uB=1.5,shdpct=as.numeric(args[2]),
                       randomseed=sseed,mixedpar=list(nbin=nbin,avchildren=0.5,par1=as.numeric(args[3])/10,par2=as.numeric(args[4])/10, dist="b"))
  
  newmixt<-mixttest
  colnames(mixttest$data)<-paste("V",1:(n+nbin),sep="")
  rownames(mixttest$data)<-paste("S",1:(ss*k),sep="")
  nodesB<-paste("V",1:nbin,sep="")
  topnodes<-paste("V",topNodes(mixttest,as.numeric(args[1])),sep="")
  
  res_local<-list()
  nint_local<-list()
  
  res_local$kmeans<-acckmeans(mixttest,abs=FALSE,npca=5,PCA=TRUE)[[2]]
  res_local$mclust<-accmclust(mixttest,abs=FALSE,npca=5,PCA=TRUE)[[2]]
  res_local$hclust<-acchclust(mixttest,abs=FALSE,npca=5,PCA=TRUE)[[2]]
  
  aMOFA<-try(accMOFA(mixttest,abs=FALSE))
  if(is.error(aMOFA)) {
    res_local$MOFA<-(-1)
  } else {
    res_local$MOFA<-aMOFA[[2]]
    mofares<-TRUE
  }
  res_local$iclust<-acciclust(mixttest,abs=FALSE)[[2]]
  res_local$CIMLR<-accCIMLR(mixttest,abs=FALSE)[[2]]
  res_local$CIMLRco<-accCIMLRco(mixttest,abs=FALSE)[[2]]
  
  commonDAG<-1*Reduce('|',mixttest$DAGs)
  colnames(commonDAG)<-rownames(commonDAG)<-paste("V",1:(n+nbin),sep="")
  colnames(mixttest$data)<-colnames(commonDAG)
  
  bnfit<-clustSubset2(mixttest,nodesB,topnodes,commonDAG,k)
  nint_local$top<-bnfit[[1]]
  res_local$bn_top<-bnfit[[2]]
  
  MOFAobject<-try(accMOFA(mixttest,abs=FALSE,accuracy=FALSE))
  if(is.error(MOFAobject)) {
    MOFAtop<-NULL
    glob$nvarMOFA[i]<-0
    nint_local$MOFA<-0
    res_local$bn_MOFA<-(-1)
  } else {
    MOFAtop<-getTopFeats(MOFAobject,"T","all",as.numeric(args[1]))
    glob$nvarMOFA[i]<-length(MOFAtop)
    bnfit<-clustSubset2(mixttest,nodesB,MOFAtop,commonDAG,k)
    nint_local$MOFA<-bnfit[[1]]
    res_local$bn_MOFA<-bnfit[[2]]
  }
  
  #bnclustOmics results using hybridNodes
  hybnodes<-hybridNodes(mixttest,nodesMOFA=MOFAtop, topn=as.numeric(args[1]),plus=FALSE)
  bnfit<-clustSubset2(mixttest,nodesB,hybnodes,commonDAG,k)
  nint_local$hyb<-bnfit[[1]]
  res_local$bn_hyb<-bnfit[[2]]
  
  glob$res<-rbind(glob$res,res_local)
  glob$nint<-rbind(glob$nint,nint_local)
  saveRDS(glob,file=paste(path,base,".rds",sep=""))
  
}
