#THIS SCRIPT CAN BE USED TO REPLICATE Figure 2D
#This script can be run in a parallelized environment

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#path to save the result
path<-""
base<-"featureSelection"


#parameter rep replicate number
#arguments are passed via a vector of characters args
#args[1] number of continuous variables selected from 1000 in the full network
#args[2] shadcpt SHD between  between networks in the mixture  as a percentage of edges in the fist network multiplied by 100
#args[3] deltamu a percentage different conditional means between networks in the mixture multiplied by 100

#example
#for a quick test, very few features selected (30) -> runs faster
#args<-c("30","30","10")
#options used in the manuscript
#args<-c("150","10","0")
#args<-c("150","20","3")
#args<-c("150","20","4")
#args<-c("150","30","6")
#rep<-1
featureSelectionCore<-function(rep,args){

  #restart R
  #Then run the scripts, e.g. in R
  library('reticulate')

  # Set conda environment to use
  use_condaenv("r-reticulate", required = TRUE)
  library(MOFA)
  library(MultiAssayExperiment)

  library(CIMLR)
  library(SIMLR)
  library(clue)
  library(BiDAG)
  library(gRbase)
  library(pcalg)
  library(mclust)
  library(iClusterPlus)
  library(bnClustOmics)
  #
  source("simulations_code/generate.R")
  source("simulations_code/otheralgos.R")
  source("simulations_code/helpfns.R")
  source("simulations_code/clustfns.R")
  source("simulations_code/comparemodels.R")
  source("simulations_code/simclust.R")

  ss<-20
  k<-3
  n<-1000
  nbin<-100

  mofares<-FALSE
  sseed<-100*rep
  set.seed(sseed)
  mixttest<-genMixture(k=k,type="mixed",centersignal="medium",sigma0=0.3, eqval=0,
                       ssvec=rep(ss,k),n=n,avpar=1,deltamu=as.numeric(args[3]),lB=0.5, uB=1.5,shdpct=as.numeric(args[2]),
                       randomseed=sseed,mixedpar=list(nbin=nbin,avchildren=0.5,par1=0.01,par2=7, dist="b"))

  newmixt<-mixttest
  colnames(mixttest$data)<-paste("V",1:(n+nbin),sep="")
  rownames(mixttest$data)<-paste("S",1:(ss*k),sep="")
  nodesB<-paste("V",1:nbin,sep="")
  topnodes<-paste("V",topNodes(mixttest,as.numeric(args[1])),sep="")

  res_local<-list()

  res_local$kmeans<-acckmeans(mixttest,abs=FALSE,npca=5,PCA=TRUE)
  res_local$mclust<-accmclust(mixttest,abs=FALSE,npca=5,PCA=TRUE)
  res_local$hclust<-acchclust(mixttest,abs=FALSE,npca=5,PCA=TRUE)

  aMOFA<-try(accMOFA(mixttest,abs=FALSE))
  if(is.error(aMOFA)) {
    res_local$MOFA<-data.frame(precision=0,recall=0,F1=0,ARI=0)
  } else {
    res_local$MOFA<-aMOFA
    mofares<-TRUE
  }
  res_local$iclust<-acciclust(mixttest,abs=FALSE)
  res_local$CIMLR<-accCIMLR(mixttest,abs=FALSE)
  res_local$CIMLRco<-accCIMLRco(mixttest,abs=FALSE)

  commonDAG<-1*Reduce('|',mixttest$DAGs)
  colnames(commonDAG)<-rownames(commonDAG)<-paste("V",1:(n+nbin),sep="")
  colnames(mixttest$data)<-colnames(commonDAG)

  #MOFA
  MOFAobject<-try(accMOFA(mixttest,abs=FALSE,accuracy=FALSE))
  if(is.error(MOFAobject)) {
    MOFAtop<-NULL
    res_local$mofan<-0
  } else {
    MOFAtop<-getTopFeats2(MOFAobject,"T","all",as.numeric(args[1]))
    res_local$mofan<-length(MOFAtop)
  }

  #hybridNodes
  hybnodes<-hybridNodes(mixttest,nodesMOFA=MOFAtop, topn=as.numeric(args[1]),plus=FALSE)#
  bnfit<-clustSubset2(mixttest,nodesB,hybnodes,commonDAG,k)
  res_local$bn_hyb<-bnfit

  #new code for selected features for other algorithms
  mixttest_fs<-selectFeatures(mixttest,nodesB,hybnodes)
  res_local$kmeans_fs<-acckmeans(mixttest_fs,abs=FALSE,PCA=TRUE)
  res_local$mclust_fs<-accmclust(mixttest_fs,abs=FALSE,PCA=TRUE)
  res_local$hclust_fs<-acchclust(mixttest_fs,abs=FALSE,PCA=TRUE)

  aMOFA<-try(accMOFA(mixttest_fs,abs=FALSE))
  if(is.error(aMOFA)) {
    res_local$MOFA_fs<-data.frame(precision=0,recall=0,F1=0,ARI=0)
  } else {
    res_local$MOFA_fs<-aMOFA
    mofares<-TRUE
  }
  res_local$iclust_fs<-acciclust(mixttest_fs,abs=FALSE)
  res_local$CIMLR_fs<-accCIMLR(mixttest_fs,abs=FALSE)
  res_local$CIMLRco_fs<-accCIMLRco(mixttest_fs,abs=FALSE)
  res_local$SHD<-as.numeric(args[2])
  res_local$mu<-as.numeric(args[3])
  res_local$rep<-rep
  res_local$nf<-as.numeric(args[1])

  return(res_local)
}


#this code can be used to run multiple replicates of simulations in parallel
library(parallel)
nrep<-2 #in the manuscript nrep=50
rep<-c(1:nrep) #number or replicates
cl <- makeCluster(nrep+1) #number of cores
outputClApply <- parallel::clusterApply(cl, rep, featureSelectionCore,args)
stopCluster(cl)
res<-outputClApply

#ADD LINE TO SAVE THE RESULT!
saveRDS(res,paste(path,base,args[1],"SHD",args[2],"MU",args[3],".rds",sep=""))
