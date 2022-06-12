#THIS SCRIPT CAN BE USED TO REPLICATE Figure 2A:C
#This script can be run in a parallelized environment

#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

#path to save the result
path<-""
base<-"featureSelection"
nrep<-as.numeric(args[7])

path<-""
base<-paste(paste("DMU",args[c(1,2,3,4,5,6)],collapse="",sep=""),sep="")
print(base)


#parameter rep replicate number
#arguments are passed via a vector of characters args

#args[1] n number of continuous variables in the simulated Bayesian networks
#args[2] k number of clusters (mixture components)
#args[3] ss number of generated observations per cluster
#args[4] deltamu a percentage different conditional means between networks in the mixture multiplied by 100
#args[5] shdpct SHD between  between networks in the mixture  as a percentage of edges in the fist network multiplied by 100
#args[6] algorithm  "bnclust" for bnClustOmics, otherwise multiple other algorithms will be run
#args[7] nrep number or replicates

#example args<-c("100", "3", "20", "20", "30", "bnclust", "3")

bnclustSimCore<-function(rep,args) {
  library('reticulate')

  # Set conda environment to use
  use_condaenv("r-reticulate", required = TRUE)

  library(MOFA)
  library(MOFAdata)
  library(MultiAssayExperiment)


  library(clue)
  library(BiDAG)
  library(gRbase)
  library(pcalg)
  library(mclust)
  library(CIMLR)
  library(SIMLR)
  library(iClusterPlus)
  library(bnClustOmics)
  #
  source("simulations_code/R/generate.R")
  source("simulations_code/R/otheralgos.R")
  source("simulations_code/R/helpfns.R")
  source("simulations_code/R/clustfns.R")
  source("simulations_code/R/comparemodels.R")
  source("simulations_code/R/simclust.R")
  #
  #
  res<-simBNclust(nrep=1,
                  type="mixed",
                  n=as.numeric(args[1]),
                  k=as.numeric(args[2]),
                  avpar=1,
                  ssvec=rep(as.numeric(args[3]),as.numeric(args[2])),
                  centersignal="medium",
                  sigma0=0.3,
                  deltamu=as.numeric(args[4]),
                  lB=0.5, uB=1.5, shdpct=as.numeric(args[5]), eqval=0,
                  randomseed=100+rep,
                  mixedpar=list(nbin=20,avchildren=0.5,par1=0.1,par2=7),
                  algorithm=args[6],
                  hardlimit=12,maxEM=6,plus1it=5,
                  p=0.5,
                  startpoint="mclustPCA",ROC=FALSE,
                  edgep=FALSE,savedata=FALSE,
                  path=NULL,base=NULL,recimlr=FALSE)
  res$accuracy$rep<-rep
  return(res)
}

#this code can be used to run 50 replicates of simulations in parallel
library(parallel)
rep<-c(1:nrep)
cl <- makeCluster(nrep+1)
outputClApply <- parallel::clusterApply(cl, rep, bnclustSimCore,args)
stopCluster(cl)
res<-outputClApply
res<-list()
res$info<-outputClApply[[1]]$info
res$accuracy<-Reduce('rbind',lapply(outputClApply,function(x)x$accuracy))
rownames(res$accuracy)<-1:nrow(res$accuracy)

#ADD LINE TO SAVE THE RESULT!
saveRDS(res,file=paste(path,base,".rds",sep=""))
