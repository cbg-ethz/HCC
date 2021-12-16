library(BiDAG)
library(gRbase)
library(pcalg)
library(mclust)
library(bnclustOmics)
#
source("bnclustSimCode/R/generate.R")
source("bnclustSimCode/R/otheralgos.R")
source("bnclustSimCode/R/simclust.R")
source("bnclustSimCode/R/helpfns.R")
source("bnclustSimCode/R/clustfns.R")
source("bnclustSimCode/R/comparemodels.R")
#
#
usepm<-TRUE #use or not penalization matrix
base<-"strfit"
path<-""
simsim<-simBNclust(nrep=50,
                   type="mixed",
                   n=100,
                   k=4, avpar=1,
                   ssvec=c(150,100,50,20),
                   centersignal="medium",
                   deltamu=20,
                   lB=0.5, uB=1.5, shdpct=30, eqval=0,
                   randomseed=100,
                   mixedpar=list(nbin=20,avchildren=0.5,par1=0.1,par2=7),
                   algorithm="mcmcMAP",
                   hardlimit=10,maxEM=6,plus1it=5,
                   p=0.5, startpoint="mclustPCA",ROC=TRUE,addalgos=FALSE,
                   onlyother=FALSE,edgep=usepm,savedata=FALSE,
                   path=path,base=base,errRate=0.1)
