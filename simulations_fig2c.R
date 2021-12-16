library(pcalg)
library(MOFA)
library(MultiAssayExperiment)
library(clue)
library(BiDAG)
library(gRbase)
library(pcalg)
library(mclust)
library(CIMLR)
library(SIMLR)
library(iClusterPlus)
library(bnclustOmics)
#
source("simulations_code/generate.R")
source("simulations_code/otheralgos.R")
source("simulations_code/simclust.R")
source("simulations_code/helpfns.R")
source("simulations_code/clustfns.R")
source("simulations_code/comparemodels.R")
#
k<-3 #change the number of clusters here
ss<-20
#this is fixed in fig2_c
deltamu<-20 # xi*10
#this is fixed in fig2_c
shdpct<-30 # eta*10
algorithm<-"mcmcMAP" #AO for other algorithms
benchmarkalgos<-FALSE #true for other algorithms
path<-""#where to save results
base<-"varyk"#basename of results file
res<-simBNclust(nrep=10,
                type="mixed",
                n=100,
                k=k,
                avpar=1,
                ssvec=rep(ss,k),
                centersignal="medium",
                sigma0=0.3,
                deltamu=20,
                lB=0.5, uB=1.5, shdpct=shdpct, eqval=0,
                randomseed=100,
                mixedpar=list(nbin=20,avchildren=0.5,par1=0.1,par2=7),
                algorithm=algorithm,
                hardlimit=10,maxEM=6,plus1it=4,
                p=0.5,
                startpoint="mclustPCA",ROC=FALSE,addalgos=benchmarkalgos,
                onlyother=benchmarkalgos,edgep=FALSE,savedata=FALSE,
                path=path,base=base,recimlr=FALSE)

#results
res$accuracy
