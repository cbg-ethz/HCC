#THIS SCRIPT CAN BE USED TO REPLICATE Figure 3A:D
#This script can be run in a parallelized environment

##########################################
#If simulations are run via command line #
##########################################
#two arguments
#first argument is the number replicates
#second argument is the number of nodes in one time slice of a DBN
args = commandArgs(trailingOnly=TRUE)
#example: args<-c("10","10")
#when run via command line the first argument is FPrate*100, e.g.
#is the proportion of added false positive edges
#(as a percentage of all true positives) is 200%
#then args[1] should be 200
FPrate<-as.numeric(args[1])/100
#when run via command line the first argument is FNrate*100, e.g.
#is the proportion of removed or false negative edges
#(as a percentage of all true positives) is 50%
#then args[2] should be 50
FNrate<-as.numeric(args[2])/100


#define path to save final result
path<-""

#parameter rep replicate number
#parameter n number of continuous variables
#parameter k number of clusters
#parameter ssevc vector of length k, number of observations per cluster
#parameter algorithm "bnclust" for bnClustOmics, otherwise multiple other algorithms will be run
#parameter savedata logical, when false only generated datasets will be stored (e.g. to test outside of R)
#parameter FPrate the proportion of added false positive edges (as a percentage of all true positives) added to the simulated database
#parameter FNrate the proportion of false negative edges (as a percentage of all true positives) removed from the simulated database

simBNclustcore<-function(rep,n=100,
                         k=4,ssvec=c(150,100,50,20),
                         algorithm="bnclust",savedata=FALSE,
                         FPrate=FPrate,FNrate=FNrate,type="mixed") {
  library(BiDAG)
  library(gRbase)
  library(pcalg)
  library(mclust)
  library(bnClustOmics)
  library(clue)
  #
  source("simulations_code/generate.R")
  source("simulations_code/otheralgos.R")
  source("simulations_code/helpfns.R")
  source("simulations_code/clustfns.R")
  source("simulations_code/comparemodels.R")
  #
  mixedpar=list(nbin=20,avchildren=0.5,par1=0.1,par2=7)
  if(algorithm!="mcmcMAP") onlyother<-FALSE else onlyother<-TRUE
  res<-list()
  res$ROC<-NULL
  res$accuracy<-NULL
  sseed<-100+rep
  set.seed(sseed)
  dflocal<-NULL
  #generate mixture
  bnmixt<-genMixture(k=k, type="mixed",centersignal="medium",
                     sigma0=0.3, ssvec=ssvec, n=100, avpar=1, deltamu=20, lB=0.5, uB=1.5,
                     shdpct=30, randomseed=sseed, eqval=0,
                     mixedpar=mixedpar, signalroots=TRUE)
  metainfo<-list(rep=rep,seed=sseed,N=sum(bnmixt$ssvec)/k,dmu=20,shdpct=30,k=4)
  if(savedata) {
    res<-bnmixt
  } else {
    res$info<-bnmixt$info
    if(type=="mixed") {
      res$info$nbin<-bnmixt$nbin
    }
    edgepmat<-simedgepmat2(bnmixt,pf=2,FNrate=FNrate,FPrate=FPrate)
    if(!onlyother) {
      #learn mixture
      set.seed(sseed)
      databn<-makeOmicsObject(bnmixt,n,mixedpar$nbin)
      omicsobj<-bnInfo(databn,types=c("b","c"),omics=c("M","T"))
      if(algorithm=="bnclust") {
        bnres<-bnClustOmics::bnclustOmics(databn,omicsobj, blacklist=NULL, edgepmat=edgepmat,kclust=k,
                                          maxEM=10,startpoint = "mclustPCA",baseprob=3/(k+2),plus1it=4,epmatrix=TRUE)

      }
      relab<-checkmembership(k,bnmixt$membership,bnres$memb)$relabel
      bnres$memb<-relabMembership(bnres$memb, relab)
      bnres$lambdas<-relabLambdas(bnres$lambdas,relab)
      bnres$DAGs<-relabDAGs(bnres$DAGs,relab)
      if(!is.null(bnres$ep)) {
        bnres$ep<-relabDAGs(bnres$ep,relab)
      }

      res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="MAP",rep=rep,seed=sseed))
      res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="p",p=c(0.3,0.5,0.7,0.9,0.95,0.99)))

      dflocal<-rbind(dflocal,data.frame(clustaccuracy(bnmixt$membership,bnres$memb,k,ss=nrow(bnmixt$data)),algorithm=algorithm, metainfo))
    }

    res$accuracy<-rbind(res$accuracy,dflocal)
  }
  return(res)
}

#

#this code can be used to run multiple replicates of simulations in parallel
library(parallel)
nrep<-2 #numbr of replicates, in the manuscript nrep=50
rep<-c(1:nrep)
cl <- makeCluster(nrep+1)
outputClApply <- parallel::clusterApply(cl, rep, simBNclustcore,
                                        n=100,
                                        k=4,
                                        ssvec=c(150,100,50,20),
                                        algorithm="bnclust",
                                        savedata=FALSE,
                                        FPrate=FPrate,FNrate=FNrate)
stopCluster(cl)

#save the result
saveRDS(outputClApply,file=paste(path,"DBsim.rds",sep=""))

