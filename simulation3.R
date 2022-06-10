#THIS SCRIPT CAN BE USED TO REPLICATE XX
#A LINE WITH A PATH TO SOVE THE RESULT MUST BE ADDED IN THE END

# rep=1
# type="mixed"
# n=20
# k=4
# avpar=1
# ssvec=c(150,100,50,20)
# centersignal="medium"
# deltamu=20
# lB=0.5
# uB=1.5
# shdpct=30
# eqval=0
# randomseed=100
# mixedpar=list(nbin=20,avchildren=0.5,par1=0.1,par2=7)
# algorithm="mcmcMAP"
# hardlimit=10
# maxEM=6
# plus1it=5
# p=0.5
# startpoint="mclustPCA"
# ROC=TRUE
# addalgos=FALSE
# onlyother=FALSE
# edgep=TRUE
# savedata=FALSE
# path=NULL
# base=NULL
# sigma0=0.3
#FPrate<-0.1
#FNrate<-0.5


##########################################
#If simulations are run via command line #
##########################################
#first argument is the number replicates
#second argument is the number of nodes in one time slice of a DBN
args = commandArgs(trailingOnly=TRUE)
FPrate<-as.numeric(args[1])/100
FNrate<-as.numeric(args[2])/100


#define path to save final result
path<-"bnclustsim/"
base<-"DB"

simBNclustcore<-function(rep,type, n=50, k=3, avpar=1, ssvec=rep(100,k), centersignal="strong",
                     sigma0=0.3, deltamu=20,lB=0.5, uB=1.5, shdpct=20, eqval=0,
                     randomseed=100,mixedpar=list(nbin,avchildren,freqmin=0.1,freqmax=0.5),
                     algorithm=c("mcmcMAP","mcmcsample","ges"),hardlimit=10,maxEM=6,plus1it=3,p=0.5,
                     startpoint=c("random","mclust","mclustPCA"),ROC=TRUE,addalgos=TRUE,basename=NULL,path=NULL,
                     onlyother=FALSE,edgep=FALSE,savedata=FALSE,recimlr=FALSE,FPrate=0,FNrate=0) {
  library(BiDAG)
  library(gRbase)
  library(pcalg)
  library(mclust)
  library(bnclustOmics)
  library(clue)
  #
  source("/Users/polinasuter/Downloads/HCC/submit/simulations_code/generate.R")
  source("/Users/polinasuter/Downloads/HCC/submit/simulations_code/otheralgos.R")
  source("/Users/polinasuter/Downloads/HCC/submit/simulations_code/helpfns.R")
  source("/Users/polinasuter/Downloads/HCC/submit/simulations_code/clustfns.R")
  source("/Users/polinasuter/Downloads/HCC/submit/simulations_code/comparemodels.R")
  #

  res<-list()
  res$ROC<-NULL
  res$accuracy<-NULL
  if(!is.null(basename) & !is.null(path)) saveRDS(res,paste(path,basename,".rds",sep=""))
    sseed<-randomseed+rep
    set.seed(sseed)
    dflocal<-NULL
    #generate mixture
    bnmixt<-genMixture(k=k, type="mixed",centersignal=centersignal,
                       sigma0=sigma0, ssvec=ssvec, n=n, avpar=avpar, deltamu=deltamu, lB=lB, uB=uB,
                       shdpct=shdpct, randomseed=sseed, eqval=eqval,
                       mixedpar=mixedpar, signalroots=TRUE)
    metainfo<-list(rep=rep,seed=sseed,N=sum(bnmixt$ssvec)/k,dmu=deltamu,es=sigma0,k=k, r2=mean(bnmixt$r2))
    if(savedata) {
      res<-bnmixt
    } else {
      res$info<-bnmixt$info
      if(type=="mixed") {
        res$info$nbin<-bnmixt$nbin
      }
      if(edgep) {
        edgepmat<-simedgepmat2(bnmixt,pf=2,FNrate=FNrate,FPrate=FPrate)
      } else {
        edgepmat=NULL
      }
      if(!onlyother) {
        #learn mixture
        set.seed(sseed)
        databn<-makeOmicsObject(bnmixt,n,mixedpar$nbin)
        omicsobj<-bnInfo(databn,types=c("b","c"),omics=c("M","T"))
        if(algorithm=="mcmcMAP") {
          bnres<-bnclustOmics::bnclustOmics(databn,omicsobj, blacklist=NULL, edgepmat=edgepmat,kclust=k,
                                            maxEM=maxEM,startpoint = startpoint,baseprob=3/(k+2),plus1it=plus1it,epmatrix=ROC)

        }
        relab<-checkmembership(k,bnmixt$membership,bnres$memb)$relabel
        bnres$memb<-relabMembership(bnres$memb, relab)
        bnres$lambdas<-relabLambdas(bnres$lambdas,relab)
        bnres$DAGs<-relabDAGs(bnres$DAGs,relab)
        if(!is.null(bnres$ep)) {
          bnres$ep<-relabDAGs(bnres$ep,relab)
        }

        if(algorithm=="mcmcMAP") res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="MAP",rep=rep,seed=sseed)) else
          res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="cons",rep=rep,seed=sseed))
        #get network comparisons
        if(ROC) {
          res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="p",p=c(0.3,0.5,0.7,0.9,0.95,0.99)))
        }

        dflocal<-rbind(dflocal,data.frame(clustaccuracy(bnmixt$membership,bnres$memb,k,ss=nrow(bnmixt$data)),algorithm=algorithm, metainfo))
      }

      res$accuracy<-rbind(res$accuracy,dflocal)
      if(!is.null(basename) & !is.null(path)) saveRDS(res,paste(path,basename,".rds",sep=""))
    }
  return(res)
}

#

#this code can be used to run 50 replicates of simulations in parallel
library(parallel)
rep<-c(1:2)
cl <- makeCluster(3)
outputClApply <- parallel::clusterApply(cl, rep, simBNclustcore,
                                        type="mixed",
                                        n=30,
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
                                        onlyother=FALSE,edgep=TRUE,savedata=FALSE,
                                        path=NULL,base=NULL,FPrate=FPrate,FNrate=FNrate)
stopCluster(cl)

#ADD LINE TO SAVE THE RESULT!
saveRDS(outputClApply,paste(path,base,"May22","FP",args[1],"FN",args[2],".rds",sep=""))


# out11<-outputClApply
# out51<-outputClApply #FP FN
# out15<-outputClApply
#
# out11[[1]]$ROC[1:4,c(4,8)]
# out51[[1]]$ROC[1:4,c(4,8)]
# out15[[1]]$ROC[1:4,c(4,8)]
#
# out11[[2]]$ROC[1:4,c(4,8)]
# out51[[2]]$ROC[1:4,c(4,8)]
# out15[[2]]$ROC[1:4,c(4,8)]
#
#
