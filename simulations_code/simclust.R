#make sure that seed fixes structure and data
simBNclust<-function(nrep,type, n=50, k=3, avpar=1, ssvec=rep(100,k), centersignal="strong",
                     sigma0=0.3, deltamu=20,lB=0.5, uB=1.5, shdpct=20, eqval=0,
                     randomseed=100,mixedpar=list(nbin,avchildren,par1=0.1,par2=7),
                     algorithm=c("bnclust","allother"),hardlimit=11,maxEM=6,plus1it=4,p=0.5,
                     startpoint=c("random","mclust","mclustPCA"),ROC=TRUE,addalgos=TRUE,basename=NULL,path=NULL,
                     onlyother=FALSE,edgep=FALSE,savedata=FALSE,recimlr=FALSE,errRate=0) {
  res<-list()
  res$ROC<-NULL
  res$accuracy<-NULL
  if(!is.null(basename) & !is.null(path)) saveRDS(res,paste(path,basename,".rds",sep=""))
  for(i in 1:nrep) {
    print("Rep")
    print(i)
    sseed<-randomseed+i
    set.seed(sseed)
    dflocal<-NULL
    #generate mixture
    bnmixt<-genMixture(k=k, type="mixed",centersignal=centersignal,
                     sigma0=sigma0, ssvec=ssvec, n=n, avpar=avpar, deltamu=deltamu, lB=lB, uB=uB,
                     shdpct=shdpct, randomseed=sseed, eqval=eqval,
                     mixedpar=mixedpar, signalroots=TRUE)
    metainfo<-list(rep=i,seed=sseed,N=sum(bnmixt$ssvec)/k,dmu=deltamu,es=sigma0,k=k, r2=mean(bnmixt$r2))
    if(savedata) {
      res[[i]]<-bnmixt
      } else {
    res$info<-bnmixt$info
    if(type=="mixed") {
      res$info$nbin<-bnmixt$nbin
    }
    if(edgep) {
      edgepmat<-simedgepmat(bnmixt,pf=2,errRate=errRate)
    } else {
      edgepmat=NULL
    }
    if(!onlyother) {
    #learn mixture
    set.seed(sseed)
    databn<-makeOmicsObject(bnmixt,n,mixedpar$nbin)
    omicsobj<-bnInfo(databn,types=c("b","c"),omics=c("M","T"))
    if(k==3) baseprob<-0.4 else baseprob<-3/(k+2)
    if(algorithm=="bnclust") {
      bnres<-bnClustOmics::bnclustOmics(databn,omicsobj, blacklist=NULL, edgepmat=edgepmat,kclust=k, hardlim=hardlimit,chixi=0.5,
                                        maxEM=maxEM,startpoint = startpoint,baseprob=baseprob,plus1it=plus1it,epmatrix=ROC)
    }

    relab<-checkmembership(k,bnmixt$membership,bnres$memb)$relabel
    bnres$memb<-relabMembership(bnres$memb, relab) #last was changeto
    bnres$lambdas<-relabLambdas(bnres$lambdas,relab) #last was changeto
    bnres$DAGs<-relabDAGs(bnres$DAGs,relab)
    if(!is.null(bnres$ep)) {
      bnres$ep<-relabDAGs(bnres$ep,relab)
    }

    if(algorithm=="bnclust") res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="MAP",rep=i,seed=sseed))
    else res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="cons",rep=i,seed=sseed))
    #get network comparisons
    if(ROC) {
      res$ROC<-rbind(res$ROC,compareMixt(bnres,bnmixt,dag="p",p=c(0.3,0.5,0.7,0.9,0.95,0.99)))
    }

    #save the result
    #acc<-bnres$corrvec[length(bnres$corrvec)]/sum(ssvec)
    dflocal<-rbind(dflocal,data.frame(clustaccuracy(bnmixt$membership,bnres$memb,k,ss=nrow(bnmixt$data)),algorithm=algorithm, metainfo))
    }
    if(algorithm!="bnclust") {
      print("mclust")
      dflocal<-rbind(dflocal,data.frame(accmclust(bnmixt,PCA=TRUE, npca=k+2),algorithm="mclustPCA", metainfo))
      print("hclust")
      dflocal<-rbind(dflocal,data.frame(acchclust(bnmixt,PCA=TRUE, npca=k+2),algorithm="hclustPCA", metainfo))
      print("kmeans")
      dflocal<-rbind(dflocal,data.frame(acckmeans(bnmixt,PCA=TRUE, npca=k+2),algorithm="kmeansPCA", metainfo))
      #print("iclust")
      if((sum(ssvec)/k)<30) dflocal<-rbind(dflocal,data.frame(acciclust(bnmixt),algorithm="iclust", metainfo))
      dflocal<-rbind(dflocal,data.frame(accMOFA(bnmixt),algorithm="MOFA", metainfo))
      dflocal<-rbind(dflocal,data.frame(accCIMLRco(bnmixt),algorithm="CIMLRco", metainfo))
      dflocal<-rbind(dflocal,data.frame(accCIMLR(bnmixt),algorithm="CIMLR", metainfo))
    }
    res$accuracy<-rbind(res$accuracy,dflocal)
    if(!is.null(basename) & !is.null(path)) saveRDS(res,paste(path,basename,".rds",sep=""))
    }

  }
  return(res)
}



