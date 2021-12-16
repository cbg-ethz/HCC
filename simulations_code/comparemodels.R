compareMixt<-function(bnres,bnmixt,dag=c("MAP","cons","p"),p=c(0.5,0.7,0.9),rep=1,seed=100) {
if(bnmixt$info$type=="mixed") {
  truecpdags<-lapply(bnmixt$DAGs,mixed2cp,bnmixt$nbin,bnmixt$info$n)
} else {
  truecpdags<-lapply(bnmixt$DAGs,dag2cpdag)
}
  compdf<-NULL
  if(dag=="MAP") {
    if(bnmixt$info$type=="mixed") {
      estcpdags<-lapply(bnres$DAGs,mixed2cp,bnmixt$nbin,bnmixt$info$n)
      for(i in 1:length(estcpdags)) {
        compdf<-rbind(compdf,compareDAGs(estcpdags[[i]],truecpdags[[i]],cpdag=FALSE))
      }
    } else {
      for(i in 1:length(truecpdags)) {
        compdf<-rbind(compdf,compareDAGs(bnres$DAGs[[i]],truecpdags[[i]],cpdag=FALSE))
      }
    }
    compdf<-cbind(compdf,bnmixt$ssvec)
    compdf<-cbind(compdf,-1)
    compdf<-cbind(compdf,c(1:bnmixt$info$k))
  } else if (dag=="cons") {
    for(i in 1:length(truecpdags)) {
      compdf<-rbind(compdf,compareDAGs(bnres$consmodel[[i]],truecpdags[[i]],cpdag=FALSE))
    }
    compdf<-cbind(compdf,bnmixt$ssvec)
    compdf<-cbind(compdf,bnres$p)
    compdf<-cbind(compdf,c(1:bnmixt$info$k))
  } else {
    for(i in 1:length(truecpdags)) {
      models<-modellist(bnres$ep[[i]],p)
      compres<-Reduce(rbind,lapply(models,compareDAGs,truegraph=truecpdags[[i]]))
      compdf<-rbind(compdf,compres)
    }
    lp<-length(p)
    compdf<-cbind(compdf,rep(bnmixt$ssvec,lp))
    compdf<-cbind(compdf,rep(p,bnmixt$info$k))
    compdf<-cbind(compdf,as.vector(sapply(1:bnmixt$info$k,rep,lp)))
  }
  compdf<-cbind(as.data.frame(compdf),dag)
  compdf<-cbind(compdf,bnres$algorithm)
  compdf<-cbind(compdf,rep)
  compdf<-cbind(compdf,seed)
  colnames(compdf)[9:15]<-c("N","p","cluster","model","algorithm","rep","seed")
  rownames(compdf)<-c(1:nrow(compdf))
  return(compdf)
}



