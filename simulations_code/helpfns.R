getModels<-function(ep, thresholds) {
  DAGlist<-list()
  n<-ncol(ep)
  labels<-colnames(ep)
  res<-list()
  for(i in 1:length(thresholds)) {
    incidence <- matrix(rep(0, n * n), nrow = n, ncol = n)
    colnames(incidence)<-rownames(incidence)<-labels
    incidence[which(ep > thresholds[i])] <- 1
    res[[i]]<-m2graph(incidence)
  }
  return(res)
}
simedgepmat2<-function(BNmixt,pf=2,blacklist=FALSE,FNrate=0,FPrate=0){
  edgepmat<-1*Reduce('|',BNmixt$DAGs)
  if(FNrate>0) {
    ones<-which(edgepmat==1)
    lo<-length(ones)
    totsamp<-ceiling(FNrate*lo)
    edgepmat[ones[sample.int(lo,ceiling(totsamp/2))]]<-0
  }

  if(FPrate>0) {
    zeros<-which(edgepmat==0)
    ones<-which(edgepmat==1)
    lo<-length(ones)
    lz<-length(zeros)
    totsamp<-ceiling(FPrate*lo)
    edgepmat[zeros[sample.int(lz,ceiling(totsamp/2))]]<-1
  }

  if(blacklist) {
    return(edgepmat)
  } else {
    edgepmat<-(pf-1)*(!edgepmat)
    return(edgepmat+1)
  }
}
simedgepmat<-function(BNmixt,pf=2,blacklist=FALSE,errRate=0){
  edgepmat<-1*Reduce('|',BNmixt$DAGs)
  if(errRate>0) {
    ones<-which(edgepmat==1)
    zeros<-which(edgepmat==0)
    lo<-length(ones)
    lz<-length(zeros)
    totsamp<-ceiling(errRate*lo)
    edgepmat[ones[sample.int(lo,ceiling(totsamp/2))]]<-0
    edgepmat[zeros[sample.int(lz,ceiling(totsamp/2))]]<-1
  }
  if(blacklist) {
    return(edgepmat)
  } else {
    edgepmat<-(pf-1)*(!edgepmat)
    return(edgepmat+1)
  }
}
dividedata<-function(BNmixt) {
  data<-list()
  data$bin<-BNmixt$data[,1:BNmixt$nbin]
  data$cont<-BNmixt$data[,1:BNmixt$info$n+BNmixt$nbin]
  return(data)
}
FSbyVarmy<-function(dataf,value) {
  myvars <- apply(dataf,1, var,na.rm=TRUE)
  myvars <- sort(myvars,decreasing=TRUE)
  myvars <- myvars[1:value]
  return(dataf[names(myvars),])
}
gettopWeightsMOFA<-function(MOFAobject,fact=c(1:5),perc=100*rep(1/length(fact),length(fact))) {
  topfeat<-c()
  k<-1
  for(i in fact) {
    MOFAweights <- getWeights(
      MOFAobject,
      views = "all",
      factors = i,
      as.data.frame = TRUE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
    )
    MOFAweights<-cbind(MOFAweights,abs(MOFAweights$value))
    colnames(MOFAweights)[5]<-"abs"
    MOFAorder<-order(MOFAweights$abs,decreasing =TRUE)
    MOFAweights<-MOFAweights[MOFAorder,]
    topfeat<-unique(c(topfeat,MOFAweights$feature[1:perc[k]]))
    k<-k+1
  }
  return(topfeat)
}
getTopFeats<-function(MOFAobject,views,facts,ntop) {
  MOFAweights <- getWeights(
      MOFAobject,
      views = views,
      factors = facts,
      as.data.frame = FALSE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
    )
  MOFAweights<-lapply(MOFAweights, abs)
  MOFAweights<-lapply(MOFAweights,rowSums)
  allfeats<-c()
  for(i in 1:length(MOFAweights)) {
    allfeats<-c(allfeats,unlist(MOFAweights[[i]]))
  }
  feats<-names(sort(allfeats,decreasing=TRUE)[1:ntop])
  print(feats)
  return(feats)
}
getTopFeats2<-function(MOFAobject,views,facts,ntop) {
  MOFAweights <- get_weights(
    MOFAobject,
    views = views,
    factors = facts,
    as.data.frame = FALSE    # if TRUE, it outputs a long dataframe format. If FALSE, it outputs a wide matrix format
  )
  MOFAweights<-lapply(MOFAweights, abs)
  MOFAweights<-lapply(MOFAweights,rowSums)
  allfeats<-c()
  for(i in 1:length(MOFAweights)) {
    allfeats<-c(allfeats,unlist(MOFAweights[[i]]))
  }
  feats<-names(sort(allfeats,decreasing=TRUE)[1:ntop])
  print(feats)
  return(feats)
}

is.error <- function(x) inherits(x, "try-error")
topNodes<-function(mixttest,topn=100,plus=FALSE){
  sds<-apply(mixttest$data,2,sd)
  means<-apply(abs(mixttest$data),2,mean)
  neib<-means/sds
  ordered<-setdiff(order(neib,decreasing = TRUE),c(1:mixttest$nbin))
  if(plus) {
    commonDAG<-1*Reduce('|',mixttest$DAGs)
    npar<-apply(commonDAG,2,sum)
    nch<-apply(commonDAG,1,sum)
    neib<-(npar+nch)[ordered]
    ordered<-ordered[which(neib>0)]
  }
    return(ordered[1:topn])
}
hybridNodes<-function(mixttest,nodesMOFA=NULL,topn=100,plus=FALSE) {
  if(is.null(nodesMOFA)) {
    varint<-paste("V",topNodes(mixttest,topn=topn,plus=plus),sep="")
    return(c(varint))
  } else {
  varint<-paste("V",topNodes(mixttest,topn=ceiling(topn/2),plus=plus),sep="")
  mofasel<-setdiff(nodesMOFA,varint)[1:floor(topn/2)]
  return(c(varint,mofasel))
  }
}
clustSubset<-function(mixttest,nodesB,nodesC,commonDAG,k) {
  newmixt<-mixttest
  newmixt$data<-mixttest$data[,c(nodesB,nodesC)]
  nint<-sum(commonDAG[c(nodesB,nodesC),c(nodesB,nodesC)])
  newmixt$info$n<-length(nodesC)
  mclustfit<-mclustPCA(newmixt$data, k, k+2)
  bnfit<-bnclust(newmixt$data,kclust=k,type="mixed",compare=TRUE,truememb=newmixt$membership,MAP=TRUE,
                 startmemb=mclustfit,
                 plus1it=3, bgnodes=1:mixttest$nbin,maxEM=4,nbin=mixttest$nbin)
  acc<-clustaccuracy(newmixt$membership,bnfit$memb,k=k,ss=nrow(newmixt$data),abs=FALSE)$ARI
  return(list(nint,acc))
}
mclustPCA<-function(datas,k,npca=k+2) {
      var0<-which(apply(datas,2,sd)==0)
      if(length(var0)>0) pcadata<-datas[,-var0] else pcadata<-datas
      pca_res <- prcomp(pcadata, scale. = TRUE)
      startmemb<-Mclust(pca_res$x[,1:(npca)],G=k)$classification

return(startmemb)
}
clustSubset2<-function(mixttest,nodesB,nodesC,commonDAG,k) {
  newmixt<-mixttest
  newmixt$data<-list()
  newmixt$data[["M"]]<-mixttest$data[,nodesB]
  newmixt$data[["T"]]<-mixttest$data[,nodesC]
  nint<-sum(commonDAG[c(nodesB,nodesC),c(nodesB,nodesC)])
  newmixt$info$n<-length(nodesC)
  omicsobj<-bnInfo(newmixt$data,types=c("b","c"),omics=c("M","T"))
  bnfit<-bnclustOmics::bnclustOmics(newmixt$data,omicsobj, blacklist=NULL, edgepmat=NULL,kclust=k,
                                    maxEM=4,startpoint = "mclustPCA",baseprob=0.4,plus1it=4,epmatrix=FALSE)
  acc<-clustaccuracy(newmixt$membership,bnfit$memb,k=k,ss=nrow(newmixt$data[[1]]),abs=FALSE)$ARI
  return(list(nint,acc))
}
makeOmicsObject<-function(mixttest,n,nbin){
  datanew<-list()
  colnames(mixttest$data)<-paste("V",1:(n+nbin),sep="")
  rownames(mixttest$data)<-paste("S",1:(nrow(mixttest$data)),sep="")
  nodesB<-paste("V",1:nbin,sep="")
  nodesC<-paste("V",1:n+nbin,sep="")
  datanew[["M"]]<-mixttest$data[,nodesB]
  datanew[["T"]]<-mixttest$data[,nodesC]
  return(datanew)
}
selectFeatures<-function(mixttest,nodesB,nodesC){
  newmixt<-mixttest
  newmixt$data<-list()
  newmixt$data[["M"]]<-mixttest$data[,nodesB]
  newmixt$data[["T"]]<-mixttest$data[,nodesC]
  newmixt$data<-cbind(newmixt$data[["M"]],newmixt$data[["T"]])
  newmixt$info$n<-length(nodesC)
  return(newmixt)
}

