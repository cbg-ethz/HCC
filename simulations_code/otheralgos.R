plotPCA<-function(BNmixt,colk=NULL,...) {

  var0<-which(apply(BNmixt$data,2,sd)==0)
  if(length(var0>0)) PCAdata<-BNmixt$data[,-var0] else PCAdata<-BNmixt$data
  pca_res <- prcomp(PCAdata, scale. = TRUE)
  pchn<-c(15,16,17)
  if(is.null(colk)) {
    plot(pca_res$x[,1],pca_res$x[,2],col=factor(BNmixt$membership),pch=pchn[BNmixt$membership],
         xlab="PC1",ylab="PC2", ...)
  } else {

    plot(pca_res$x[,1],pca_res$x[,2],col=col3[BNmixt$membership],pch=pchn[BNmixt$membership],
         xlab="PC1",ylab="PC2", ...)
  }
}
accSIMLR<-function(BNmixt,abs=FALSE) {
  k<-length(BNmixt$DAGs)
  ss<-nrow(BNmixt$data)
  simlrfit<-SIMLR(t(BNmixt$data),c=k,cores.ratio = 0)
  return(clustaccuracy(BNmixt$membership,simlrfit$y$cluster,k,ss,abs=abs))
}
accCIMLR<-function(BNmixt,abs=FALSE) {
  k<-length(BNmixt$DAGs)
  ss<-nrow(BNmixt$data)
  nbin<-BNmixt$nbin
  ncont<-BNmixt$info$n
  datalist<-list()
  datalist[[1]]<-t(BNmixt$data[,1:nbin])
  datalist[[2]]<-t(BNmixt$data[,1:ncont+nbin])
  cimlrfit<-CIMLR::CIMLR(datalist,c=k,cores.ratio = 0)
  return(clustaccuracy(BNmixt$membership,cimlrfit$y$cluster,k,ss,abs=abs))
}
accCIMLRco<-function(BNmixt,abs=FALSE) {
  k<-length(BNmixt$DAGs)
  ss<-nrow(BNmixt$data)
  nbin<-BNmixt$nbin
  ncont<-BNmixt$info$n
  datalist<-list()
  datalist[[1]]<-t(BNmixt$data[,1:ncont+nbin])#
  cimlrfit<-CIMLR::CIMLR(datalist,c=k,cores.ratio = 0)
  return(clustaccuracy(BNmixt$membership,cimlrfit$y$cluster,k,ss,abs=abs))
}
accmclust<-function(BNmixt,PCA=FALSE,npca=5,abs=FALSE) {
  k<-length(BNmixt$DAGs)
  ss<-nrow(BNmixt$data)
  if(!PCA) {
    mclustfit<-Mclust(BNmixt$data,G=k)
    return(clustaccuracy(BNmixt$membership,mclustfit$classification,k,ss,abs=abs))
  } else {
    var0<-which(apply(BNmixt$data,2,sd)==0)
    if(length(var0>0)) PCAdata<-BNmixt$data[,-var0] else PCAdata<-BNmixt$data
    pca_res <- prcomp(PCAdata, scale. = TRUE)
    mclustfit<-Mclust(pca_res$x[,1:npca],G=k)
    return(clustaccuracy(BNmixt$membership,mclustfit$classification,k,ss,abs=abs))
  }
}
acchclust<-function(BNmixt,abs=FALSE,npca=5,PCA=FALSE) {
  k<-length(BNmixt$DAGs)
  ss<-nrow(BNmixt$data)
  if(!PCA) {
  dist_mat <- dist(BNmixt$data, method = 'euclidean')
  hclust_avg <- hclust(dist_mat, method = 'ward.D2')
  cut_avg <- cutree(hclust_avg, k = k)

  } else {
    var0<-which(apply(BNmixt$data,2,sd)==0)
    if(length(var0>0)) PCAdata<-BNmixt$data[,-var0] else PCAdata<-BNmixt$data
    pca_res <- prcomp(PCAdata, scale. = TRUE)
    dist_mat <- dist( pca_res$x[,1:npca], method = 'euclidean')
    hclust_avg <- hclust(dist_mat, method = 'ward.D2')
    cut_avg <- cutree(hclust_avg, k = k)
  }
  return(clustaccuracy(BNmixt$membership,cut_avg,k,ss,abs=abs))
}
acckmeans<-function(BNmixt,abs=FALSE,npca=5,PCA=FALSE) {
  k<-length(BNmixt$DAGs)
  ss<-nrow(BNmixt$data)
  if(!PCA) {
  kmeansfit <- kmeans(BNmixt$data, k)
  } else {
    var0<-which(apply(BNmixt$data,2,sd)==0)
    if(length(var0>0)) PCAdata<-BNmixt$data[,-var0] else PCAdata<-BNmixt$data
    pca_res <- prcomp(PCAdata, scale. = TRUE)
    kmeansfit <- kmeans( pca_res$x[,1:npca], k)
  }
  return(clustaccuracy(BNmixt$membership,kmeansfit$cluster,k,ss,abs=abs))
}
acciclust<-function(BNmixt,abs=FALSE) {
  k<-length(BNmixt$DAGs)
  dd<-dividedata(BNmixt)
  ss<-nrow(BNmixt$data)
  var0<-which(apply(dd$bin,2,sd)==0)
  if(length(var0>0)) dd$bin<-dd$bin[,-var0] else dd$bin<-dd$bin
  iclustfit<-iClusterPlus(dt1=dd$bin,dt2=dd$cont,type=c("binomial","gaussian"),K=k-1, maxiter=10)
  #print(clustaccuracy(BNmixt$membership,iclustfit$clusters,k,ss,abs=TRUE))
  #print(clustaccuracy(BNmixt$membership,iclustfit$clusters,k,ss,abs=FALSE))
  return(clustaccuracy(BNmixt$membership,iclustfit$clusters,k,ss,abs=abs))
}
accMOFA<-function(BNmixt,abs=FALSE,accuracy=TRUE) {
  k<-length(BNmixt$DAGs)
  dd<-dividedata(BNmixt)
  ss<-nrow(BNmixt$data)
  rownames(dd$bin)<-rownames(dd$cont)<-paste("S",1:ss,sep="")
  var0<-which(apply(dd$bin,2,sd)==0)
  if(length(var0>0)) dd$bin<-dd$bin[,-var0] else dd$bin<-dd$bin

  HCCDI<-list()
  HCCDI[["M"]]<-t(dd$bin)
  HCCDI[["T"]]<-t(dd$cont)
  MOFAobject <- createMOFAobject(HCCDI)

  mae_HCC <- MultiAssayExperiment(
    experiments = HCCDI
  )
  MOFAobject <- createMOFAobject(mae_HCC)
  MOFAobject

  DataOptions <- getDefaultDataOptions()
  DataOptions

  ModelOptions <- getDefaultModelOptions(MOFAobject)
  ModelOptions$numFactors <- 6
  ModelOptions

  TrainOptions <- getDefaultTrainOptions()
  TrainOptions$DropFactorThreshold <- 0.01
  TrainOptions$seed <- 200
  TrainOptions

  MOFAobject <- prepareMOFA(
    MOFAobject,
    DataOptions = DataOptions,
    ModelOptions = ModelOptions,
    TrainOptions = TrainOptions
  )

  MOFAobject <- runMOFA(MOFAobject)
  if(accuracy) {
  clusters <- clusterSamples(MOFAobject, k=k, factors=1:MOFAobject@Dimensions$K)
  #print(clustaccuracy(BNmixt$membership,iclustfit$clusters,k,ss,abs=TRUE))
  #print(clustaccuracy(BNmixt$membership,iclustfit$clusters,k,ss,abs=FALSE))
  return(clustaccuracy(BNmixt$membership,clusters,k,ss,abs=abs))
  } else  {
    return(MOFAobject)
  }
}
clustaccuracy<-function(truememb,estmemb,k,ss,abs=FALSE) {
  if(abs) {
    return((checkmembership(k,truememb,estmemb)$ncorr)/ss)
  } else {
    return(data.frame(ABS=(checkmembership(k,truememb,estmemb)$ncorr)/ss, ARI=adjustedRandIndex(truememb, estmemb)))
  }
}





