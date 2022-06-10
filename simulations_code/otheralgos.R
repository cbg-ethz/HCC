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
clustaccuracy<-function(truememb,estmemb,k,ss,abs=FALSE,prec=FALSE) {
  if(abs) {
    return((checkmembership(k,truememb,estmemb)$ncorr)/ss)
  } else {
    if(prec) {
      return(data.frame(precision=precision_clusters(estmemb,truememb), ARI=adjustedRandIndex(truememb, estmemb)))
    } else {
      return(data.frame(ABS=(checkmembership(k,truememb,estmemb)$ncorr)/ss, ARI=adjustedRandIndex(truememb, estmemb)))
    }
  }
}
precision_clusters<-function(clusters,true_clusters) {
  cluster_labels<-unique(clusters)
  n_clust<-length(cluster_labels)
  N<-length(clusters)
  
  cluster_index<-lapply(cluster_labels,single_clusters,clusters)
  cluster_bins<-lapply(cluster_index,make_bins,true_clusters)
  
  tot_pairs<-N*(N-1)/2
  tot_pos<-sum(unlist(lapply(cluster_bins,total_pos_pairs)))
  tot_neg<-total_neg_pairs(tot_pairs,tot_pos)
  
  TP<-sum(unlist(lapply(cluster_bins,pairs_TP)))
  FP<-tot_pos-TP
  FN<-pairs_FP(cluster_bins,cluster_labels)
  TN<-tot_neg-FP
  
  Pr<-TP/(TP+FP)
  Rec<-TP/(TP+FN)
  F1<-2*(Pr*Rec)/(Pr+Rec)
  return(c(Pr,Rec,F1))
}
single_clusters<-function(i,clusters) {
  return(which(clusters==i))
}
make_bins<-function(index,true_clusters){
  return(true_clusters[index])
}
total_pos_pairs<-function(cluster_bin){
  return(choose(length(cluster_bin),2))
}
pairs_TP<-function(cluster_bin){
  taby<-table(cluster_bin)
  return(sum(sapply(taby[which(taby>1)],choose,2)))
}
pairs_FP<-function(cluster_bins,cluster_labels){
  mm_matrix<-matrix(nrow=length(cluster_labels),
                    ncol=length(cluster_bins))
  for(i in 1:length(cluster_bins)) {
    for(j in cluster_labels) {
      mm_matrix[i,j]<-length(which(cluster_bins[[i]]==j))
    }
  }
  
  mm_tot<-0
  n_row<-nrow(mm_matrix)
  
  for(i in 1:ncol(mm_matrix)) {
    for(j in 1:(nrow(mm_matrix)-1)) {
      mm_tot<-mm_tot+mm_matrix[j,i]*sum(mm_matrix[(j+1):n_row,i])
    }
  }
  return(mm_tot)
}
total_neg_pairs<-function(tot_pairs,tot_pos) {
  return(tot_pairs-tot_pos)
}





