library(mclust)
library(CIMLR)
library(iClusterPlus)
library(MOFA2)
library(factoextra)
#clinical data for fitting Cox model
HCC_surv_df<-readRDS("/Users/polinasuter/Downloads/HCC/submit/HCCinputs/info.rds") #load clinical data
source("/Users/polinasuter/Downloads/HCC/submit/helpresfns.R")

clustres<-matrix(ncol=50,nrow=10)
rownames(clustres)<-c("mclust","kmeans_3","kmeans_2",
                      "hclust_3","hclust_2","CIMLR",
                      "CIMLR_co","iclust","MOFA")
FSbyVarmy<-function(dataf,value) {
  myvars <- apply(dataf,1, var,na.rm=TRUE)
  myvars <- sort(myvars,decreasing=TRUE)
  myvars <- myvars[1:value]
  return(dataf[names(myvars),])
}
HCCfullz<-readRDS("/Users/polinasuter/Downloads/HCC/submit/HCCinputs/HCCfullz.rds")
dim(HCCfullz$M)
dim(HCCfullz$CN)
dim(HCCfullz$T)
dim(HCCfullz$P)
dim(HCCfullz$PP)

HCCcontz<-Reduce('rbind',HCCfullz[c("T","P","PP")])
HCCcontz<-t(HCCcontz)
HccMCN<-Reduce('rbind',HCCfullz[c("M","CN")])
HccfullstackedPCAz<-cbind(t(HccMCN),HCCcontz)
dim(HccfullstackedPCAz)

#need PCA for general clustering methods
npca<-5
k<-3
var0<-which(apply(HccfullstackedPCAz,2,sd)==0)
if(length(var0>0)) HccfullstackedPCAz<-HccfullstackedPCAz[,-var0] else HccfullstackedPCAz<-HccfullstackedPCAz
pca_res <- prcomp(HccfullstackedPCAz, scale. = TRUE)

########
#mclust#
########

set.seed(100)
BIC <- mclustBIC(pca_res$x[,1:npca])
mclustfit <- Mclust(pca_res$x[,1:npca], x = BIC)
#optimal number of clusters is based on BIC
membfull <- mclustfit$classification
length(unique(mclustfit$classification))
nonadjustedCox(membfull,HCC_surv_df)
adjustedCox(membfull,HCC_surv_df)

########
#hclust#
########

set.seed(100)
dist_mat <- dist( pca_res$x[,1:npca], method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'ward.D2')
membfull <- cutree(hclust_avg, k = k)

#optimal number of clusters 2, wss method
fviz_nbclust(pca_res$x[,1:npca], hcut, method = "wss") +
  labs(subtitle = "wss")
#perform clustering for k=6
membfull <- cutree(hclust_avg, k = 6)
nonadjustedCox(membfull,HCC_surv_df)
adjustedCox(membfull,HCC_surv_df)

#########
#k-means#
#########

library(factoextra)
library(NbClust)
set.seed(100)
kmeansfit<-kmeans( pca_res$x[,1:npca], k)
membfull<-kmeansfit$cluster

#optimal number of clusters 2, ww method
fviz_nbclust(pca_res$x[,1:npca], kmeans, method = "wss") +
  labs(subtitle = "wss")

#perform clustering for k=6
kmeansfit<-kmeans(pca_res$x[,1:npca], 6)
membfull<-kmeansfit$cluster
nonadjustedCox(membfull,HCC_surv_df)
adjustedCox(membfull,HCC_surv_df)

#################
#CIMLR all omics#
#################

#optimal number of clusters is 3
numk<-2:7
optk<-CIMLR::CIMLR_Estimate_Number_of_Clusters(HCCfullz,numk)
which.min(optk$K1)
which.min(optk$K2)

#perform clustering for k=3
set.seed(100)
cimplrfit<-CIMLR::CIMLR(HCCfullz,c=k,cores.ratio = 0)
membfull<-cimplrfit$y$cluster
nonadjustedCox(membfull,HCC_surv_df)
adjustedCox(membfull,HCC_surv_df)

#################
#CIMLR cont only#
#################

#optimal number of clusters is 3
numk<-2:7
optk<-CIMLR::CIMLR_Estimate_Number_of_Clusters(HCCfullz[c("T","P","PP")],numk)
which.min(optk$K1)
which.min(optk$K2)

#perform clustering for k=3
set.seed(100)
cimplrfit<-CIMLR::CIMLR(HCCfullz[c("T","P","PP")],c=k,cores.ratio = 0) #
membfull<-cimplrfit$y$cluster
nonadjustedCox(membfull,HCC_surv_df)
adjustedCox(membfull,HCC_surv_df)

######
#MOFA#
######

 library('reticulate')

 # Run these once before running the batch jobs to prepare the python environment
 conda_create("r-reticulate")
 use_condaenv("r-reticulate", required = TRUE)
 py_install('mofapy2', pip_options = '--user', pip = TRUE)

topn<-c(500,2000,2000,1000,1000)
 for(i in 1:5){
  HCCfullz[[i]]<-as.matrix(HCCfullz[[i]])
  colnames(HCCfullz[[i]])<-colnames(HCCfullz[[1]])
  #select features with higher variance as recommended by the authors
  HCCfullz[[i]]<-FSbyVarmy(HCCfullz[[i]],topn[i])
 }

 MOFAobject <- create_mofa(HCCfullz)
 data_opts <- get_default_data_options(MOFAobject)
 model_opts <- get_default_model_options(MOFAobject)
 model_opts$num_factors<-6
 train_opts <- get_default_training_options(MOFAobject)
 train_opts$drop_factor_threshold<-0.01
 MOFAobject <- prepare_mofa(
   object = MOFAobject,
   data_options = data_opts,
   model_options = model_opts,
   training_options = train_opts
 )
 outfile = file.path("/Users/polinasuter/Downloads/MOFAresults/","model1.hdf5")
 model <- run_mofa(MOFAobject, outfile)
 mofaclusters <- cluster_samples(model, k=3, factors="all")

 nonadjustedCox(mofaclusters$cluster,HCC_surv_df)
 adjustedCox(mofaclusters$cluster,HCC_surv_df)

 set.seed(100)
 wss<-vector()
 for(i in 1:7) {
   clusters <- cluster_samples(model, k=i, factors="all")
   wss[i]<-sum(clusters$withinss)
 }
 #elbow method
 plot(c(1:7),wss,type="b",col="blue")
#optimal number of clusters is 4
 mofaclusters <- cluster_samples(model, k=4, factors="all")
 nonadjustedCox(mofaclusters$cluster,HCC_surv_df)
 adjustedCox(mofaclusters$cluster,HCC_surv_df)

 ##############
 #iClusterPlus#
 ##############
 #We use a reduced set of features, same as for MOFA
 HCCcontz<-Reduce('rbind',HCCfullz[c("T","P","PP")])
 HCCcontz<-t(HCCcontz)

 set.seed(100)
 iclustfit<- iClusterPlus(dt1=t(HCCfullz$M),dt2=t(HCCfullz$CN),dt3=HCCcontz,
                          type=c("binomial","gaussian","gaussian"),K=k-1, maxiter=20)
 membfull<-iclustfit$clusters
 nonadjustedCox(membfull,HCC_surv_df)
 adjustedCox(membfull,HCC_surv_df)

 #identifying the optimal number of clusters, runtime >20hrs
 for(k in 1:5) {
   cv.fit = tune.iClusterPlus(cpus=2,dt1=t(HCCfullz$M),dt2=t(HCCfullz$CN),dt3=HCCcontz,
                              type=c("binomial","gaussian","gaussian"),K=k,maxiter=20)
   save(cv.fit, file=paste("cv.fit.k",k,".Rdata",sep=""))
 }
 output=alist()
 files=grep("cv.fit",dir())
 for(i in 1:length(files)){
   load(dir()[files[i]])
   output[[i]]=cv.fit
 }
 nLambda = nrow(output[[1]]$lambda)
 nK = length(output)
 BIC = getBIC(output)
 devR = getDevR(output)

 minBICid = apply(BIC,2,which.min)
 devRatMinBIC = rep(NA,nK)
 for(i in 1:nK){
   devRatMinBIC[i] = devR[minBICid[i],i]
 }

 clusters=getClusters(output)
 rownames(clusters)=rownames(gbm.exp)
 colnames(clusters)=paste("K=",2:(length(output)+1),sep="")
 best.fit=output[[k]]$fit[[which.min(BIC[,k])]]

 plot(1:(nK+1),c(0,devRatMinBIC),type="b",xlab="Number of clusters (K+1)",
      ylab="%Explained Variation")
