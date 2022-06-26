library(mclust)
library(CIMLR)
library(iClusterPlus)
library(MOFA2)
library(factoextra)
library(NbClust)
#clinical data for fitting Cox model
HCC_surv_df<-readRDS("HCCinputs/info.rds") #load clinical data
source("helpresfns.R")

#load the dataset
HCCfullz<-readRDS("HCCinputs/HCCfullz.rds")
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
#optimal number of clusters
length(unique(mclustfit$classification))
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

membfull <- cutree(hclust_avg, k = 3)
nonadjustedCox(membfull,HCC_surv_df)
adjustedCox(membfull,HCC_surv_df)

#########
#k-means#
#########


set.seed(100)
kmeansfit<-kmeans( pca_res$x[,1:npca], k)
membfull<-kmeansfit$cluster

#optimal number of clusters 2, ww method
fviz_nbclust(pca_res$x[,1:npca], kmeans, method = "wss") +
  labs(subtitle = "wss")

#perform clustering for k=6
set.seed(100)
kmeansfit<-kmeans(pca_res$x[,1:npca], 6)
membfull<-kmeansfit$cluster
nonadjustedCox(membfull,HCC_surv_df)
adjustedCox(membfull,HCC_surv_df)

#k=3
set.seed(100)
kmeansfit<-kmeans(pca_res$x[,1:npca], 3)
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


HCCmofafs<-readRDS("/Users/polinasuter/Downloads/HCC/HCC/HCCinputs/HCCfull.rds")
stdev<-c(0.1,0.5,2,1,2)
#pick only genes mutated more than in one sample
#HCCmofafs[[1]]<-HCCmofafs[[1]][which(apply(HCCmofafs[[1]],1,sum)>1),]
 for(i in 2:5){
   HCCmofafs[[i]]<-as.matrix(HCCmofafs[[i]])
  colnames(HCCmofafs[[i]])<-colnames(HCCmofafs[[1]])
  #select features with higher variance as recommended by the authors
  HCCmofafs[[i]]<-as.matrix(FSbyVar(HCCmofafs[[i]],stdev[i]))
 }

 MOFAobject <- create_mofa(HCCmofafs)
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
 model <- run_mofa(MOFAobject, "model1.hdf5")
 mofaclusters <- cluster_samples(model, k=3, factors="all")

 nonadjustedCox(mofaclusters$cluster,HCC_surv_df)
 adjustedCox(mofaclusters$cluster,HCC_surv_df)

 set.seed(200)
 wss<-vector()
 for(i in 1:10) {
   clusters <- cluster_samples(model, k=i, factors="all")
   wss[i]<-sum(clusters$withinss)
 }
 #elbow method
 plot(c(1:10),wss,type="b",col="blue")
#optimal number of clusters is 4
 mofaclusters <- cluster_samples(model, k=6, factors="all")
 nonadjustedCox(mofaclusters$cluster,HCC_surv_df)
 adjustedCox(mofaclusters$cluster,HCC_surv_df)

 ##############
 #iClusterPlus#
 ##############
 #We use a reduced set of features, same as for MOFA
 HCCcontz<-Reduce('rbind',HCCfullz[c("T","P","PP")])
 HCCcontz<-t(HCCcontz)


 set.seed(100)
 iclustfit2<- iClusterPlus(dt1=t(HCCfullz$M),dt2=t(HCCfullz$CN),dt3=HCCcontz,
                           type=c("binomial","gaussian","gaussian"),K=1, maxiter=20)
 set.seed(100)
 iclustfit3<- iClusterPlus(dt1=t(HCCfullz$M),dt2=t(HCCfullz$CN),dt3=HCCcontz,
                           type=c("binomial","gaussian","gaussian"),K=2, maxiter=20)
 set.seed(100)
 iclustfit4<- iClusterPlus(dt1=t(HCCfullz$M),dt2=t(HCCfullz$CN),dt3=HCCcontz,
                          type=c("binomial","gaussian","gaussian"),K=3, maxiter=20)
 set.seed(100)
 iclustfit5<- iClusterPlus(dt1=t(HCCfullz$M),dt2=t(HCCfullz$CN),dt3=HCCcontz,
                           type=c("binomial","gaussian","gaussian"),K=4, maxiter=20)


 alliclust<-list(iclustfit2,iclustfit3,iclustfit4,iclustfit5)
 which.min(unlist(lapply(alliclust,function(x)x$BIC)))

 membfull<-iclustfit3$clusters
 nonadjustedCox(membfull,HCC_surv_df)
 adjustedCox(membfull,HCC_surv_df)

 membfull<-iclustfit2$clusters
 nonadjustedCox(membfull,HCC_surv_df)
 adjustedCox(membfull,HCC_surv_df)
