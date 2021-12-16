library(mclust)
library(CIMLR)
library(iClusterPlus)
library(MOFA)

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

#PCA
npca<-5
k<-3
var0<-which(apply(HccfullstackedPCAz,2,sd)==0)
if(length(var0>0)) HccfullstackedPCAz<-HccfullstackedPCAz[,-var0] else HccfullstackedPCAz<-HccfullstackedPCAz
pca_res <- prcomp(HccfullstackedPCAz, scale. = TRUE)


#mclust
set.seed(100)
mclustfit<-Mclust(pca_res$x[,1:npca],G=k)
membfull <- mclustfit$classification
#hclust
set.seed(100)
dist_mat <- dist( pca_res$x[,1:npca], method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'ward.D2')
membfull <- cutree(hclust_avg, k = k)
#kmeans
set.seed(100)
kmeansfit<-kmeans( pca_res$x[,1:npca], k)
membfull<-kmeansfit$cluster
#CIMLR all omics
set.seed(100)
cimplrfit<-CIMLR::CIMLR(HCCdataz,c=k,cores.ratio = 0) 
membfull<-cimplrfit$y$cluster
#CIMLR cont only
set.seed(100)
cimplrfit<-CIMLR::CIMLR(HCCdataz[c("T","P","PP")],c=k,cores.ratio = 0) #
membfull<-cimplrfit$y$cluster
#iclusterplus
set.seed(100)
iclustfit<- iClusterPlus(dt1=t(curdata$M),dt2=t(curdata$CN),dt3=HCCcontz,
                         type=c("binomial","gaussian","gaussian"),K=k-1, maxiter=20)
membfull<-iclustfit$clusters

#for survival analysis 
memb2check<-membfull #see Cox_fit.R

