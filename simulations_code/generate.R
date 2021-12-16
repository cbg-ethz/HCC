
#bnclust functions
#bnclust: MAP/sample option, structure learning GES/MCMC, cont, binary, mixed (all bge, bge+bde)

#bnbench functions
#-generate DAG continuous genDAG(type="cont, bin")
#-generate DAG mixed genDAGmixed(contpar=, binpar=)
#-generate data genDatacont, genDatabinm genDatamixed
#-generate mixture, limited number of simulation settings, reasonable ones
#-bnbench, clustbench function (nrep, ....)


#reasonable setting noise edgesignal, nodesignal, "strong/medium/weak", pctshd,k, ss,n


#omicsbn function
#black list penalization
#omiclearning class
#omicbn class
#plotting methods omicbn
#comparison methods omic/bn

paramfit<-function(dag,data) {
  res<-list()
  ordery<-orderdag(dag)
  if(!is.data.frame(data)) {
    data<-as.data.frame(data)
  }
  colnames(dag)<-colnames(data)
  for(i in ordery) {
    pari<-which(dag[,i]==1)
    vars<-colnames(dag)[pari]
    if(length(pari)>0) {
      parents<-paste(vars, collapse="+")
      fm <- as.formula(paste(paste("V", i, sep = ""), "~", parents,sep=""))
      res[[i]]<-lm(fm, data)
      #res[[i]]<-lm(as.formula(paste(paste("V", i, sep = ""), "~", paste(vars,sep="+"), sep="")), data)
    } else {
      res[[i]]<-lm(as.formula(paste(paste("V", i, sep = ""), "~ 1")), data)

    }
  }
  res
}
extractRsquared<-function(lmfit) {
  if(length(lmfit$coefficients)>1) {
    return(summary(lmfit)$r.squared)
  } else return(NA)
}
avRsquared<-function(bnfull, all=FALSE) {
  if(all) {
    return(mean(unlist(lapply(bnfull,function(x)summary(x)$r.squared))))
  } else {
    return(mean(unlist(lapply(bnfull,extractRsquared)), na.rm=TRUE))
  }
}
getmodel<-function(ep, p) {
  n<-ncol(ep)
  labels<-colnames(ep)
  incidence <- matrix(rep(0, n * n), nrow = n, ncol = n)
  colnames(incidence)<-rownames(incidence)<-labels
  incidence[which(ep > p)] <- 1
  return(incidence)
}
sortdag<-function(DAG,n) {
  adj<-t(graph2m(DAG))
  sort<-as.numeric(tsort(DAG))
  sortback<-vector()
  sortedadj<-matrix(rep(0,n*n),nrow=n,ncol=n)
  newedgeweights<-edgeWeights(DAG)
  from<-c()
  to<-c()
  for (i in 1:n) {sortback[i]<-which(sort==i)}

  for (i in 1:n){
    for (j in 1:n){
      if(adj[i,j]>0){
        colindex<-which(sort==j)
        rowindex<-which(sort==i)
        sortedadj[rowindex,colindex]<-1
      }
    }
  }
  return(m2graph(t(sortedadj),nodes=c(1:n)))
}
orderdag<-function(adj) {
  n<-ncol(adj)
  allnodes<-c(1:n)
  curnodes<-c(1)
  order<-c()
  cntr<-1
  while(length(curnodes)<n & cntr<n) {
    npar<-apply(adj,2,sum)
    curnodes<-which(npar==0)
    order<-c(setdiff(curnodes,order),order)
    adj[curnodes,]<-0
    cntr<-cntr+1
  }

  if(sum(adj)==0) return(order)
    else stop("not a DAG")

}
getrepr<-function(pdag,dag) {
  n<-ncol(pdag)
  for(i in 1:n) {
    for(j in 1:n) {
      if(pdag[i,j]==pdag[j,i] & pdag[i,j]==1) {
        pdag[i,j]<-dag[i,j]
        pdag[j,i]<-dag[j,i]
      }
    }
  }
  if(is.DAG(m2graph(pdag))) {
    return(pdag)
  } else {
    warning("not possible to resolve cycles, using MAP model!")
    return(dag)
  }
}
modifydag<-function(dag,shd) {
  n<-ncol(dag)
  ord<-orderdag(dag)
  edges<-which(dag==1)
  nedges<-length(edges)
  delE<-min(nedges,floor(shd/2))
  addE<-shd-delE

  deledges<-sample(edges,delE)
  addedges<-sample(ord[1:(n-1)],addE,replace=TRUE)
  for(i in addedges) {
    pos<-which(ord==i)
    pari<-ord[(pos+1):n]
    if(length(pari)==1) {
      newpar<-pari
      dag[newpar,i]<-1
    } else {
      newpar<-sample(pari,1)
      if(dag[newpar,i]!=1) {
        dag[newpar,i]<-1 } else {
          newpar<-sample(pari,1)
          dag[newpar,i]<-1
        }
    }
  }
  dag[deledges]<-0
  return(dag)

}
genDAG<-function(d, n, lB=0.1, uB=1, wmpct=0, wmmin=1,wmmax=3, shdpct=0, wm=NULL) {
  res<-list()
  if(shdpct>0) wmpct<-0
  if(is.null(wm)) {
     dag<-graph2m(pcalg::randomDAG(n, d*2/n, V=as.character(1:n)))
     edges<-which(dag==1)
     nedges<-sum(dag)
     wm<-dag
     edgeWeights<-runif(nedges, lB, uB)
     wm[edges]<-edgeWeights
  } else {
     dag<-matrix(0,nrow=n,ncol=n)
     edges<-which(abs(wm)>lB)
     dag[edges]<-1
     nedges<-length(edges)

      if(wmpct>0) {
        nme<-ceiling(wmpct*nedges/100)
        newedges<-sample(edges,nme)
        nweights<-wm[edges]
        nweights[newedges]<-sapply(nweights, function(x)x^(sample(c(1,-1),1)*runif(1,min=wmmin,max=wmmax)))
        wm[edges]<-sapply(nweights,max,lB)
      } else if(shdpct>0) {
        shdy<-max(ceiling(shdpct*nedges/100),2)
        newdag<-modifydag(dag,shdy) #get new dag with a defined shd with the start one
        newedges<-which(newdag==1) #edges of the new dag
        commonedges<-which((1*(dag & newdag))==1) #edges which are common in old/new dags
        newedges<-setdiff(newedges,commonedges)
        newn<-length(newedges)
        newwm<-matrix(0,nrow=n,ncol=n)
        newwm[commonedges]<-wm[commonedges]
        newwm[newedges]<-runif(newn, lB, uB)
        res$dag<-newdag
        res$wm<-newwm
        return(res)
      }
  }
  res$dag<-dag
  res$wm<-wm
  return(res)
}
genData<-function(dag, wm, ss, meanpct=0,uneqnodes=c("roots","random"),shiftmin=1,shiftmax=2,eqval=0, meaninit=NULL,
                  sigma=1, sigmasd=0, sigmas=NULL, type="cont") {
  n<-ncol(dag)
  npar<-apply(dag,2,sum)

  if(is.null(meaninit)) {
    means<-rep(eqval,n)
    orderx<-orderdag(dag)
    if(meanpct>0) {
      if(uneqnodes=="roots") {
        nm<-floor(meanpct*n/100)
        initmeans<-orderx[1:(n-nm)]
        newmeans<-orderx[(n-nm+1):n]
        means[initmeans]<-eqval
        means[newmeans]<-eqval+sample(c(1,-1),nm,replace=TRUE)*runif(nm,min=shiftmin,max=shiftmax)
      } else {
        nm<-floor(meanpct*n/100)
        newmeans<-sample(n,nm)
        initmeans<-c(1:n)[-newmeans]
        means[initmeans]<-eqval
        means[newmeans]<-eqval+sample(c(1,-1),nm,replace=TRUE)*runif(nm,min=shiftmin,max=shiftmax)
      }
    }
  } else {
    orderx<-orderdag(dag)
    if(meanpct>0) {
      if(uneqnodes=="roots") {
        nm<-floor(meanpct*n/100)
        initmeans<-orderx[1:(n-nm)]
        newmeans<-orderx[(n-nm+1):n]
        means[initmeans]<-init[initmeans]
        means[newmeans]<-sapply(init[newmeans],function(x)x+sample(c(1,-1),1)*runif(1,min=shiftmin,max=shiftmax))
        }else {
          newmeans<-sample(n,floor(meanpct*n/100))
          initmeans<-c(1:n)[-newmeans]
          means[initmeans]<-init[initmeans]
          means[newmeans]<- sapply(meaninit[newmeans],function(x)x+sample(c(1,-1),1)*runif(1,min=shiftmin,max=shiftmax))
        }

    } else means<-meaninit
  }
  if(is.null(sigmas)) sigmas<-abs(rnorm(n,sigma,sigmasd))

  #first generate weight matrix
  edges<-which(dag!=0)
  nedges<-length(edges)
  #define order of a dag
  ordery<-rev(orderx)


  #generate data
    datas<-matrix(nrow=ss, ncol=n)
    for(i in ordery) {
      if(npar[i]==0) {
        datas[,i]<-rnorm(ss,mean=means[i])
        } else {
            pari<-as.vector(which(dag[,i]!=0))
            datas[,i]<-wm[pari,i] %*% t(datas[,pari])+rnorm(ss,mean=means[i],sd=sigmas[i])
        }
    }
    if (type=="bin") {
    condmeans<-conditionalmeans(means, dag, ordery, npar)
    for(i in 1:n) {
      datas[,i]<-1*(datas[,i]>condmeans[i])
    }
  }

  res<-list()
  res$means<-means
  res$sigmas<-sigmas
  res$data<-datas

  res
}
#avpar - average number of binary parents for continuous nodes
genDAGmixed<-function(dagcont, dagbin, avchildren=0.5, lB=0.1, uB=1, childnodes=NULL, parentnodes=NULL) {
  res<-list()
  nbin<-ncol(dagbin$dag)
  ncont<-ncol(dagcont$dag)
  n<-nbin+ncont
  jointdag<-matrix(0,nrow=n,ncol=n)
  jointdag[1:nbin,1:nbin]<-dagbin$dag
  jointdag[1:ncont+nbin,1:ncont+nbin]<-dagcont$dag
  jointwm<-matrix(0,nrow=n,ncol=n)
  jointwm[1:nbin,1:nbin]<-dagbin$wm
  jointwm[1:ncont+nbin,1:ncont+nbin]<-dagcont$wm

  if(is.null(childnodes)) {
  nedges<-max(2,floor(avchildren*nbin))
  childnodes<-sample(1:ncont,nedges,replace=TRUE)
  parentnodes<-sample(1:nbin,nedges,replace=TRUE)
  } else {
    nedges<-length(childnodes)
  }

  for(i in 1:nedges) {
    jointdag[parentnodes[i],childnodes[i]+nbin]<-1
    jointwm[parentnodes[i],childnodes[i]+nbin]<-runif(1,min=lB,max=uB)
  }
  res$dag<-jointdag
  res$wm<-jointwm
  res$dagcont<-dagcont
  res$dagbin<-dagbin
  res$childnodes<-childnodes
  res$parentnodes<-parentnodes

  return(res)
}
genBinempty<-function(nbin, ss, par1=0.1, par2=0.5, dist="b") {
  if(dist=="u") freqs<-runif(nbin, min=par1, max=par2) else {
    freqs<-vector()
    freqs[1]<-rbeta(1,0.5,1)
    freqs[2:nbin]<-rbeta(nbin-1,par1,par2)
  }
  datas<-matrix(0,ncol=nbin,nrow=ss)
  for(i in 1:nbin) {
    datas[,i]<-rbinom(ss,1,freqs[i])
  }
  res<-list()
  res$datas<-datas
  res$freqs<-freqs
  return(res)
}
genDatamixed<-function(mixeddag, bindata, ss, meanpct=0,uneqnodes=c("roots","random"),shiftmin=1,shiftmax=2,eqval=0, meaninit=NULL,
                  sigma=1, sigmasd=0, sigmas=NULL, freqs=NULL) {
  dag<-mixeddag$dag
  wm<-mixeddag$wm
  n<-ncol(dag)
  nbin<-ncol(mixeddag$dagbin$dag)
  ncont<-ncol(mixeddag$dagcont$dag)
  #generate data
  datas<-matrix(nrow=ss, ncol=n)
  datas[,1:nbin]<-bindata
  dag[1:nbin,1:nbin]<-0
  wm[1:nbin,1:nbin]<-0

  npar<-apply(dag,2,sum)

  if(is.null(meaninit)) {
    meanscont<-rep(eqval,ncont)
    ordercont<-orderdag(mixeddag$dagcont$dag)
    if(meanpct>0) {
      if(uneqnodes=="roots") {
        nm<-floor(meanpct*ncont/100)
        initmeans<-ordercont[1:(ncont-nm)]
        newmeans<-ordercont[(ncont-nm+1):ncont]
        meanscont[initmeans]<-eqval
        meanscont[newmeans]<-eqval+sample(c(1,-1),nm,replace=TRUE)*runif(nm,min=shiftmin,max=shiftmax)
      } else {
        nm<-floor(meanpct*ncont/100)
        newmeans<-sample(ncont,nm)
        initmeans<-c(1:ncont)[-newmeans]
        meanscont[initmeans]<-eqval
        meanscont[newmeans]<-eqval+sample(c(1,-1),nm,replace=TRUE)*runif(nm,min=shiftmin,max=shiftmax)
      }
    }
  } else {
    ordercont<-orderdag(mixeddag$dagcont$dag)
    if(meanpct>0) {
      if(uneqnodes=="roots") {
        nm<-floor(meanpct*ncont/100)
        initmeans<-orderx[1:(ncont-nm)]
        newmeans<-orderx[(ncont-nm+1):ncont]
        meanscont[initmeans]<-init[initmeans]
        meanscont[newmeans]<-sapply(init[newmeans],function(x)x+sample(c(1,-1),1)*runif(1,min=shiftmin,max=shiftmax))
      }else {
        newmeans<-sample(ncont,floor(meanpct*ncont/100))
        initmeans<-c(1:ncont)[-newmeans]
        meanscont[initmeans]<-init[initmeans]
        meanscont[newmeans]<- sapply(meaninit[newmeans],function(x)x+sample(c(1,-1),1)*runif(1,min=shiftmin,max=shiftmax))
      }

    } else means<-meaninit
  }
  if(is.null(sigmas)) sigmas<-abs(rnorm(ncont,sigma,sigmasd))

  means<-c(freqs,meanscont)
  sigmas<-c(freqs,sigmas)
  absmeans<-vector()
  abssigmas<-vector()

  #first generate weight matrix
  edges<-which(dag!=0)
  nedges<-length(edges)
  #define order of a dag
  ordery<-rev(orderdag(dag))
  binnodes<-c(1:nbin)


  for(i in ordery) {
    if(i%in%binnodes) {
    } else if(npar[i]==0) {
      datas[,i]<-rnorm(ss,mean=means[i])
      absmeans[i]<-means[i]
      abssigmas[i]<-sigmas[i]
    } else {
      pari<-as.vector(which(dag[,i]!=0))
      datas[,i]<-wm[pari,i] %*% t(datas[,pari])+rnorm(ss,mean=means[i],sd=sigmas[i])
      absmeans[i]<-means[i]+sum(means[pari])
      abssigmas[i]<-sqrt(sigmas[i]^2+sum(sigmas[pari]^2))
    }
  }

  res<-list()
  res$means<-means
  res$sigmas<-sigmas
  res$absmeans<-absmeans
  res$abssigmas<-abssigmas
  res$data<-datas

  res
}
conditionalmeans<-function(means, dag, ordery, npar) {
  condmeans<-means
  for(i in ordery) {
    if(npar[i]!=0) {
      pari<-as.vector(which(dag[,i]!=0))
      condmeans[i]<-sum(condmeans[pari])
    }
  }
  return(condmeans)
}
learnmixed<-function(data, ncont, nbin, bgbin=NULL, bgcont=NULL, learnbin=FALSE, ...) {

  if(learnbin) {
    scorebin<-scoreparameters("bde",data[,1:nbin], bgnodes=bgbin)
    itfitbin<-iterativeMCMC(scorebin, ...)
  }

  scorecont<-scoreparameters("bge",data, bgnodes=c(1:nbin,bgcont))
  itfitcont<-iterativeMCMC(scorecont, ...)

  res<-list()
  res$dag<-itfitcont$DAG
  res$cpdag<-itfitcont$CPDAG

  if(learnbin) {
    res$dag[1:nbin,1:nbin]<-itfitbin$DAG
    res$cpdag[1:nbin,1:nbin]<-itfitbin$CPDAG
  }

  return(res)
}
samplemixed<-function(data, ncont, nbin, bgbin=NULL, bgcont=NULL, optimize=TRUE, learnbin=FALSE,...) {
  if(optimize) {

    if(learnbin) {
      scorebin<-scoreparameters("bde",data[,1:nbin], bgnodes=bgbin)
      itfitbin<-iterativeMCMC(scorebin, ...)
      samplebin<-orderMCMC(scorebin, chainout=TRUE, startspace = itfitbin$endspace,...)
      edgebin<-edgep(samplebin, pdag=TRUE)
    }

    scorecont<-scoreparameters("bge",data, bgnodes=c(1:nbin,bgcont))
    itfitcont<-iterativeMCMC(scorecont, ...)
    samplecont<-orderMCMC(scorecont, chainout=TRUE, startspace = itfitcont$endspace,...)
    edgecont<-edgep(samplecont, pdag=TRUE)
    edgecont[,1:nbin]<-0
    if(learnbin) edgecont[1:nbin,1:nbin]<-edgebin

    return(edgecont)
  } else {

    if(learnbin) {
      scorebin<-scoreparameters("bde",data[,1:nbin], bgnodes=bgbin)
      samplebin<-orderMCMC(scorebin, chainout=TRUE, ...)
      edgebin<-edgep(samplebin, pdag=TRUE)
    }

    scorecont<-scoreparameters("bge",data, bgnodes=c(1:nbin,bgcont))
    samplecont<-orderMCMC(scorecont, chainout=TRUE, ...)
    edgecont<-edgep(samplecont, pdag=TRUE)

    edgecont[,1:nbin]<-0
    edgecont[1:nbin,1:nbin]<-0
    if(learnbin) edgecont[1:nbin,1:nbin]<-edgebin

    return(edgecont)

  }

}
modellist<-function(pedges, p=c(0.5,0.6,0.7,0.8,0.9,0.99)) {
  res<-list()
  n<-ncol(pedges)
  k<-1
  for(i in p) {
    res[[k]]<-matrix(0,nrow=n,ncol=n)
    ones<-which(pedges>i)
    res[[k]][ones]<-1
    k<-k+1
  }
  return(res)
}
mixed2cp<-function(dag,nbin,ncont){
bindag<-dag[1:nbin,1:nbin]
bindag<-graph2m(dag2cpdag(m2graph(bindag)))
contdag<-dag
contdag[,1:nbin]<-0
contdag<-graph2m(dag2cpdag(m2graph(contdag)))
contdag[,1:nbin]<-0
contdag[1:nbin,1:nbin]<-bindag
return(contdag)
}
getmixed<-function(dag,nbin){
  dag[1:nbin,1:nbin]<-0
  return(dag)
}

#first for continuous variables only
genMixture<-function(k=3,type=c("cont","mixed"),centersignal=c("strong","medium","weak"),
                     sigma0,ssvec=rep(100,k),n=50,avpar=1, deltamu, lB, uB,
                     shdpct, randomseed=100, eqval=0,
                     mixedpar=list(nbin,avchildren,par1=0.1,par2=7,dist="b"),signalroots=TRUE,
                     returnpar=FALSE) {
  if(is.null(mixedpar$dist)) mixedpar$dist<-"b"
  dagl<-list()
  datal<-list()
  ss<-sum(ssvec)

  if (centersignal=="strong") {
    shiftmin<-1
    shiftmax<-2
  } else if (centersignal=="medium") {
    shiftmin<-0.5
    shiftmax<-1.5
  } else {
    shiftmin<-0.2
    shiftmax<-1
  }

 res<-list()
 pars<-list()
 if(type=="cont") {
  for(i in 1:k) {
    if(i==1) {
      dagl[[i]]<-genDAG(avpar,n,lB=lB,uB=uB)
    } else {
      dagl[[i]]<-genDAG(d,n,wm=dagl[[1]]$wm,shdpct=shdpct,lB=lB,uB=uB)
    }
    newdata<-genData(dagl[[i]]$dag,dagl[[i]]$wm,ss=ssvec[i],meanpct=deltamu, eqval=eqval, shiftmin=shiftmin,shiftmax=shiftmax,
                        meaninit=NULL, sigma=sigma0,sigmasd=0.2,uneqnodes="roots",type=type)
    datal[[i]]<-newdata$data
    if(returnpar) {
      pars[[i]]<-newdata[c("means","sigmas")]
    }
  }
 } else if(type=="mixed") {
   for(i in 1:k) {
      print(i)
     if(i==1) {
       bindag<-genDAG(0,mixedpar$nbin)
       bind<-genBinempty(mixedpar$nbin,ssvec[i], par1=mixedpar$par1, par2=mixedpar$par2, dist=mixedpar$dist)
       bindata<-bind$datas
       dagl[[i]]<-genDAG(avpar,n,lB=lB,uB=uB)
       dagl[[i]]<-genDAGmixed(dagl[[i]],bindag,avchildren=mixedpar$avchildren)

     } else {
       bind<-genBinempty(mixedpar$nbin,ssvec[i], par1=mixedpar$par1, par2=mixedpar$par2, dist=mixedpar$dist)
       bindata<-bind$datas
       dagl[[i]]<-genDAG(d,n,wm=dagl[[1]]$dagcont$wm,shdpct=shdpct,lB=lB,uB=uB)
       dagl[[i]]<-genDAGmixed(dagl[[i]],bindag,avchildren=mixedpar$avchildren,parentnodes=dagl[[1]]$parentnodes,
                              childnodes=dagl[[1]]$childnodes)

     }
     newdata<-genDatamixed(dagl[[i]],bindata=bindata,ss=ssvec[i],meanpct=deltamu, eqval=eqval, shiftmin=shiftmin,shiftmax=shiftmax,
                         meaninit=NULL, sigma=sigma0,sigmasd=0.2,uneqnodes="roots",freqs=bind$freqs)
     datal[[i]]<-newdata$data
     if(returnpar) {
       pars[[i]]<-newdata[c("means","sigmas","absmeans","abssigmas")]
     }
   }


 }



  if(returnpar) res$par<-pars
  res$membership<-as.vector(unlist(sapply(1:k,function(x)rep(x,ssvec[x]))))
  res$data<-Reduce('rbind',datal)
  res$r2<-vector()
  for(i in 1:k) {
    lmlocal<-paramfit(dagl[[i]]$dag,datal[[i]])
    res$r2[i]<-avRsquared(lmlocal)
  }

  res$DAGs<-list()
  for(i in 1:k) {
    res$DAGs[[i]]<-dagl[[i]]$dag
  }
  res$nbin<-mixedpar$nbin
  res$ssvec<-ssvec
  res$info<-list()
  res$info$n<-n
  res$info$avpar<-avpar
  res$info$shdpct<-shdpct
  res$info$lB<-lB
  res$info$uB<-uB
  res$info$type<-type
  res$info$deltamu<-deltamu
  res$info$randomseed<-randomseed
  res$info$k<-k

  return(res)
}


#bnomics package
#bnomics object, all placklists, penalizations, omics types, possible/not possible relationships
#clustering function based on bnclust, with parallelization
#plotting function
#getting various subnetworks
#getting marginals for T, CNA (graph?)
#comparing to a list of interations from a DB
#fitting parameters, getting list of high R2 nodes



