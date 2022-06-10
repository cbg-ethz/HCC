checkduplicated<-function(memb2check,pats,addex=c()) {
  dupind<-which(duplicated(pats))
  exindex<-c()

  if(length(dupind)>0) {
    pdup<-pats[dupind]

    for(i in pdup) {
      indlocal<-which(pats==i)
      dupval<-memb2check[indlocal]
      if(sd(dupval)!=0) exindex<-c(exindex,indlocal) else exindex<-c(exindex,indlocal[2])
    }
  }

  if(length(exindex)>0 | length(addex)>0) {
    return(c(1:length(memb2check))[-unique(c(exindex,addex))])
  } else return(c(1:length(memb2check)))

}
plot1nodecircle<-function(localint, node, y1=20, thresh=0.3,rmult=7,stringcheck=TRUE,cex=0.5, addint=NULL) {
  old.par<-par()
  on.exit(par(old.par))
  nint<-nrow(localint)
  bidir<-rep(0,nint)

  nodeint<-localint[which(localint$from==node | localint$to==node)[1],]
  if(nodeint$from==node) nodetype=nodeint$type1 else nodetype=nodeint$type2

  if(1>2) { #node%in%localint$to
    egint<-localint[which(localint$to==node),]
    logicDE<-as.character(egint[1,c(14:16)]<0.05)
    logicDE<-sapply( logicDE, function(x)paste(substr(strsplit(x, " ")[[1]], 1, 1), collapse='') , USE.NAMES=FALSE)
    logicDE<-paste(logicDE,collapse='')
    nodebord<-edgecol[logicDE]
  } else {
    nodebord<-typecols[nodetype]
  }

  comnodes<-setdiff(intersect(localint$from,localint$to),node)
  deleteind<-c()
  if(length(comnodes)>0) {
    for(i in comnodes) {
      cind<-which(localint$from==i | localint$to==i)
      bidir[cind]<-1
      pclmax<-apply(localint[cind,8:10],2,max)
      localint[cind[1],8:10]<-pclmax
      localint[cind[2],8:10]<-pclmax
      deleteind<-c(deleteind,cind[2])
    }
    localint<-localint[-deleteind,]
    bidir<-bidir[-deleteind]
  }

  nint<-nrow(localint)
  namesint<-c()
  typesint<-c()
  dirint<-c()
  logint<-c()
  loge<-c()

  ws<-c(1,1,1,2,2,2,3)
  names(ws)<-c("TFF","FTF","FFT","TTF","TFT","FTT","TTT")
  edgecol<-c("#F8766D", "#00BA38", "#619CFF","#800026","#88419d","#02818a","black")
  names(edgecol)<-c("TFF","FTF","FFT","TTF","TFT","FTT","TTT")

  for(i in 1:nrow(localint)) {
    if(localint$from[i]==node) {
      namesint<-c(namesint,localint$gene2[i])
      typesint<-c(typesint,localint$type2[i])
      dirint<-c(dirint,"f")
      logint<-c(logint,localint$db | (localint$gene1[i]==localint$gene2[i]))
      logicy<-as.character(localint[i,c(8:10)]>thresh)
      logicy<-sapply( logicy, function(x)paste(substr(strsplit(x, " ")[[1]], 1, 1), collapse='') , USE.NAMES=FALSE)
      logicy<-paste(logicy,collapse='')
      loge<-c(loge,logicy)



      if(localint$type2[i]=="PP") namesint[length(namesint)]<-paste(localint$gene2[i],sub('.*_',"\n",localint$to[i]),sep="")
    } else {
      namesint<-c(namesint,localint$gene1[i])
      typesint<-c(typesint,localint$type1[i])
      dirint<-c(dirint,"t")
      #logint<-c(logint,localint$string[i] | localint$TF[i] | localint$ks[i] | (localint$gene1[i]==localint$gene2[i])) #localint$TF[i] |
      logint<-c(logint,localint$db[i] | (localint$gene1[i]==localint$gene2[i])) #localint$TF[i] |
      logicy<-as.character(localint[i,c(8:10)]>thresh)
      logicy<-sapply( logicy, function(x)paste(substr(strsplit(x, " ")[[1]], 1, 1), collapse='') , USE.NAMES=FALSE)
      logicy<-paste(logicy,collapse='')
      loge<-c(loge,logicy)
      if(localint$type1[i]=="PP") namesint[length(namesint)]<-paste(localint$gene1[i],sub('.*_',"\n",localint$from[i]),sep="")
    }
  }

  if(localint$from[1]==node) {
    node<-localint$gene1[1]
    if(localint$type1[1]=="PP") node<-paste(node,sub('.*_',"\n",localint$from[1]),sep="")
    nodetype<-localint$type1[1]
  } else {
    node<-localint$gene2[1]
    if(localint$type2[1]=="PP") node<-paste(node,sub('.*_',"\n",localint$to[1]),sep="")
    nodetype<-localint$type2[1]
  }
  par(mar=rep(1,4))

  typecols<-c("#fbb4ae", "#fed9a6", "#ffffcc", "#b3cde3", "#decbe4")
  names(typecols)<-c("M","CN","T","P","PP")
  nint<-nrow(localint)
  deltar<-r*0.2
  x<-17.5
  y<-17.5
  r<-1.5
  if(nint>10) {
    R<-rmult*r
    Ri<-(rmult-1.2*r)*r
  } else {
    R<-rmult*r
    Ri<-rmult*r
  }


  #cex<-0.7
  plot(0,0, type="n",xaxt="n", yaxt="n", xlab="", ylab="",
       xlim=c(0,35), ylim=c(0,35),bty="n")

  draw.circle(x=x, y=y, r=r,col=typecols[nodetype], border=nodebord)
  text(x=x,y=y,node,col="red",cex=cex)
  centers<-matrix(nrow=nint,ncol=2)
  colnames(centers)<-c("x","y")
  if(nint<10) {
    deltaphi<-629/(nint)
  } else {
    deltaphi<-629/(nint+1)
  }
  fromRx<-x+2*r
  fromRy<-y+2*r
  toRx<-x+R-r
  toRy<-y+R-r
  for(i in 1:nint) {
    if(i%%2==0) {
      centers[i,"x"]<-R*cos(deltaphi*i/100)
      centers[i,"y"]<-R*sin(deltaphi*i/100)
    } else {
      centers[i,"x"]<-Ri*cos(deltaphi*i/100)
      centers[i,"y"]<-Ri*sin(deltaphi*i/100)
    }
    if(1>2){#dirint[i]=="f"
      egint<-localint[i,]
      logicDE<-as.character(egint[1,c(14:16)]<0.05)
      logicDE<-sapply( logicDE, function(x)paste(substr(strsplit(x, " ")[[1]], 1, 1), collapse='') , USE.NAMES=FALSE)
      logicDE<-paste(logicDE,collapse='')
      bordcol<-edgecol[logicDE]
    } else {
      bordcol<-typecols[typesint[i]]
    }

    draw.circle(x=x+centers[i,"x"], y=y+centers[i,"y"], r=r,col=typecols[typesint[i]], border=bordcol)
    textt(x=x+centers[i,"x"],y=y+centers[i,"y"],namesint[i],F,cex=cex)
    slopey<-(centers[i,2]-y)/(centers[i,1]-x)
    if(dirint[i]=="f") {
      sz<-ifelse(bidir[i]==1,0,0.5)
      wh<-ifelse(bidir[i]==1,0,0.8)
      iArrows(r*cos(deltaphi*i/100)+x,r*sin(deltaphi*i/100)+y,
              -r*cos(deltaphi*i/100)+x+centers[i,"x"],-r*sin(deltaphi*i/100)+y+centers[i,"y"],
              sh.lty=ifelse(logint[i],1,2),
              sh.lwd=ws[loge[i]], sh.col=edgecol[loge[i]], size=sz,width=wh)
    } else {
      sz<-ifelse(bidir[i]==1,0,0.5)
      wh<-ifelse(bidir[i]==1,0,0.8)
      iArrows(-r*cos(deltaphi*i/100)+x+centers[i,"x"],-r*sin(deltaphi*i/100)+y+centers[i,"y"],
              r*cos(deltaphi*i/100)+x,r*sin(deltaphi*i/100)+y,
              sh.lty=ifelse(logint[i],1,2),
              sh.lwd=ws[loge[i]], sh.col=edgecol[loge[i]], size=sz,width=wh)
    }
  }
  if(!is.null(addint)) {
    i<-which(localint$gene1==addint$from | localint$gene2==addint$from)[1]
    j<-4
    print(i)
    print(j)
    iArrows(r+x+centers[i,"x"],y+centers[i,"y"],
            x+centers[j,"x"],r+y+centers[j,"y"],
            sh.lty=1,
            sh.lwd=ws["TFF"], sh.col=edgecol["TFF"], size=0.5,width=0.8)
  }

}

FSbyVarmy<-function(dataf,value) {
  myvars <- apply(dataf,1, var,na.rm=TRUE)
  myvars <- sort(myvars,decreasing=TRUE)
  myvars <- myvars[1:value]
  return(dataf[names(myvars),])
}

nonadjustedCox<-function(memb,survdf) {
  pats<-survdf$patientID
  curru<-checkduplicated(memb,pats) #check patients with duplicated biopsies (same cluster or not, if not - exclude both, if same include one)
  survobj<-survdf[curru,]
  survobj<-cbind(survobj,as.factor(memb[curru]))
  nrow(survobj) #how many samples were included
  colnames(survobj)[ncol(survobj)]<-"cluster"
  surv_object <- Surv(time = survobj$survival_time, event = survobj$death)
  fitmod <- coxph(surv_object ~  cluster, survobj)
  #here is the result for non-adjusted model
  return(summary(fitmod)$logtest[3])

}
adjustedCox<-function(memb,survdf) {
  pats<-survdf$patientID
  #for adjusted model, exclude patient 41, group '0'
  curru<-checkduplicated(memb,pats,addex=c(41)) #check patients with duplicated biopsies (same cluster or not, if not - exclude both, if same include one)
  survobj<-survdf[curru,]
  survobj<-cbind(survobj,as.factor(memb[curru]))
  nrow(survobj) #how many samples were included
  colnames(survobj)[ncol(survobj)]<-"cluster"
  surv_object <- Surv(time = survobj$survival_time, event = survobj$death)
  fitbase <- coxph(surv_object ~  bclc, survobj)
  fitbn <- update(fitbase, . ~ . + cluster)
  #here is the result for adjusted model
  fit<-anova(fitbase, fitbn)
  return(fit$`P(>|Chi|)`[2])
}
