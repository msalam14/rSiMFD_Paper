iter_hpsim<-function(iter,n,nc,tS,mi_pop,corr_model="spherical",err_var=c(1,1)){
  set.seed(iter+n+nc+300)
  mi<-sample(mi_pop,(n+n),replace = TRUE)
  trn_sub<-1:((2*n)-(2*tS))
  ts_sub<-((2*n)-(2*tS)+1):(2*n)
  mi[ts_sub]<-5
  Y<-do.call(c,lapply(mi, function(u){
    rbinom(u,size=1,prob=0.5)
  }))   #rep("cancer",sum(mi))
  mFUN1<-outer(rep(1,sum(mi)*nrow(xRowCol)),mfun_dat[[1]])
  mFUN2<-outer(rep(1,sum(mi)*nrow(xRowCol)),mfun_dat[[2]])
  
  Zeta<-apply(sapply(vzeta, function(u){
    rnorm((n+n),mean = 0,sd=u)
  }),2,function(w){rep(w,mi*nrow(xRowCol))})
  
  sptRF<-sapply(seq_len(nEF), function(k){
    as.numeric(sapply(Y, function(u){
      if(u==1){
        mvrnorm(n=1,mu = rep(0,nrow(xRowCol)),Sigma = canVM[[k]]) 
      } else{
        mvrnorm(n=1,mu = rep(0,nrow(xRowCol)),Sigma = ncanVM[[k]]) 
      }
    }))
  })
  
  #cDR$spatMEAN[[k]]
  
  genFUN1<-mFUN1+(Zeta%*%t(peso_dcom$SubEigenF[[1]]))+(sptRF%*%t(peso_dcom$EigenF[[1]]))+matrix(rnorm(prod(dim(mFUN1)),mean = 0,sd=sqrt(err_var[1])),ncol=length(r))
  genFUN2<-mFUN2+(Zeta%*%t(peso_dcom$SubEigenF[[2]]))+(sptRF%*%t(peso_dcom$EigenF[[2]]))+matrix(rnorm(prod(dim(mFUN2)),mean = 0,sd=sqrt(err_var[2])),ncol=length(r))
  
  ## Training and and Test formation
  tr_size<-sum(mi[trn_sub])*nrow(xRowCol)
  tot_size<-sum(mi)*nrow(xRowCol)
  gFUN1<-genFUN1[1:tr_size,]
  gFUN2<-genFUN2[1:tr_size,]
  TgFUN1<-genFUN1[(tr_size+1):tot_size,]
  TgFUN2<-genFUN2[(tr_size+1):tot_size,]
  trY<-ifelse(Y[1:sum(mi[trn_sub])]==1,"cancer","non-cancer")
  tsY<-ifelse(Y[(sum(mi[trn_sub])+1):sum(mi)]==1,"cancer","non-cancer")
  # Centering of the data
  gfM<-lapply(list(gFUN1,gFUN2), colMeans)
  cgenDAT<-lapply(list(gFUN1,gFUN2),scale,scale=FALSE) 
  TcgenDAT<-list(TgFUN1,TgFUN2)
  tcgenDAT<-lapply(seq_len(length(cgenDAT)),function(u){
    as.matrix(TcgenDAT[[u]]-outer(rep(1,nrow(TcgenDAT[[u]])),gfM[[u]]))
  }) 
  sdrFIT<-dcompSIMVFD_PD(fDATA=list(gFUN1,gFUN2),
                      indxJ=mi[1:((2*n)-(2*tS))],
                      spat_indx=as.matrix(xRowCol),
                      funARG=r,
                      PVE=c(0.95,0.95,0.95),
                      corrM = corr_model,
                      sub_center=TRUE,
                      nknots=5,
                      order_penalty=2,
                      basisT=c("ps","ps"),
                      testFDATA = list(TgFUN1,TgFUN2),
                      test_indxJ = mi[((2*n)-(2*tS)+1):(2*n)], 
                      control = lmeControl(maxIter = 300,msMaxIter = 300,niterEM = 300))

  trDAT<-data.frame("Y"=ifelse(trY=="cancer",1,0),sdrFIT$trainRD)
  # model building
  # logistic
  gfit<-glm(Y~.,data = trDAT,family = binomial(link="logit"),control = glm.control(maxit=100))
  trPY<-predict(gfit,type="response")
  trGOF<-classSUM(actual=trDAT$Y,pred_score = trPY)
  # gam with tensor product smoothing
  
  
  # test set performance
  tsDAT<-data.frame("Y"=ifelse(tsY=="cancer",1,0),sdrFIT$testRD)
  tsPY<-predict(gfit,type="response",newdata=tsDAT)
  tsGOF<-classSUM(actual=tsDAT$Y,pred_score = tsPY)
  #c(trGOF,tsGOF)
  trAUC<-as.numeric(pROC::auc(trDAT$Y,trPY))
  tsAUC<-as.numeric(pROC::auc(tsDAT$Y,tsPY))
  c(trGOF,trAUC,tsGOF,tsAUC)
}