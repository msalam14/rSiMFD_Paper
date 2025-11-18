# Required packages
library(tidyverse)
library(dplyr)
library(MASS)
library(refund)
library(matrixcalc)
library(boot)
library(mgcv)
library(e1071)
library(doParallel)
#library(Rmpi) # it is needed for running this script on cluster computer with 
# distributed memory
library(foreach)

# Required source file
source("repeated_simvfd.R")

# Required data input
infodat<-read.csv("peso_testset_mapping.csv")
r<-seq(10,90,7.25)/90
xq<-seq(1,2501,180)-1
xQ<-xq[-c(1,length(xq))]
yq<-seq(1,2501,180)-1
yQ<-xq[-c(1,length(xq))]
xRowCol<-data.frame(expand.grid(xQ,yQ))/2500
distM<-as.matrix(dist(xRowCol))
Lind<-which(lower.tri(distM),arr.ind = TRUE)
disV<-round(distM[Lind],2)
uD<-unique(disV)
bsz<-300
LocS<-250
#RGBmean<-read.table("RGBmean.txt",header = TRUE)

# Input of Q functions
QEfun<-as.matrix(read.csv("QEfun12.csv"))
QLfun<-as.matrix(read.csv("QLfun12.csv"))

QE_fun<-apply(QEfun,2,function(w){
  w<-ifelse(w<1e-2,1e-2,w)
  w<-ifelse(w>(1-1e-2),(1-1e-2),w)
  logit(w)
})

QL_fun<-apply(QLfun,2,function(w){
  w<-ifelse(w<1e-2,1e-2,w)
  w<-ifelse(w>(1-1e-2),(1-1e-2),w)
  logit(w)
})

c_opt<-TRUE


train_test_iter<-function(iter){
  ntrain<-30
  WIM<-rep(infodat$id,each=nrow(xRowCol))
  sID<-sample(unique(infodat$slide),ntrain)
  InfoDat<-infodat %>%filter(slide %in% sID)
  wsi<-InfoDat$id
  TsPOP<-infodat %>%filter(!slide %in%sID)
  tsPOP<-TsPOP 
  twsi<-tsPOP$id
  ########## 
  sld<-InfoDat$slide
  Sld<-rep(sld,each=nrow(xRowCol))
  Tsld<-tsPOP$slide
  TSld<-rep(Tsld,each=nrow(xRowCol))
  ##########
  trs<-c(60,80,120)
  do.call(rbind,lapply(trs, function(i){
    wind<-do.call(c,lapply(wsi[1:i],function(w){
      which(WIM==w)
    }))
    
    twind<-do.call(c,lapply(twsi,function(w){
      which(WIM==w)
    }))
    
    
    canS<-InfoDat$type[1:i]
    tcanS<-tsPOP$type
    
    #####################
    # Training data fromation
    QEfuN<-QE_fun[wind,]
    QLfuN<-QL_fun[wind,]
    trdata<-list(QEfuN,QLfuN) # training data
    #Test data formation
    TQEfuN<-QE_fun[twind,]
    TQLfuN<-QL_fun[twind,]
    tsdata<-list(TQEfuN,TQLfuN)
    ###############
    # Computation of Mean
    trMN<-lapply(lapply(trdata, colMeans),function(u){rbind(u,u)})
    # Centered training and test
    ctrDAT<-centerPHI(fDATA=trdata,findx = 2,meanF = trMN)
    ctsDAT<-centerPHI(fDATA=tsdata,findx = 2,meanF = trMN)
    
    fitD<-dcompSiMVFD(fDATA=ctrDAT$cfDATA,
                      indxJ=rep(4,i/4),
                      spat_indx=as.matrix(xRowCol),
                      funARG=r,
                      PVE=c(0.95,0.95),
                      corrM="spherical",
                      sub_center=c_opt,
                      nknots=5,
                      order_penalty=2,
                      basisT=c("ps","ps"),
                      testFDATA=ctsDAT$cfDATA,
                      test_indxJ=rep(4,40-ntrain),
                      control = lmeControl(maxIter = 300,msMaxIter = 300,tolerance = 1e-5,msTol = 1e-5))
    
    
    ### Performance of the predictors 
    trDAT<-data.frame(fitD$trainRD) %>% 
      mutate(Y=ifelse(canS=="cancer",1,0))
    #mutate(Rmean=RGBmean[InfoDat$id,1],
    #       Gmean=RGBmean[InfoDat$id,2],
    #       Bmean=RGBmean[InfoDat$id,3])
    tsDAT<-data.frame(fitD$testRD) %>% 
      mutate(Y=ifelse(tcanS=="cancer",1,0)) %>%
      set_names(colnames(trDAT))
    #  mutate(Rmean=RGBmean[InfoDat$id,1],
    #         Gmean=RGBmean[InfoDat$id,2],
    #         Bmean=RGBmean[InfoDat$id,3])
    
    gfit<-glm(Y~.,data = trDAT,family = binomial(link="logit"))
    # Training
    trPY<-predict(gfit,type="response")
    trGOF<-classSUM(actual=trDAT$Y,pred_score = trPY)
    # Test
    tsPY<-predict(gfit,newdata=tsDAT,type="response")
    tsGOF<-classSUM(actual=tsDAT$Y,pred_score = tsPY)
    trAUC<-as.numeric(pROC::auc(trDAT$Y,trPY))
    tsAUC<-as.numeric(pROC::auc(tsDAT$Y,tsPY))
    lg_res<-c(trGOF,trAUC,tsGOF,tsAUC,i,NA)
    
    ## gam modeling
    nknots<-c(5,8,10)
    gam_res<-t(sapply(nknots,function(h){
      fm<-"Y~"
      n_prd<-(ncol(trDAT)-1)
      for(b in 1:n_prd){
        if(b==1){
          fm<-paste(fm,"s(",colnames(trDAT)[b],",k=",h,",bs=",'"ps"',",m=2",")",sep="")
        } else{
          fm<-paste(fm,"+s(",colnames(trDAT)[b],",k=",h,",bs=",'"ps"',",m=2",")",sep="")
        }
      }
      frm<-formula(fm)
      gam_fit<-gam(frm,data = trDAT,family = binomial(link="logit"),method = "REML")
      rm(frm); rm(fm)
      gtrPY<-gam_fit$fitted.values
      gtrGOF<-classSUM(actual=trDAT$Y,pred_score = gtrPY)
      # Test
      gtsPY<-predict(gam_fit,newdata=tsDAT,type="response")
      gtsGOF<-classSUM(actual=tsDAT$Y,pred_score = gtsPY)
      gtrAUC<-as.numeric(pROC::auc(trDAT$Y,gtrPY))
      gtsAUC<-as.numeric(pROC::auc(tsDAT$Y,gtsPY))
      c(gtrGOF,gtrAUC,gtsGOF,gtsAUC,i,h)
    }))
    
    ## GAM: Tensor product
    tnknots<-c(4,6,8)
    tgam_res<-t(sapply(tnknots,function(h){
      fm<-"Y~"
      n_prd<-(ncol(trDAT)-1)/2
      for(b in 1:n_prd){
        if(b==1){
          fm<-paste(fm,"te(",colnames(trDAT)[b],",",colnames(trDAT)[b+1],",k=",h,",bs=",'"ps"',",m=2",")",sep="")
        } else{
          fm<-paste(fm,"+te(",colnames(trDAT)[(2*(b-1))+1],",",colnames(trDAT)[(2*(b-1))+2],",k=",h,",bs=",'"ps"',",m=2",")",sep="")
        }
      }
      frm<-formula(fm)
      tgam_fit<-gam(frm,data = trDAT,family = binomial(link="logit"),method = "REML")
      rm(frm); rm(fm)
      tgtrPY<-tgam_fit$fitted.values
      tgtrGOF<-classSUM(actual=trDAT$Y,pred_score = tgtrPY)
      # Test
      tgtsPY<-predict(tgam_fit,newdata=tsDAT,type="response")
      tgtsGOF<-classSUM(actual=tsDAT$Y,pred_score = tgtsPY)
      tgtrAUC<-as.numeric(pROC::auc(trDAT$Y,tgtrPY))
      tgtsAUC<-as.numeric(pROC::auc(tsDAT$Y,tgtsPY))
      c(tgtrGOF,tgtrAUC,tgtsGOF,tgtsAUC,i,h)
    }))
    
    
    ### SVM
    svm_basis<-c(1,2)
    svm_res<-t(sapply(svm_basis, function(b){
      svm_trDAT<-trDAT %>%
        mutate(Y=as.factor(ifelse(Y==0,-1,1)))
      svm_fit<-svm(Y~.,data =svm_trDAT,type="C-classification",kernel=ifelse(b==1,"linear","radial"),probability=TRUE)
      tr_svm_pred<-predict(svm_fit,newdata=svm_trDAT,probability=TRUE)
      strPY<-attr(tr_svm_pred,"probabilities")[,2]
      strGOF<-classSUM(actual=trDAT$Y,pred_score = strPY)
      ts_svm_pred<-predict(svm_fit,newdata=tsDAT,probability=TRUE)
      stsPY<-attr(ts_svm_pred,"probabilities")[,2]
      stsGOF<-classSUM(actual=tsDAT$Y,pred_score = stsPY)
      strAUC<-as.numeric(pROC::auc(trDAT$Y,strPY))
      stsAUC<-as.numeric(pROC::auc(tsDAT$Y,stsPY))
      c(strGOF,strAUC,stsGOF,stsAUC,i,b)
    }))
    
    cbind(rbind(c(lg_res,1),
                cbind(gam_res,2),
                cbind(tgam_res,3),
                cbind(svm_res,4)),c_opt)
  }))
}




cls<-parallel::makeCluster( (mpi.universe.size()-1) , type='MPI' )
registerDoParallel(cls)

Nsim<-100
par_res<-foreach(iter=1:Nsim,
                 .packages = c("tidyverse",
                               "dplyr","MASS",
                               "refund","nlme",
                               "boot","mgcv","e1071"),
                 .export = c("r",
                             "xRowCol",
                             "infodat",
                             "QE_fun",
                             "QL_fun",
                             "dcompSIMVFD",
                             "smth_covar",
                             "classSUM"),
                 .combine = rbind) %dopar%
  train_test_iter(iter)

parallel::stopCluster(cls)
mpi.exit()

dirP<- "PESO_Results/rSiMFD_PESO_DA/" #"/share/astaicu/malam3/HPImage/" needs to set for NCSU HPC

write.table(par_res,file=paste(dirP,"peso_dat_res12_",ifelse(c_opt,1,2),".txt",sep=""),col.names = FALSE,row.names = FALSE)

