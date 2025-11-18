library(tidyverse)
library(dplyr)
library(MASS)
library(refund)
library(nlme)
library(boot)
library(parallel)
library(snow)
library(doSNOW)
library(doParallel)
library(foreach)
library(rSiMFD)
library(e1071)
library(doRNG)
# Information on slides available for prostate cancer patients
infodat<-read.csv("Data/peso_testset_mapping.csv")

# Q functions summarized from available images
QE_fun<-apply(as.matrix(read.csv("Data/QEfun.csv")),2,function(w){
  w<-ifelse(w<1e-2,1e-2,w)
  w<-ifelse(w>(1-1e-2),(1-1e-2),w)
  logit(w)
})
QL_fun<-apply(as.matrix(read.csv("Data/QLfun.csv")),2,function(w){
  w<-ifelse(w<1e-2,1e-2,w)
  w<-ifelse(w>(1-1e-2),(1-1e-2),w)
  logit(w)
})
# functional arguments
r<-seq(10,125,10)/125

# Spatial locations where RF is samples  
xq<-seq(1,2501,250)-1
xQ<-xq[-c(1,length(xq))]
yq<-seq(1,2501,250)-1
yQ<-xq[-c(1,length(xq))]
xRowCol<-as.matrix(expand.grid(xQ,yQ))/2500
rm(xq); rm(xQ); rm(yq); rm(yQ)

# Decomposition of PESO data

fun_dat<-list(QE_fun,QL_fun)
mfun_dat<-lapply(fun_dat,colMeans)
cfun_dat<-lapply(fun_dat,scale,scale=FALSE)
peso_indxJ<-infodat %>%
  group_by(slide) %>% 
  summarise(Cont=n()) %>% 
  dplyr::select(Cont)
peso_dcom<-MFPCA_HG_dcompSiMVFD(fDATA=cfun_dat,
                                indxJ=peso_indxJ$Cont,
                                spat_indx=as.matrix(xRowCol),
                                funARG=r,
                                PVE=c(0.95,0.95,0.95),
                                corrM="spherical",
                                center=FALSE,
                                nknots=5,
                                order_penalty=2,
                                basisT=c("ps","ps"),
                                sameEF="FALSE",
                                testFDATA=NULL,
                                test_indxJ=NULL)
save(peso_dcom,file="SimulationStudy_Results/peso_dcom.RData")

peso_dcomPD<-dcompSiMVFD(fDATA=fun_dat,
                         indxJ=peso_indxJ$Cont,
                         spat_indx=as.matrix(xRowCol),
                         funARG=r,
                         PVE=c(0.95,0.95,0.95),
                         corrM="spherical",
                         sub_center=TRUE,
                         nknots=5,
                         order_penalty=2,
                         basisT=c("ps","ps"),
                         testFDATA=NULL,
                         test_indxJ=NULL)

save(peso_dcomPD,file = "SimulationStudy_Results/peso_dcomPD.RData")


# Setting of true parameters
spSCR<-as.data.frame(peso_dcom$XiS) %>%
  mutate(Subject=rep(infodat$slide,each=nrow(xRowCol)),
         IndexJ=rep(seq_len(nrow(infodat)),each=nrow(xRowCol)),
         Row=rep(xRowCol[,1],times=nrow(infodat)),
         Col=rep(xRowCol[,2],times=nrow(infodat)),
         Type=rep(infodat$type,each=nrow(xRowCol)))

## Cancer
cspSCR<-spSCR %>%
  filter(Type=="cancer")
ncspSCR<-spSCR %>%
  filter(Type=="non-cancer")

corrM<-"spherical"

cspatFIT<-lapply(seq_len(ncol(peso_dcom$XiS)), function(u){
  sdata<-as.data.frame(cbind("Y"=cspSCR[,u],cspSCR[,-c(1:ncol(peso_dcom$XiS))]))
  gls(Y~1,data=sdata,correlation = corSpatial(form=~Row+Col|IndexJ,nugget = TRUE,type = corrM),control = lmeControl(maxIter = 200,msMaxIter = 200,tolerance = 1e-5,msTol = 1e-5))
})

save(cspatFIT,file = "SimulationStudy_Results/cspatFIT.RData")



ncspatFIT<-lapply(seq_len(ncol(peso_dcom$XiS)), function(u){
  sdata<-as.data.frame(cbind("Y"=ncspSCR[,u],ncspSCR[,-c(1:ncol(peso_dcom$XiS))]))
  gls(Y~1,data=sdata,correlation = corSpatial(form=~Row+Col|IndexJ,nugget = TRUE,type = corrM),control = lmeControl(maxIter = 200,msMaxIter = 200,tolerance = 1e-5,msTol = 1e-5))
})

save(ncspatFIT,file = "SimulationStudy_Results/ncspatFIT.RData")

