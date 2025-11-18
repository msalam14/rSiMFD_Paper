# Source codes and necessary R libraries
source("SimulationStudy_RCodes/iter_hpathD.R")
library(tidyverse)
library(dplyr)
library(MASS)
library(refund)
library(mgcv)
library(pROC)
library(nlme)
library(boot)
source("SimulationStudy_RCodes/rSiMFD.R") #library(rSiMFD)
library(e1071)
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

load("SimulationStudy_Results/peso_dcom.RData")

load("SimulationStudy_Results/peso_dcomPD.RData")

# Setting of true parameters
vzeta<-apply(peso_dcom$Zeta,2,sd)
nEF<-ncol(peso_dcom$XiS)

load("SimulationStudy_Results/cspatFIT.RData")

csp_par<-lapply(cspatFIT, function(u){
  list(as.numeric(u$coefficients),
       simplify2array(getVarCov(u,type = "conditional"))[1:81,1:81],u$sigma,
       as.numeric(coef(u$modelStruct$corStruct,unconstrained = F)))
})

load("SimulationStudy_Results/ncspatFIT.RData")

ncsp_par<-lapply(ncspatFIT, function(u){
  list(as.numeric(u$coefficients),
       simplify2array(getVarCov(u,type = "conditional"))[1:81,1:81],u$sigma,
       as.numeric(coef(u$modelStruct$corStruct,unconstrained = F)))
})

# Adjustment of estimated parameters
ncsp_par[[1]][[4]][1]<-0.65
ncsp_par[[2]][[4]][1]<-1.25
ncsp_par[[1]][[3]]<-1.9
#ncsp_par[[2]][[3]]<-0.65

disD<-as.matrix(dist(xRowCol,diag = TRUE,upper = TRUE))


canVM<-lapply(csp_par, function(u){
  apply(disD,2,function(v){sapply(v,COVspher,Cov0=u[[3]]^2,range=u[[4]][1],nugget=u[[4]][2])})
})

ncanVM<-lapply(ncsp_par, function(u){
  apply(disD,2,function(v){sapply(v,COVspher,Cov0=u[[3]]^2,range=u[[4]][1],nugget=u[[4]][2])})
})


# Running the simulation
n<-15
nc<-15
tS<-10

pt<-proc.time()
par_res<-iter_hpsim(iter,
                    n=n,
                    nc=nc,
                    tS=tS,
                    mi_pop=5:10,
                    corr_model="spherical",
                    err_var=peso_dcomPD$Sigma2) 


npt<-proc.time()-pt
print(paste("Spherical case is done at",as.numeric(npt[3])))

par_resG<-iter_hpsim(iter,
                     n=n,
                     nc=nc,
                     tS=tS,
                     mi_pop=5:10,
                     corr_model="gaussian",
                     err_var=peso_dcomPD$Sigma2)

ngpt<-proc.time()-npt
print(paste("Gaussian case is done at",as.numeric(ngpt[3])))

par_resE<-iter_hpsim(iter,
                     n=n,
                     nc=nc,
                     tS=tS,
                     mi_pop=5:10,
                     corr_model="exponential",
                     err_var=peso_dcomPD$Sigma2)

nept<-proc.time()-ngpt
print(paste("Exponential case is done at",as.numeric(nept[3])))

# Storing Simulation Results
dirP<-"SimulationStudy_Results/"

res<-rbind(cbind(par_res,1,npt[3]),cbind(par_resG,2,ngpt[3]),cbind(par_resE,3,nept[3]))

write.table(res,file=paste(dirP,"himage_sim_",n-tS,"_",iter,".txt",sep=""),col.names = FALSE,row.names = FALSE)
