# Source codes and necessary R libraries
source("repeated_simvfd.R")
source("iter_hpathD.R")
library(tidyverse)
library(dplyr)
library(MASS)
library(refund)
library(nlme)
library(boot)
library(parallel)
#library(Rmpi)
library(snow)
library(doSNOW)
library(doParallel)
library(foreach)
#library(InformationValue)
# Information on slides available for prostate cancer patients
infodat<-read.csv("peso_testset_mapping.csv")

# Q functions summarized from available images
QE_fun<-apply(as.matrix(read.csv("QEfun.csv")),2,function(w){
  w<-ifelse(w<1e-2,1e-2,w)
  w<-ifelse(w>(1-1e-2),(1-1e-2),w)
  logit(w)
})
QL_fun<-apply(as.matrix(read.csv("QLfun.csv")),2,function(w){
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
peso_dcom<-dcompSIMVFD(fDATA=cfun_dat,
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

peso_dcomPD<-dcompSIMVFD_PD(fDATA=fun_dat,
                            indxJ=peso_indxJ$Cont,
                            spat_indx=as.matrix(xRowCol),
                            funARG=r,
                            PVE=c(0.95,0.95,0.95),
                            corrM="spherical",
                            center_opt=1,
                            nknots=5,
                            order_penalty=2,
                            basisT=c("ps","ps"),
                            testFDATA=NULL,
                            test_indxJ=NULL)

# Setting of true parameters
vzeta<-apply(peso_dcom$Zeta,2,sd)
nEF<-ncol(peso_dcom$XiS)
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

csp_par<-lapply(cspatFIT, function(u){
  list(as.numeric(u$coefficients),
       simplify2array(getVarCov(u,type = "conditional"))[1:81,1:81],u$sigma,
       as.numeric(coef(u$modelStruct$corStruct,unconstrained = F)))
})

#csp_par[[1]][[2]]<-2.25
#csp_par[[1]][[3]]<-2.12; csp_par[[2]][[3]]<-0.545
#csp_par[[1]][[4]][1]<-0.450


ncspatFIT<-lapply(seq_len(ncol(peso_dcom$XiS)), function(u){
  sdata<-as.data.frame(cbind("Y"=ncspSCR[,u],ncspSCR[,-c(1:ncol(peso_dcom$XiS))]))
  gls(Y~1,data=sdata,correlation = corSpatial(form=~Row+Col|IndexJ,nugget = TRUE,type = corrM),control = lmeControl(maxIter = 200,msMaxIter = 200,tolerance = 1e-5,msTol = 1e-5))
})

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


# Cluster Construction
cls<-parallel::makeCluster(25)
registerDoParallel(cls)

# Running the simulation
Nsim<-100
n<-15;nc<-15;tS=10


pt<-proc.time()
par_res<-foreach(iter=1:Nsim,
                 .packages = c("tidyverse",
                               "dplyr","MASS",
                               "refund","nlme",
                               "boot"),
                 .export = c("r",
                             "xRowCol",
                             "mfun_dat",
                             "peso_dcom",
                             "vzeta",
                             "nEF",
                             "csp_par",
                             "canVM",
                             "ncsp_par",
                             "ncanVM"),
                 .combine = rbind) %dopar%
  iter_hpsim(iter,
             n=n,
             nc=nc,
             tS=tS,
             mi_pop=5:10,
             corr_model="spherical",
             err_var=peso_dcomPD$Sigma2)


npt<-proc.time()-pt
print(paste("Spherical case is done at",as.numeric(npt[3])))

par_resG<-foreach(iter=1:Nsim,
                 .packages = c("tidyverse",
                               "dplyr","MASS",
                               "refund","nlme",
                               "boot"),
                 .export = c("r",
                             "xRowCol",
                             "mfun_dat",
                             "peso_dcom",
                             "vzeta",
                             "nEF",
                             "csp_par",
                             "canVM",
                             "ncsp_par",
                             "ncanVM"),
                 .combine = rbind) %dopar%
  iter_hpsim(iter,
             n=n,
             nc=nc,
             tS=tS,
             mi_pop=5:10,
             corr_model="gaussian",
             err_var=peso_dcomPD$Sigma2)

ngpt<-proc.time()-npt
print(paste("Gaussian case is done at",as.numeric(ngpt[3])))

par_resE<-foreach(iter=1:Nsim,
                  .packages = c("tidyverse",
                                "dplyr","MASS",
                                "refund","nlme",
                                "boot"),
                  .export = c("r",
                              "xRowCol",
                              "mfun_dat",
                              "peso_dcom",
                              "vzeta",
                              "nEF",
                              "csp_par",
                              "canVM",
                              "ncsp_par",
                              "ncanVM"),
                  .combine = rbind) %dopar%
  iter_hpsim(iter,
             n=n,
             nc=nc,
             tS=tS,
             mi_pop=5:10,
             corr_model="exponential",
             err_var=peso_dcomPD$Sigma2)

nept<-proc.time()-ngpt
print(paste("Exponential case is done at",as.numeric(nept[3])))

#Stopping cluster 
parallel::stopCluster(cls)
#mpi.exit()
# Storing Simulation Results
dirP<-"/hpc/home/ma521/HPImage/"

res<-rbind(cbind(par_res,1),cbind(par_resG,2),cbind(par_resE,3))

write.table(res,file=paste(dirP,"himage_sim",n-tS,".txt",sep=""),col.names = FALSE,row.names = FALSE)
