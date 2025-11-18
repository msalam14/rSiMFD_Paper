#' This function generates repeated spatially indexed multivariate functional 
#' where the spatial domain is [0,1]^2. It considers a spacial case of the general
#' situation described in Alam and Staicu (20xx). Specifically, this functions generates
#' data where identical spatial index is used for every subject and functions are
#' sampled at identical points.
#' 
#' @param n number of subject for which data would be generated
#' @param mi is a vector of length n consisting of number of occasions each 
#' subject would have data
#' @param spat_indx is a L by 2 matrix of spatial index
#' @param funARG is a numerical vector of functional argument
#' @param meanFUN is a  function that returns a p-dimensional mean vector
#' @param basisFUN is a list of p functions. Each function provides values of basis vectors. 
#' @param zetaVAR is a vector of length K and is consist of the variance of zeta.
#' @param corrFUN is a surrogate function for a spatial correlation function. 
#' It takes distance as input and yields value of the correlation.
#' @param spatVAR is scalar represents the variance parameter for the spatially indexed scores xi
#' @param sigma2 is a scalar represents the error variance
#' @example examples/example_rsimvFD.R
#' @export
rsimvFD<-function(n,mi,spat_indx,funARG,meanFUN,basisFUN,zetaVAR,corrFUN,spatVAR,sigma2){
  M<-sum(mi)
  p<-length(basisFUN)
  K<-length(zetaVAR)
  L<-nrow(spat_indx)
  H<-length(funARG)
  zCOEF<-sapply(zetaVAR, function(sk){rep(rnorm(n = n,mean = 0,sd = sqrt(sk)),mi*L)})
  mnVEC<-sapply(funARG,meanFUN)
  xiK<-sapply(1:K, function(k){
    cvM<-spatVAR[k]*(corrFUN[[k]](as.matrix(dist(spat_indx))))
    mnV<-rep(0,nrow(spat_indx))
    as.numeric(t(rmvnorm(M,mean=mnV,sigma=cvM)))
  })
  basisV<-lapply(1:p, function(j){
    basisFUN[[j]](funARG)
  })
  gen_data<-lapply(1:p,function(j){
    outer(rep(1,M*L),mnVEC[j,])+zCOEF%*%t(basisV[[j]])+xiK%*%t(basisV[[j]])+matrix(rnorm(M*L*H),ncol = H)
  })
  list("subject_index"=rep(1:n,mi*L),
       "patch_index"=sapply(1:n,function(i){rep(1:mi[i],each=L)}),
       "data"=gen_data)
}


#' Smooth covariance and cross-covariance for a zero-mean multivariate functional data
#' 
#' @param funDATA a list with elements as univariate observed functional data
#' @param funARG arguments of the observed function
#' @param pve a scalar represents the threshold for percentage of variation 
#' explained in the univariate FPCA for smoothing auto-covariance
#' @param nKnot a scalar representing the number of knots to be used in fpca.face and tensor product smoothing for cross-covariance
#' @param basisT is a vector of two elements contains the types of basis be used in tensor product smoothing
#' @param tps_m order of the penaly for tensor-product smoothing for cross-covariance
#' @export
smth_covar<-function(funDATA,funARG,pve,nKnot,basisT,tps_m){
  p<-length(funDATA)
  smth_cov<-Matrix::Matrix(0,nrow=p*length(funARG),ncol=p*length(funARG))
  for(u in 1:p){
    for(w in u:p){
      cvMAT<-(1/nrow(funDATA[[u]]))*(t(funDATA[[u]])%*%funDATA[[w]])
      a1<-((u-1)*length(funARG))+1
      a2<-u*length(funARG)
      b1<-((w-1)*length(funARG))+1
      b2<-w*length(funARG)
      if(u==w){
        ufpc<-fpca.face(Y=funDATA[[u]],argvals = funARG,pve = pve,knots = nKnot)
        smth_cov[a1:a2,b1:b2]<-Reduce(`+`,lapply(seq_len(ufpc$npc), function(l){
          ((1/sqrt(length(funARG)))*ufpc$evalues[l])*outer(as.matrix(ufpc$efunctions)[,l],as.matrix(ufpc$efunctions*sqrt(length(funARG)))[,l])
        }))
      } else{
        posMC<-as.matrix(expand.grid(1:length(funARG),1:length(funARG)))
        datMC<-data.frame("row"=funARG[posMC[,1]],
                          "col"=funARG[posMC[,2]],
                          "FunV"=cvMAT[posMC])
        frm<-as.formula(paste("FunV~te(row,col,k=c(nKnot,nKnot),bs=basisT,m=tps_m)"))
        smFIT<-mgcv::bam(frm,data = datMC)
        smth_cov[a1:a2,b1:b2]<-matrix(smFIT$fitted.values,nrow = length(funARG))
        smth_cov[b1:b2,a1:a2]<-t(matrix(smFIT$fitted.values,nrow = length(funARG)))
      }
    }
  }
  as.matrix(smth_cov)
}


#' Calculates accuracy, sensitivity and specificity in a binary classification.
#' 
#' @param actual a vector actual class in 1 and 0 representing class of interest
#' and the other.
#' @param pred_score is the predicted probability for class 1
#' @param thresh is a cut off for class 1, the default is 0.5
#' @return a vector of three values, representing accuracy, sensitivity,
#'  and specificity
#' @example examples/example_classSUM.R
#' @export
classSUM<-function(actual,pred_score,thresh=0.5){
  pred<-ifelse(pred_score>=thresh,1,0)
  res<-as.numeric(abs(actual-pred)==0)
  c(mean(res),mean(res[actual==1]),mean(res[actual==0]))
}

#' Calculate the spherical covariance at a given distance when range, nugget and
#' total variance are known.
#' 
#' @param h the spatial distance where covariance to be calculated
#' @param Cov0 a scalar representing the total variance
#' @param range a scalar representing the range parameter
#' @param nugget a scalar representing the nugget parameter
#' @return a scalar representing the value of spatial covariance
#' @example examples/example_COVspher.R
#' @export
COVspher<-function(h,Cov0,range,nugget){
  if(h<0)
    stop("h must be non-negative")
  if(h==0){
    Cov0
  } else if(h<=range){
    crv<-(1-(1.5*(h/range))+(0.5*(h/range)^3))
    crv*(1-nugget)*Cov0
  } else{
    0
  }
}

#' Calculate different spatial autocorrelation functions, Spherical, Gaussian 
#' and Exponential, at a given distance, assuming range and nugget parameters are known. 
#' 
#' @param r the spatial distance where covariance to be calculated
#' @param d a scalar representing the range parameter
#' @param nugget a scalar representing the nugget parameter
#' @param type a character, spherical or gaussian, that specifies which type of 
#' spatial correlation needs to be calculated. If anything else is chosen it
#' computes the exponential correlation.  
#' @return a scalar representing the value of spatial correlation
#' @example examples/example_corSpatialVal.R
#' @export
corSpatialVal<-function(r,d,nugget=0,type="spherical"){
  if(nugget>0){
    if(type=="spherical"){
      ifelse(r>d,0,(1-nugget)*(1-1.5*(r/d)+0.5*(r/d)^3))
    } else if(type=="gaussian"){
      (1-nugget)*exp(-(r/d)^2)
    } else {
      (1-nugget)*exp(-r/d)
    }
  } else{
    if(type=="spherical"){
      ifelse(r>d,0,1-1.5*(r/d)+0.5*(r/d)^3)
    } else if(type=="gaussian"){
      exp(-(r/d)^2)
    } else {
      exp(-r/d)
    } 
  }
}

#' Calculate different spatial autocovariance functions, Spherical, Gaussian 
#' and Exponential, at a given distance, assuming range and nugget parameters are known. 
#' 
#' @param r the spatial distance where covariance to be calculated
#' @param d a scalar representing the range parameter
#' @param nugget a scalar representing the nugget parameter
#' @param sill_cov a scalar representing the sill variance
#' @param type a character, spherical/gaussian/exponential, that specifies which type of 
#' spatial correlation needs to be calculated. If anything else is chosen it
#' computes the exponential correlation.  
#' @return a scalar representing the value of spatial covariance
#' @example examples/example_covSpatialVal.R
#' @export
covSpatialVal<-function(r,d,nugget=0,sill_cov,type="spherical"){
  if(type=="spherical"){
    ifelse(r==0,sill_cov,ifelse(r>d,0,sill_cov*(1-nugget)*(1-1.5*(r/d)+0.5*(r/d)^3)))
  } else if(type=="gaussian"){
    ifelse(r==0,sill_cov,sill_cov*(1-nugget)*exp(-(r/d)^2))
  } else if(type=="exponential"){
    ifelse(r==0,sill_cov,sill_cov*(1-nugget)*exp(-r/d))
  } else{
    stop("type should any of spherical, gaussian or exponential")
  }
}


#' Performs multivariate functional principal component analysis on subject-
#' specific means for a given iMVFD. This is a function used in the simulation 
#' study only for estimating the variances of scores involved with subject-
#' specific deviation. 
#' 
#' @param fDATA is a list of observed functional data. Each element represents 
#' one component of the multivariate function
#' @param indxJ is a vector of length n that has elements represent m_i
#' @param spat_indx is a L by 2 matrix contains spatial locations
#' @param funARG is vector of sampled points where functional data were observed. 
#' This can be generalized for sparse functional data. However this version 
#' is only for dense functional data. Also all the functions are observed at 
#' identical sampled point.
#' @param PVE is a vector of length 2 for percentage of variation explained. 
#' The first is used for fpca of each functional variable. The second is for 
#' eigen decomposition of scores from each of the 
#' @param center is a logical vector of length 1 indicates whether the 
#' functional data need centering by the sample mean
#' @return a list with following objects:
#' \itemize{
#'  \item n : number of subjects in the data
#'  \item p : number of features in the data
#'  \item funARG : functional grid where FD were observed
#'  \item mFDATA : mean of each feature where center=TRUE
#'  \item SubEigenF : estimated eigenfunctions
#'  \item Zeta : projection of subject means on to the estimated eigenfunctions
#' }
#' @example examples/example_mfpca_HG.R
#' @export
mfpca_HG<-function(fDATA,indxJ,spat_indx,funARG,PVE,center=FALSE){
  n<-length(indxJ) 
  p<-length(fDATA)
  L<-nrow(spat_indx)
  if(center){
    fDATA<-lapply(fDATA, function(u){
      scale(u,center = TRUE,scale = FALSE)
    })
  }
  # Subject effect: basis and scores
  subID<-rep(1:n,indxJ*L)
  subMEAN<-lapply(fDATA, function(u){
    apply(u,2,function(v){
      as.numeric(tapply(v,subID,mean))
    })
  })
  
  # Univariate fpca of subject-specific means
  umFIT<-lapply(seq_len(length(subMEAN)), function(u){
    fpca.sc(Y=subMEAN[[u]],argvals = funARG,pve = PVE[1])
  })
  
  # estimation of univariate scores
  uzeta_vec<-do.call(cbind,lapply(seq_len(length(subMEAN)),function(u){
    (1/length(funARG))*((subMEAN[[u]]-outer(rep(1,n),umFIT[[u]]$mu))%*%umFIT[[u]]$efunctions)
  }))
  
  spd_uzeta<-eigen(cov(uzeta_vec))
  
  npc_uz<-which((cumsum(spd_uzeta$values)/sum(spd_uzeta$values))>PVE[2])[1]
  
  CMat<-spd_uzeta$vectors[,1:npc_uz]
  
  # subject-effect eigen-function
  sub_basisF<-split.data.frame(
    as.matrix(Matrix::bdiag(lapply(umFIT,function(u){u$efunctions}))%*%CMat),
    rep(1:p,each=length(funARG)))
  
  # Subject-specific scores
  Zeta<-Reduce(`+`,lapply(1:p,function(u){
    (1/length(funARG))*((subMEAN[[u]]-outer(rep(1,n),umFIT[[u]]$mu))%*%sub_basisF[[u]])
  }))
  
  mfDATA<-if(center){
    lapply(fDATA,colMeans)
  } else{
    NULL
  }
  
  
  retOBJ<-list(
    "n" = n,
    "p" = p,
    "funARG" = funARG,
    "L" = nrow(spat_indx),
    "mFDATA"=mfDATA,
    "SubEigenF"=sub_basisF,
    "Zeta" = Zeta)
  retOBJ
}



#' Train a classifier based multivariate FPCA (separate basis) of spatially indexed multivariate functional data
#' 
#' @param fDATA is a list of functional data, length is equal to the number of function
#' @param indxJ is a vector of length n that has elements represent m_i
#' @param spat_indx is a L by 2 matrix contains spatial locations
#' @param funARG is vector of sampled points where functional data were observed. 
#' This can be generalized for sparse functional data. However this version is only for
#' dense functional data. Also all the functions are observed at identical sampled point.
#' @param PVE is a vector of length 2 for percentage of variation explained. 
#' The first is used for fpca of each functional variable. The second is for 
#' eigen decomposition of scores from each of the 
#' @param corrM is a character specifies the name of spatial correlation
#' @param sub_center is a logical vector of length 1 indicates whether the 
#' functional data would be centered using subject-specific means. If FALSE, 
#' grand mean will be subtracted.
#' @param testFDATA is a list of functional data from the test for which we want 
#' to obtain the projection onto the reduced dimension
#' @param test_indxJ is similar to indxJ, however, it corresponds to the test data provided
#' @return a list with following objects:
#' \itemize{
#'  \item n : number of subjects in the data
#'  \item p : number of features in the data
#'  \item funARG : functional grid where FD were observed
#'  \item mFDATA : mean of each feature where center=TRUE
#'  \item SubEigenF : estimated eigenfunctions
#'  \item Zeta : projection of subject means on to the estimated eigenfunctions
#' }
#' @example examples/example_dcompSiMVFD.R
#' @export
dcompSiMVFD<-function(fDATA,indxJ,spat_indx,funARG,PVE,corrM,sub_center=TRUE,nknots=5,
                         order_penalty=2,basisT=c("ps","ps"),
                         testFDATA=NULL,test_indxJ=NULL, ...){
  if(!is.null(testFDATA)&is.null(test_indxJ))
    stop("For prediction test_indxJ is required")
  n<-length(indxJ) 
  p<-length(fDATA)
  L<-nrow(spat_indx)
  if(sub_center){
    subID<-rep(1:n,indxJ*L)
    cfDATA<-lapply(fDATA, function(u){
      apply(u,2,function(v){
        v-rep(as.numeric(tapply(v,subID,mean)),indxJ*L)
      })
    })
  } else{
    cfDATA<-lapply(fDATA, function(u){
      scale(u,center = TRUE,scale = FALSE)
    })
  }
  
  if(!sub_center){
    ctrMEAN<-lapply(fDATA, function(u){
      as.numeric(colMeans(u))
    })
  }
  
  us_var<-as.numeric(do.call(c,lapply(cfDATA, function(u){
    diag(t(u)%*%u)/(nrow(u)-1)
  })))
  
  
  smth_cov<-smth_covar(funDATA=cfDATA,
                       funARG=funARG,
                       pve=PVE[1],
                       nKnot=nknots,
                       basisT=basisT,
                       tps_m=order_penalty)
  
  # Spectral decomposition of smooth covariance
  Wmat<-diag(rep(1/sqrt(length(funARG)),2*length(funARG)))
  IWmat<-diag(rep(sqrt(length(funARG)),length(funARG)))
  spd_smth_covar<-eigen((Wmat%*%(smth_cov%*%Wmat)))
  plamd<-spd_smth_covar$values[spd_smth_covar$values>0]
  sp_nbasis<-which((cumsum(plamd)/sum(plamd))>PVE[2])[1]
  sp_basis<-lapply(split.data.frame((spd_smth_covar$vectors[,1:sp_nbasis]),rep(1:p,each=length(funARG))),function(w){IWmat%*%w})
  sp_score<-(do.call(cbind,fDATA)%*%do.call(rbind,sp_basis))*(1/length(funARG))
  
  # Fitting of the projected data
  
  mfpcaSCR<-as.data.frame(sp_score)
  mfpcaSCR$Subject<-rep(1:n,indxJ*L)
  mfpcaSCR$IndexJ<-as.factor(rep(1:sum(indxJ),each=L))
  mfpcaSCR$Row<-rep(spat_indx[,1],times=sum(indxJ))
  mfpcaSCR$Col<-rep(spat_indx[,2],times=sum(indxJ))
  
  
  spatFIT<-lapply(seq_len(ncol(sp_score)), function(u){
    sdata<-as.data.frame(cbind("Y"=mfpcaSCR[,u],mfpcaSCR[,-c(1:ncol(sp_score))]))
    lme(fixed=Y~1,data=sdata,random=~1|Subject,correlation = corSpatial(form=~Row+Col|Subject/IndexJ,nugget = TRUE,type = corrM), ...)
  })
  
  dFACT<-sapply(spatFIT, function(u){
    as.numeric(u$coefficients$fixed)*rep(1,nrow(sp_score))+
      rep(as.numeric(u$coefficients$random$Subject),indxJ*L)
  })
  
  ## Construction of marginal covariance
  disD<-as.matrix(dist(spat_indx,diag = TRUE,upper = TRUE))
  vcovSP<-lapply(spatFIT, function(u){
    cor_par<-as.numeric(coef(u$modelStruct$corStruct,unconstrained = F))
    apply(disD,2,function(v){sapply(v,covSpatialVal,d=cor_par[1],nugget=cor_par[2],sill_cov=u$sigma^2,type=corrM)})
  })
  

  scMAT<-do.call(cbind,lapply(seq_len(sp_nbasis),function(v){
    t(sapply(split(sp_score[,v]-dFACT[,v], mfpcaSCR$IndexJ),function(u){
      c(mean(u),mean(u^2))
    }))
  }))
  
  ## Test data
  
  if(!is.null(testFDATA)){
    # Projected test data
    tn<-length(test_indxJ)
    tsubID<-rep(1:tn,test_indxJ*L)
    tsp_score<-(do.call(cbind,testFDATA)%*%do.call(rbind,sp_basis))*(1/length(funARG))
    
    tspDAT<-split.data.frame(tsp_score,tsubID)
    
    tsREF<-apply(t(sapply(seq_len(length(tspDAT)), function(v){
      OmL<-rep(1,test_indxJ[v]*L)
      ImM<-diag(test_indxJ[v])
      sapply(seq_len(length(spatFIT)), function(u){
        tauK<-(as.numeric(attr(corMatrix(spatFIT[[u]]$modelStruct[[1]])[[1]],"stdDev")*spatFIT[[u]]$sigma))^2
        covM<-(tauK*outer(OmL,OmL))+kronecker(ImM,vcovSP[[u]])
        as.numeric(matrix(tauK*OmL,nrow=1)%*%(solve(covM)%*%matrix(tspDAT[[v]][,u]-as.numeric(spatFIT[[u]]$coefficients$fixed))))
      })  
    })),2,function(j){
      rep(j,test_indxJ*L)
    })
    
    tsDFACT<-sapply(seq_len(ncol(tsREF)),function(u){
      tsREF[,u]+as.numeric(spatFIT[[u]]$coefficients$fixed)
    })
    
    tscMAT<-do.call(cbind,lapply(seq_len(sp_nbasis),function(v){
      t(sapply(split(tsp_score[,v]-tsDFACT[,v], rep(1:sum(test_indxJ),each=L)),function(u){
        c(mean(u),mean(u^2))
      }))
    }))
  } else{
    tscMAT<-NULL
  }
  list(
    "n" = n,
    "p" = p,
    "funARG" = funARG,
    "L" = nrow(spat_indx),
    "K" = sp_nbasis,
    "EigenF"=sp_basis,
    "XiS" = sp_score,
    "fitOBJ_SRF"=spatFIT,
    "spat_corr_type" = corrM,
    "spatCOV" = lapply(seq_len(length(spatFIT)),function(u){
      as.numeric(coef(spatFIT[[u]]$modelStruct$corStruct,unconstrained = F))
    }),
    "REvar" = sapply(seq_len(length(spatFIT)),function(u){
      as.numeric(simplify2array(getVarCov(spatFIT[[u]],type = "random.effects")))
    }),
    "SigmaK" = sapply(seq_len(length(spatFIT)),function(u){
      as.numeric(spatFIT[[u]]$sigma)
    }),
    "muK" = sapply(seq_len(length(spatFIT)),function(u){
      as.numeric(spatFIT[[u]]$coefficients$fixed)
    }),
    "Sigma2" = as.numeric(tapply(us_var-diag(smth_cov),rep(1:2,each=length(funARG)),mean)),
    "trainRD" = scMAT,
    "testRD" = tscMAT
  )
}




#' Give the projection on the reduced dimension obtained from RSIMVFD decomposition
#' 
#' @param fitOBJ an object obtained by the dredSIMVFD
#' @param testFDATA is a list of functional data from the test for which we want 
#' to obtain the projection onto the reduced dimension
#' @param test_indxJ is similar to indxJ, however, it corresponds to the test data provided
#' @param test_spat_index is a L by 2 matrix contains spatial locations for the test data
#' contains information on mi
#' @return a list with following objects:
#' \itemize{
#'  \item n : number of subjects in the test data
#'  \item test_SRF : extracted spatial random fields
#'  \item test_factor_subtract : predicted subject-specific mean
#'  \item testRD : derived features from extracted SRF
#' }
#' @example examples/example_pred_dcompSiMVFD.R
#' @export
pred_dcompSiMVFD<-function(fitOBJ,testFDATA,test_indxJ,test_spat_indx){
  p<-fitOBJ$p
  L<-fitOBJ$L
  sp_basis<-fitOBJ$EigenF
  sp_nbasis<-fitOBJ$K
  funARG<-fitOBJ$funARG
  spatFIT<-fitOBJ$fitOBJ_SRF
  ## Construction of marginal covariance
  corrM<-fitOBJ$spat_corr_type
  disD<-as.matrix(dist(test_spat_indx,diag = TRUE,upper = TRUE))
  vcovSP<-lapply(spatFIT, function(u){
    cor_par<-as.numeric(coef(u$modelStruct$corStruct,unconstrained = F))
    apply(disD,2,function(v){sapply(v,covSpatialVal,d=cor_par[1],nugget=cor_par[2],sill_cov=u$sigma^2,type=corrM)})
  })
  
  tn<-length(test_indxJ)
  tsubID<-rep(1:tn,test_indxJ*L)
  tsp_score<-(do.call(cbind,testFDATA)%*%do.call(rbind,sp_basis))*(1/length(funARG))
  
  tspDAT<-split.data.frame(tsp_score,tsubID)
  
  tsREF<-apply(t(sapply(seq_len(length(tspDAT)), function(v){
    OmL<-rep(1,test_indxJ[v]*L)
    ImM<-diag(test_indxJ[v])
    sapply(seq_len(length(spatFIT)), function(u){
      tauK<-(as.numeric(attr(corMatrix(spatFIT[[u]]$modelStruct[[1]])[[1]],"stdDev")*spatFIT[[u]]$sigma))^2
      covM<-(tauK*outer(OmL,OmL))+kronecker(ImM,vcovSP[[u]])
      as.numeric(matrix(tauK*OmL,nrow=1)%*%(solve(covM)%*%matrix(tspDAT[[v]][,u]-as.numeric(spatFIT[[u]]$coefficients$fixed))))
    })  
  })),2,function(j){
    rep(j,test_indxJ*L)
  })
  
  tsDFACT<-sapply(seq_len(ncol(tsREF)),function(u){
    tsREF[,u]+as.numeric(spatFIT[[u]]$coefficients$fixed)
  })
  
  tscMAT<-do.call(cbind,lapply(seq_len(sp_nbasis),function(v){
    t(sapply(split(tsp_score[,v]-tsDFACT[,v], rep(1:sum(test_indxJ),each=L)),function(u){
      c(mean(u),mean(u^2))
    }))
  }))
  list(
    "n" = length(test_indxJ),
    "test_SRF" = tsREF,
    "test_factor_subtract" = tsDFACT,
    "testRD" = tscMAT
  )
}


#' Train a classifier based multivariate FPCA (separate basis) of spatially 
#' indexed multivariate functional data
#' 
#' @param fDATA is a list of functional data, length is equal to the number of 
#' function
#' @param indxJ is a vector of length n that has elements represent m_i
#' @param spat_indx is a L by 2 matrix contains spatial locations
#' @param funARG is vector of sampled points where functional data were observed. 
#' This can be generalized for sparse functional data. However this version 
#' is only for
#' dense functional data. Also all the functions are observed at identical 
#' sampled point.
#' @param PVE is a vector of length 2 for percentage of variation explained. 
#' The first is used for fpca of each functional variable. The second is for 
#' eigen decomposition of scores from each of the 
#' @param corrM is a character specifies the name of spatial correlation
#' @param center is a logical vector of length 1 indicates whether the 
#' functional data need centering by the sample mean
#' @param testFDATA is a list of functional data from the test for which we want 
#' to obtain the projection onto the reduced dimension
#' @param test_indxJ is similar to indxJ, however, it corresponds to the test 
#' data provided
#' @keywords internal
#' @noRd
#' @export
MFPCA_HG_dcompSiMVFD<-function(fDATA,indxJ,spat_indx,funARG,PVE,corrM,
                               center=FALSE,nknots=5,order_penalty=2,
                               basisT=c("ps","ps"),sameEF="FALSE",
                               testFDATA=NULL,test_indxJ=NULL, ...){
  if(!is.null(testFDATA)&is.null(test_indxJ))
    stop("For prediction test_indxJ is required")
  n<-length(indxJ) 
  p<-length(fDATA)
  L<-nrow(spat_indx)
  if(center){
    fDATA<-lapply(fDATA, function(u){
      scale(u,center = TRUE,scale = FALSE)
    })
  }
  
  # Dimension reduction
  if(sameEF){
    # multivariate functional principal component analysis
    # Subject-patch specific basis
    smth_cov<-smth_covar(funDATA=fDATA,
                         funARG=funARG,
                         pve=PVE[1],
                         nKnot=nknots,
                         basisT=basisT,
                         tps_m=order_penalty)
    
    # Spectral decomposition of smooth covariance
    Wmat<-diag(rep(1/sqrt(length(funARG)),2*length(funARG)))
    IWmat<-diag(rep(sqrt(length(funARG)),length(funARG)))
    spd_smth_covar<-eigen((Wmat%*%(smth_cov%*%Wmat)))
    plamd<-spd_smth_covar$values[spd_smth_covar$values>0]
    sp_nbasis<-which((cumsum(plamd)/sum(plamd))>PVE[2])[1]
    sp_basis<-lapply(split.data.frame((spd_smth_covar$vectors[,1:sp_nbasis]),rep(1:p,each=length(funARG))),function(w){IWmat%*%w})
    MfpcaSCR<-(do.call(cbind,fDATA)%*%do.call(rbind,sp_basis))*(1/length(funARG))
    
    ###### Estimation of Zeta
    subID<-rep(1:n,indxJ*L)
    subMEAN<-lapply(fDATA, function(u){
      apply(u,2,function(v){
        as.numeric(tapply(v,subID,mean))
      })
    })
    Zeta<-Reduce(`+`,lapply(seq_len(length(fDATA)), function(u){
      (subMEAN[[u]]%*%sp_basis[[u]])*(1/length(funARG))
    }))
    
    sp_score<-as.data.frame(MfpcaSCR-apply(Zeta,2,function(w){rep(w,indxJ*L)}))
    
    sub_basisF<-sp_basis
    # Test data subject-patch score
    if(!is.null(testFDATA)){
      if(center){
        tFDATA<-lapply(seq_len(p), function(u){
          testFDATA[[u]]-outer(rep(1,nrow(testFDATA[[u]])),colMeans(fDATA[[u]]))
        })
      } else{
        tFDATA<-testFDATA
      }
      
      ###### Estimation of Zeta for test data
      tn<-length(test_indxJ)
      tsubID<-rep(1:tn,test_indxJ*L)
      tsubMEAN<-lapply(testFDATA, function(u){
        apply(u,2,function(v){
          as.numeric(tapply(v,tsubID,mean))
        })
      })
      
      Tscr<-Reduce(`+`,lapply(seq_len(p), function(u){
        (tFDATA[[u]]%*%sp_basis[[u]])*(1/length(funARG))
      }))
      
      tZeta<-Reduce(`+`,lapply(seq_len(p), function(u){
        (tsubMEAN[[u]]%*%sp_basis[[u]])*(1/length(funARG))
      }))
      
      tsp_score<-Tscr-apply(tZeta,2,function(w){rep(w,test_indxJ*L)})
    }
  } else{
    # Subject effect: basis and scores
    subID<-rep(1:n,indxJ*L)
    subMEAN<-lapply(fDATA, function(u){
      apply(u,2,function(v){
        as.numeric(tapply(v,subID,mean))
      })
    })
    
    # Univariate fpca of subject-specific means
    umFIT<-lapply(seq_len(length(subMEAN)), function(u){
      fpca.sc(Y=subMEAN[[u]],argvals = funARG,pve = PVE[1])
    })
    
    # estimation of univariate scores
    uzeta_vec<-do.call(cbind,lapply(seq_len(length(subMEAN)),function(u){
      (1/length(funARG))*((subMEAN[[u]]-outer(rep(1,n),umFIT[[u]]$mu))%*%umFIT[[u]]$efunctions)
    }))
    
    spd_uzeta<-eigen(cov(uzeta_vec))
    
    npc_uz<-which((cumsum(spd_uzeta$values)/sum(spd_uzeta$values))>PVE[2])[1]
    
    CMat<-spd_uzeta$vectors[,1:npc_uz]
    
    # subject-effect eigen-function
    sub_basisF<-split.data.frame(as.matrix(Matrix::bdiag(lapply(umFIT,function(u){u$efunctions}))%*%CMat),
                                 rep(1:p,each=length(funARG)))
    
    # Subject-specific scores
    Zeta<-Reduce(`+`,lapply(1:p,function(u){
      (1/length(funARG))*((subMEAN[[u]]-outer(rep(1,n),umFIT[[u]]$mu))%*%sub_basisF[[u]])
    }))
    
    # demeaned data
    
    dfDATA<-lapply(1:p, function(u){
      fDATA[[u]]-outer(rep(1,length(subID)),umFIT[[u]]$mu)-(Zeta[subID,]%*%t(sub_basisF[[u]]))
    })
    
    # Subject-patch specific basis
    smth_cov<-smth_covar(funDATA=dfDATA,
                         funARG=funARG,
                         pve=PVE[1],
                         nKnot=nknots,
                         basisT=basisT,
                         tps_m=order_penalty)
    
    # Spectral decomposition of smooth covariance
    Wmat<-diag(rep(1/sqrt(length(funARG)),2*length(funARG)))
    IWmat<-diag(rep(sqrt(length(funARG)),length(funARG)))
    spd_smth_covar<-eigen((Wmat%*%(smth_cov%*%Wmat)))
    plamd<-spd_smth_covar$values[spd_smth_covar$values>0]
    sp_nbasis<-which((cumsum(plamd)/sum(plamd))>PVE[3])[1]
    sp_basis<-lapply(split.data.frame((spd_smth_covar$vectors[,1:sp_nbasis]),rep(1:p,each=length(funARG))),function(w){IWmat%*%w})
    sp_score<-(do.call(cbind,dfDATA)%*%do.call(rbind,sp_basis))*(1/length(funARG))
    
    # Test data subject-patch specific scores
    if(!is.null(testFDATA)){
      if(center){
        tFDATA<-lapply(seq_len(p), function(u){
          testFDATA[[u]]-outer(rep(1,nrow(testFDATA[[u]])),colMeans(fDATA[[u]]))
        })
      } else{
        tFDATA<-testFDATA
      }
      ## Estimation of Zeta for test data
      tn<-length(test_indxJ)
      tsubID<-rep(1:tn,test_indxJ*L)
      tsubMEAN<-lapply(testFDATA, function(u){
        apply(u,2,function(v){
          as.numeric(tapply(v,tsubID,mean))
        })
      })
      
      # Subject-specific scores for test data
      tZeta<-Reduce(`+`,lapply(1:p,function(u){
        (1/length(funARG))*((tsubMEAN[[u]]-outer(rep(1,tn),umFIT[[u]]$mu))%*%sub_basisF[[u]])
      }))
      
      # demeaned data for test set
      
      tdfDATA<-lapply(1:p, function(u){
        testFDATA[[u]]-outer(rep(1,length(tsubID)),umFIT[[u]]$mu)-(tZeta[tsubID,]%*%t(sub_basisF[[u]]))
      })
      
      # subject-patch specific scores for test data
      tsp_score<-(do.call(cbind,tdfDATA)%*%do.call(rbind,sp_basis))*(1/length(funARG))
    }
  }
  
  # Predictors for training data
  mfpcaSCR<-as.data.frame(sp_score)
  mfpcaSCR$Subject<-rep(1:n,indxJ*L)
  mfpcaSCR$IndexJ<-as.factor(rep(1:sum(indxJ),each=L))
  mfpcaSCR$Row<-rep(spat_indx[,1],times=sum(indxJ))
  mfpcaSCR$Col<-rep(spat_indx[,2],times=sum(indxJ))
  
  spatFIT<-lapply(seq_len(ncol(sp_score)), function(u){
    sdata<-as.data.frame(cbind("Y"=mfpcaSCR[,u],mfpcaSCR[,-c(1:ncol(sp_score))]))
    lme(fixed=Y~1,data=sdata,random=~1|IndexJ,correlation = corSpatial(form=~Row+Col|IndexJ,nugget = TRUE,type = corrM), ...)
  })
  
  rd_effect<-sapply(seq_len(length(spatFIT)), function(u){
    vcv_wme<-simplify2array(getVarCov(spatFIT[[u]],type = "marginal"))[1:L,1:L,1]
    gvcv<-MASS::ginv(vcv_wme)
    sapply(split((sp_score[,u]-as.numeric(spatFIT[[u]]$coefficients$fixed)),rep(1:sum(indxJ),each=L)),function(w){
      sum(gvcv%*%matrix(w,ncol=1))*as.numeric(simplify2array(getVarCov(spatFIT[[u]],type = "random.effects")))
    })
  })
  
  scMAT<-do.call(cbind,lapply(seq_len(length(spatFIT)), function(u){
    cbind(rd_effect[,u],log(tapply(mfpcaSCR[,u]-as.numeric(spatFIT[[u]]$coefficients$fixed)-rep(rd_effect[,u],each=L),mfpcaSCR$IndexJ,var)))
  }))
  
  
  mfDATA<-if(center){
    lapply(fDATA,colMeans)
  } else{
    NULL
  }
  
  if(!is.null(testFDATA)){
    ## Test random effects 
    trd_effect<-sapply(seq_len(length(spatFIT)), function(u){
      vcv_wme<-simplify2array(getVarCov(spatFIT[[u]],type = "marginal"))[1:L,1:L,1]
      gvcv<-MASS::ginv(vcv_wme)
      sapply(split((tsp_score[,u]-as.numeric(spatFIT[[u]]$coefficients$fixed)),rep(1:sum(test_indxJ),each=L)),function(w){
        sum(gvcv%*%matrix(w,ncol=1))*as.numeric(simplify2array(getVarCov(spatFIT[[u]],type = "random.effects")))
      })
    })
    
    
    tscMAT<-do.call(cbind,lapply(seq_len(length(spatFIT)), function(u){
      cbind(trd_effect[,u],log(tapply(tsp_score[,u]-as.numeric(spatFIT[[u]]$coefficients$fixed)-rep(trd_effect[,u],each=L),rep(1:sum(test_indxJ),each=L),var)))
    }))
    
    
    retOBJ<-list(
      "n" = n,
      "p" = p,
      "funARG" = funARG,
      "L" = nrow(spat_indx),
      "mFDATA"=mfDATA,
      "SubEigenF"=sub_basisF,
      "Zeta" = Zeta,
      "EigenF"=sp_basis,
      "XiS" = sp_score,
      "spatCOV" = lapply(seq_len(length(spatFIT)),function(u){
        simplify2array(getVarCov(spatFIT[[u]],type = "marginal"))[1:L,1:L,1]
      }),
      "REvar" = sapply(seq_len(length(spatFIT)),function(u){
        as.numeric(simplify2array(getVarCov(spatFIT[[u]],type = "random.effects")))
      }),
      "trainRD" = scMAT,
      "testRD" = tscMAT,
      "TZeta" = tZeta
    )
  } else{
    retOBJ<-list(
      "n" = n,
      "p" = p,
      "funARG" = funARG,
      "L" = nrow(spat_indx),
      "mFDATA"=mfDATA,
      "SubEigenF"=sub_basisF,
      "Zeta" = Zeta,
      "EigenF"=sp_basis,
      "XiS" = sp_score,
      "spatCOV" = lapply(seq_len(length(spatFIT)),function(u){
        simplify2array(getVarCov(spatFIT[[u]],type = "marginal"))[1:L,1:L,1]
      }),
      "REvar" = sapply(seq_len(length(spatFIT)),function(u){
        as.numeric(simplify2array(getVarCov(spatFIT[[u]],type = "random.effects")))
      }),
      "trainRD" = scMAT,
      "testRD" = NULL,
      "TZeta" = NULL
    )
  }
  retOBJ
}