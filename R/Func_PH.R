#==========================================================================#
# Cox proportional hazards model (PH) -- based on: Profile likelihood
#==========================================================================#

#==== Function for fitting Cox proportional hazards model with partial likelihood ====#
#' @title Cox proportional hazards model based on the profiled likelihood estimation procedure
#'
#' @description Fit Cox proportional hazards model based on the profiled likelihood estimation procedure.
#'
#' @aliases PH.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates
#' @param Var a logical value. The default setup is \code{TRUE}, indicating that the standard errors of the estimated regression coefficients will be calculated.
#' @param Var.Robust a logical value. The default setup is \code{TRUE}, indicating that we assume that the model is correctly specified, so that there is no need to use the robust estimator of standard errors.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{res}.
#'
#' @export PH.fit
PH.fit <- function(yobs,delta,X,Var=TRUE,Var.Robust=FALSE,eps=1e-6,maxit=5e4){
  # need: survival

  ### Preparations
  N <- length(yobs)
  pbet <- ncol(X)

  ### calculate initial values for bet
  bet.init <- rep(0,pbet)

  ### calculate MLE for beta
  numit <- 1
  bet.old <- bet.init
  repeat{
    InfoMScore <- PH.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
    Score <- InfoMScore$Score
    InfoM <- InfoMScore$InfoM
    dev <- MASS::ginv(InfoM)%*%Score # solve(InfoM,Score)
    bet <- bet.old + as.vector(dev)
    if( max(abs(bet-bet.old))>eps & numit<maxit ){
      bet.old <- bet
      numit <- numit + 1
    }else{
      break
    }
  }
  convergence <- (numit<maxit)

  ### calculate SEs or not (and tidy them)
  if(Var==TRUE){

    ## calculate SEs for bet using explicit formula !!!
    InfoMScore <- PH.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
    InfoM <- InfoMScore$InfoM
    Score <- InfoMScore$Score
    InfoM.inv <- solve(InfoM)
    if(Var.Robust==TRUE){
      Score.Influs <- PH.Influence.EE(yobs=yobs,delta=delta,X=X,bet=bet)
      Score.Influs.Cov <- t(Score.Influs)%*%Score.Influs/N
      VarCov <- InfoM.inv %*% Score.Influs.Cov %*% InfoM.inv
    }else{
      VarCov <- InfoM.inv
    }
    bet.se <- sqrt(diag(VarCov)/N)

    ### tidy the results: inference
    zvalue.bet <- bet/bet.se
    pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
    res <- data.frame(Est=bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                      row.names=colnames(X))

  }else{

    # tidy the results directly
    VarCov <- NULL
    res <- data.frame(Est=bet,row.names=colnames(X))

  }

  ### output
  out <- list(
    info = list(
      convergence = convergence,
      bet.init = bet.init
    ),
    res=res,
    Score=Score,
    InfoM=InfoM,
    VarCov=VarCov
  )
  return(out)

}


#==== Profiled log-likelihood function in Cox proportional hazards model ====#
#' @title Profiled log-likelihood function in Cox proportional hazards model
#'
#' @description Calculate profiled log-likelihood function in Cox proportional hazards model.
#'
#' @aliases PH.Beta.LogLik
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates
#'
#' @export PH.Beta.LogLik
PH.Beta.LogLik <- function(bet,yobs,delta,X){
  # for maximization

  ## prepare
  N <- length(yobs)
  Xbet <- as.vector(X %*% bet)
  log.SS0 <- sapply(yobs,function(Yi){log(mean(exp(Xbet)*(yobs>=Yi)))})

  ## calculate the partial likelihood function for beta (log form)
  val.log <- sum((Xbet-log.SS0)*delta)

  # output
  return(val.log)

}


#==== Obtain score vector and information matrix ====#
#' @title Score vector and information matrix in Cox proportional hazards model
#'
#' @description Calculate score vector and information matrix in Cox proportional hazards model.
#'
#' @aliases PH.ScoreInfoM
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param IsScore whether the score vector is calculated or not.
#' @param IsInfoM whether the information matrix is calculated or not.
#'
#' @export PH.ScoreInfoM
PH.ScoreInfoM <-  function(bet,yobs,delta,X,IsScore=FALSE,IsInfoM=TRUE){
  # Score is the first  derivative of [positive] log-likelihood
  # InfoM is the second derivative of [negative] log-likelihood

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])})/N # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))/N # [S(1)(t,bet)] at pre-specified yobs and beta
  out <- list()

  ## prepare information matrix
  if(IsInfoM==TRUE){
    SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
      yobsGYi <- yobs>=Yi
      as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
    }))/N
    I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
    I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
    InfoM <- I1-I2
    out <- c(out,list(InfoM=InfoM))
  }

  ## prepare score vector
  if(IsScore==TRUE){
    U <- (X - SS1/SS0)*delta
    Score <- apply(U,2,mean)
    out <- c(out,list(Score=Score))
  }


  ## output
  return(out)

}


#==== Obtain individual level score vector ====#
#' @title Individual level score vector in Cox proportional hazards model
#'
#' @description Calculate the individual level score vector in Cox proportional hazards model.
#'
#' @aliases PH.Score.Individual
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param Iscov whether convert the individual level form into a covariance form.
#'
#' @export PH.Score.Individual
PH.Score.Individual <-  function(bet,yobs,delta,X,Iscov=TRUE){
  # Score is the first  derivative of [positive] log-likelihood

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])})/N # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)}))/N # [S(1)(t,bet)] at pre-specified yobs and beta
  U <- (X - SS1/SS0)*delta
  if(Iscov==TRUE){
    out <- t(U)%*%U/N
  }else{
    out <- U
  }

  ## output
  return(out)

}


#==== Nonparametric component ====#
#' @title Nonparametric component in Cox proportional hazards model
#'
#' @description  Nonparametric component for the Cox proportional hazards model at specified time points.
#'
#' @aliases PH.Lam
#'
#' @param tm The time points that the nonparametric baseline function will be estimated at.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#' @param type if \code{type="right"}, right-continuous results will be calculated, and if \code{type="left"}, left-continuous results will be calculated.
#'
#' @export PH.Lam
PH.Lam <- function(tm,yobs,delta,X,bet=NULL,type="right"){

  # prepare
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  N <- length(yobs)
  expXbet <- as.vector(exp(X %*% bet))
  SS0 <- sapply(yobs,function(Yi){mean(expXbet*(yobs>=Yi))}) # [S(0)(t,bet)] at pre-specified yobs and beta

  # calculate the Lam(t) at specified time points tm
  if(type=="right"){
    Lam <- sapply(tm,function(tmj){
      sum( (yobs<=tmj)[delta==1]/SS0[delta==1] ) / N
    })
  }else if(type=="left"){
    Lam <- sapply(tm,function(tmj){
      sum( (yobs< tmj)[delta==1]/SS0[delta==1] ) / N
    })
  }


  # output
  return(Lam)

}


#==== The influences for proportional hazards model ====#
#' @title  Influence functions for the estimator of regression coefficients in the Cox proportional hazards model
#'
#' @description Calculate the influence functions for the estimator of regression coefficients in the Cox proportional hazards model.
#'
#' @aliases PH.Influence
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#'
#' @export PH.Influence
PH.Influence <- function(yobs,delta,X,bet=NULL){

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## calculate main part (except an inverse of information matrix)
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet

  ## prepare information matrix
  SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
    yobsGYi <- yobs>=Yi
    as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
  }))
  I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
  I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
  InfoM <- I1-I2

  ## final influence functions (individual level) and output
  Influs <- U%*%solve(InfoM)
  return(Influs)

}


#==== The influences for score of proportional hazards model ====#
#' @title  Influence functions for the estimator of the score equation in the Cox proportional hazards model
#'
#' @description Calculate the influence functions for the estimator of the score equation in the Cox proportional hazards model.
#'
#' @aliases PH.Influence.EE
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#'
#' @export PH.Influence.EE
PH.Influence.EE <- function(yobs,delta,X,bet=NULL){
  # The incluence function for estimating equations

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## calculate main part
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet

  ## final influence functions (individual level) and output
  return(U)

}


#==== The influences for nonparametric part of proportional hazards model ====#
#' @title  Influence functions for the estimator of the nonparametric part in the Cox proportional hazards model
#'
#' @description Calculate the influence functions for the estimator of the nonparametric part in the Cox proportional hazards model.
#'
#' @aliases PH.Influence.Lam
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet unknown parameters corresponding to the model.
#'
#' @export PH.Influence.Lam
PH.Influence.Lam <- function(tm,yobs,delta,X,bet=NULL){

  # tm <- c(0.2,0.4,0.6,0.8,1)

  ## prepare elements
  N <- length(yobs)
  pbet <- ncol(X)
  if(is.null(bet)){
    bet <- PH.fit(yobs,delta,X)$res[,1]
  }
  Xbet <- as.vector(X %*% bet)
  expXbet <- exp(Xbet)
  XexpXbet <- X*expXbet
  SS0 <- sapply(yobs,function(Yi){sum(expXbet[yobs>=Yi])}) # [S(0)(t,bet)] at pre-specified yobs and beta
  SS1 <- do.call(rbind,lapply(yobs,function(Yi){apply(XexpXbet[yobs>=Yi,,drop=F],2,sum)})) # [S(1)(t,bet)] at pre-specified yobs and beta

  ## prepare information matrix
  SS2.vec <- do.call(rbind,lapply(yobs,function(Yi){
    yobsGYi <- yobs>=Yi
    as.vector( t(X[yobsGYi,,drop=F])%*%(XexpXbet[yobsGYi,,drop=F]) )
  }))
  I1 <- matrix(apply(SS2.vec*delta/SS0,2,mean),nrow=pbet,ncol=pbet)
  I2 <- t(SS1)%*%(SS1*delta/SS0^2)/N
  InfoM <- I1-I2

  ## calculate main part 1 (same as bet)
  UU1 <- (X - SS1/SS0)*delta
  UU2 <- X * sapply(yobs,function(Yi){sum(delta*(yobs<=Yi)/SS0)})
  UU3 <- do.call(rbind,lapply(yobs,function(Yi){apply(delta*(yobs<=Yi)*SS1/SS0^2,2,sum)}))
  U <- UU1-(UU2-UU3)*expXbet
  Influs.bet <- U%*%solve(InfoM)

  ## calculate main part 2 (with tm vary)
  LL1.tm <- do.call(cbind,lapply(tm,function(tmi){(yobs<=tmi)*delta/(SS0/N)}))
  LL2.tm <- do.call(cbind,lapply(tm,function(tmi){
    sapply(yobs,function(Yj){sum(delta*(yobs<=min(tmi,Yj))/(SS0^2/N))})
  }))*expXbet
  LL3.pre <- do.call(cbind,lapply(tm,function(tmi){apply(SS1*delta*(yobs<=tmi)/SS0^2,2,sum)}))
  LL3.tm <- Influs.bet %*% LL3.pre

  ## influence function for Lam (individual level) and output
  Influs.Lam <- LL1.tm - LL2.tm - LL3.tm

  ## final influence functions
  return(Influs.Lam)

}


