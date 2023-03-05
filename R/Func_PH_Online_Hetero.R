


#==== Online estimation in Cox proportional hazards model with dynamic covariate effects (one batch, Renewable, ALasso) ====#
#' @title Online estimation in Cox proportional hazards model with dynamic covariate effects
#'
#' @description Conduct Online estimation in Cox proportional hazards model with dynamic covariate effects.
#'
#' @aliases PH.Renew.Hetero.ALasso.Batch
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param initial a logical value. The default specification is \code{TRUE}, indicating that the current data batch is the first batch and no historical summary statistics is available.
#' @param prev indicates the historical summary data. The default specification is \code{NULL} (with \code{initial=TRUE}). Otherwise, it is a list with the following five elements:
#'   \code{bet} is the historical estimated result of regression coefficients;
#'   \code{bet.batch} is the historical estimated result of regression coefficients using only the previous dataset;
#'   \code{InfoM} is the historical estimated result of the information matrix;
#'   \code{N} is the sample of historical raw data.
#' @param boots basic setup of the bootstrap procedure.
#' @param tunpara.num the number of tuning parameter values that will be fitted at.
#' @param tunpara.min.ratio Smallest value for tuning parameter, as a fraction of tunpara.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero).
#' @param tunpara A user supplied \code{tunpara} sequence. Typical usage is to have the program compute its own \code{tunpara} sequence based on \code{tunpara.num} and \code{tunpara.min.ratio}. Supplying a value of \code{tunpara} overrides this.
#' @param evatype the information criterion that will be used in the selection of optimal tuning parameter.
#' @param learnrate the learning rate in solving the optimization problem.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param threshold a threshold value. The deafult value is \code{1e-5}. If the absolute value of the each component of coefficients is smaller than \code{threshold}, it will be converted into zero.
#' @param trace a logical value. The default specification is \code{TRUE}, indicating that several useful information about the current process of fitting will be displayed.
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{resbet}.
#'
#' @examples
#' ### ==== An example for fitting the online Cox proportional hazards model (Heterogeneity exist) ==== ###
#'
#' ## ---- basic setup of the simulated dataset ---- ##
#' pbet <- 3
#' B.oT <- 5
#' Nb.oT <- 400
#' bettau.jump <- matrix(
#'   c(c( 1,-2, 0),  # bet0
#'     c(-1, 0, 2),  # bet1
#'     c( 1, 0, 0)), # bet2
#'   ncol=pbet)
#' Tnum <- nrow(bettau.jump)
#' bettau.all  <- matrix(0,nrow=B.oT*Tnum,ncol=pbet)
#' bettau.all[(0:(Tnum-1))*B.oT+1,]  <- bettau.jump
#' bet.all <- apply(bettau.all,2,cumsum)
#'
#' ## ---- fit the  Cox proportional hazards model in an online manner ---- ##
#' # prepare basic elements
#' B  <- nrow(bet.all)
#' Nb <- Nb.oT
#' N <- B*Nb
#' Res.Online <- array(0,dim=c(B,pbet,2),
#'                     dimnames=list(paste("batch",1:B,sep=""),paste("bet",1:pbet,sep=""),c("EST","SE")))
#' Res.Online_Hetero <- list(
#'   res = array(0,dim=c(B,pbet,2),
#'               dimnames=list(paste("batch",1:B,sep=""),paste("bet",1:pbet,sep=""),c("EST","SE"))),
#'   restau = array(0,dim=c(B-1,pbet,2),
#'                  dimnames=list(paste("batch",2:B,sep=""),paste("taub",1:pbet,sep=""),c("EST","SE")))
#' )
#' # Online estimation procedures
#' for(b in 1:B){ # b <- 1
#'   cat("-------------------------\nCurrent batch number: ",b,"\n")
#'
#'   ## Generate the current data batch
#'   sdata.batch  <- sdata.PH(N=Nb,bet=bet.all[b,])
#'   yobs  <- sdata.batch$sdata$yobs
#'   delta <- sdata.batch$sdata$delta
#'   X     <- as.matrix(sdata.batch$sdata$X)
#'
#'   ### ____ Online PH (homogeneous) ____ ###
#'   ## preparation: the previous elements
#'   if(b==1){ prevb <- NULL }else{
#'     prevb <- c(
#'       fitb$fit
#'     )
#'   }
#'   ## updating
#'   fitb <- PH.Renew.Batch(
#'     yobs=yobs,delta=delta,X=X,
#'     initial=(b==1),prev=prevb)
#'   Res.Online[b,,1] <- fitb$res[,1]
#'   Res.Online[b,,2] <- fitb$res[,2]
#'
#'   ### ____ Online PH (heterogeneous) ____ ###
#'   ## preparation: the previous elements
#'   if(b==1){ prevb.Hetero <- NULL }else{
#'     prevb.Hetero <- c(
#'       fitb.Hetero$fit,
#'       fit.boots=list(fitb.Hetero$fit.boots)
#'     )
#'   }
#'   ## updating
#'   fitb.Hetero <- PH.Renew.Hetero.ALasso.Batch(
#'     yobs=yobs,delta=delta,X=X,
#'     initial=(b==1),prev=prevb.Hetero,
#'     boots=list(do=TRUE,nboot=100,trace=TRUE))
#'   Res.Online_Hetero$res[b,,1] <- fitb.Hetero$resbet[,1] # Res.Online_Hetero$res[,,1]
#'   Res.Online_Hetero$res[b,,2] <- fitb.Hetero$resbet[,2]
#'   if(b>1){
#'     Res.Online_Hetero$restau[b-1,,1] <- fitb.Hetero$restau[,1] # Res.Online_Hetero$restau[,,1]
#'     Res.Online_Hetero$restau[b-1,,2] <- fitb.Hetero$restau[,2]
#'   }
#'
#' }
#' # present the fitted results
#' print(Res.Online)
#' print(Res.Online_Hetero$res)
#' print(Res.Online_Hetero$restau)
#'
#' @export PH.Renew.Hetero.ALasso.Batch
PH.Renew.Hetero.ALasso.Batch <- function(
    yobs,delta,X,initial=TRUE,prev=NULL,
    boots=list(do=FALSE,nboot=100,trace=TRUE),
    tunpara.num=50,
    tunpara.min.ratio=0.0001,
    tunpara=NULL,
    evatype=c("AIC","BIC")[2],
    learnrate=c(0.4),
    maxit=1e3,eps=1e-6,
    threshold=5e-5,
    trace=TRUE
){

  N <- length(yobs)

  ### fit the basic model
  resss <- PH.Renew.Hetero.ALasso.Batch.fit(
    yobs=yobs,delta=delta,X=X,initial=initial,prev=prev,
    tunpara.num=tunpara.num,tunpara.min.ratio=tunpara.min.ratio,tunpara=tunpara,
    evatype=evatype,learnrate=learnrate,maxit=maxit,eps=eps,threshold=threshold,
    trace=trace
  )
  fit <- resss$fit

  ### calculate SEs
  if(boots$do==FALSE){


    if(initial==TRUE){


      InfoM <- fit$InfoM
      V <- solve(InfoM)
      bet.se <- sqrt(diag(V)/N)
      tau.se <- NA


    }else{

      InfoM <- fit$InfoM
      InfoM.prev <- prev$InfoM
      N.prev <- prev$N
      N <-  length(yobs)
      pbet <- ncol(X)
      InfoM.Big <- array(0,dim=c(2*pbet,2*pbet))
      InfoM.Big[1:pbet,1:pbet] <- InfoM
      InfoM.Big[(pbet+1):(2*pbet),1:pbet] <- InfoM.Big[1:pbet,(pbet+1):(2*pbet)] <- -(N.prev/(N.prev+N))*InfoM.prev
      InfoM.Big[(pbet+1):(2*pbet),(pbet+1):(2*pbet)] <- (N.prev/(N.prev+N))*InfoM.prev

      bet <- fit$bet
      tau <- resss$tau
      bettau <- c(bet,tau)
      tunpara.opt <- resss$tunpara.opt
      idx.nozero <- c(rep(TRUE,pbet),tau!=0)
      idx.zero <- !idx.nozero

      InfoM11 <- InfoM.Big[idx.nozero,idx.nozero]
      InfoM12 <- InfoM.Big[idx.nozero,!idx.nozero]
      InfoM22 <- InfoM.Big[!idx.nozero,!idx.nozero]
      InfoM11.inv <- MASS::ginv(InfoM11)
      if(any(tau==0)){
        tunpara.opt.vector <- tunpara.opt#c(rep(0,pbet),rep(tunpara.opt,sum(idx.nozero)-pbet))
        InfoM11.tilde.inv <- MASS::ginv(InfoM11 + tunpara.opt.vector*diag(1/bettau[idx.nozero]^2))
        E.inv <- MASS::ginv(InfoM22 - t(InfoM12)%*%InfoM11.inv%*%InfoM12)
        covbet1 <- InfoM11.inv + (InfoM11.inv-InfoM11.tilde.inv)%*%InfoM12%*%E.inv%*%t(InfoM12)%*%(InfoM11.inv-InfoM11.tilde.inv)
      }else{
        covbet1 <- InfoM11.inv
      }
      bettau.se <- rep(NA,2*pbet)
      bettau.se[idx.nozero] <- sqrt(diag(covbet1)/(N+N.prev))
      bet.se <- bettau.se[1:pbet]
      tau.se <- bettau.se[-c(1:pbet)]

    }

    # can we use explicit formula
    fit.boots <- NULL
    # bet.se <- NA
    # tau.se <- NA

  }else if(boots$do==TRUE){

    ## calculate SEs for bet using bootstrap (direct or weighted)
    fit.boots <- rep(list(list()),boots$nboot)
    tau.boots <- array(0,dim=c(boots$nboot,ncol(X)))
    if(boots$trace==TRUE){cat("Do bootstrap:\n")}
    iboot <- 1
    while(iboot<=boots$nboot){

      cytry <- try({
        # bootstrap idx of samples
        idx.iboot <- sample(1:N, replace = TRUE)
        # idx.iboot <- c(sample((1:N)[delta==1], replace = TRUE),sample((1:N)[delta==0], replace = TRUE))
        # solve the current estimator
        resss.iboot <- PH.Renew.Hetero.ALasso.Batch.fit(
          yobs=yobs[idx.iboot],delta=delta[idx.iboot],X=X[idx.iboot,,drop=F],initial=initial,prev=prev$fit.boots[[iboot]],
          tunpara.num=NULL,tunpara.min.ratio=NULL,tunpara=resss$tunpara.opt,
          evatype=evatype,learnrate=learnrate,maxit=maxit,eps=eps,threshold=threshold,
          trace=FALSE)
      }, silent = T) # end try

      if(is(cytry, "try-error") == TRUE  ){
        next
      }else{
        fit.boots[[iboot]] <- resss.iboot$fit
        tau.boots[iboot,]  <- resss.iboot$tau
        if(boots$trace==TRUE){if(iboot %in% round(rev(seq(boots$nboot,1,by=-boots$nboot/5)))){cat(" - Already",round(iboot/boots$nboot*10)*10,"percent \n")}}
        iboot <- iboot + 1
      }

    }
    # for(iboot in 1:boots$nboot){
    #   # bootstrap idx of samples
    #   idx.iboot <- sample(1:N, replace = TRUE)
    #   # idx.iboot <- c(sample((1:N)[delta==1], replace = TRUE),sample((1:N)[delta==0], replace = TRUE))
    #   # solve the current estimator
    #   resss.iboot <- PH.Renew.Hetero.ALasso.Batch.fit(
    #     yobs=yobs[idx.iboot],delta=delta[idx.iboot],X=X[idx.iboot,,drop=F],initial=initial,prev=prev$fit.boots[[iboot]],
    #     tunpara.num=NULL,tunpara.min.ratio=NULL,tunpara=resss$tunpara.opt,
    #     evatype=evatype,learnrate=learnrate,maxit=maxit,eps=eps,threshold=threshold,
    #     trace=FALSE)
    #   fit.boots[[iboot]] <- resss.iboot$fit
    #   tau.boots[iboot,]  <- resss.iboot$tau
    #   if(boots$trace==TRUE){if(iboot %in% round(rev(seq(boots$nboot,1,by=-boots$nboot/5)))){cat(" - Already",round(iboot/boots$nboot*10)*10,"percent \n")}}
    # }
    bet.boots <- do.call(rbind,lapply(1:boots$nboot,function(iboot){fit.boots[[iboot]]$bet}))
    V.bet <- cov(bet.boots)*(boots$nboot-1)/boots$nboot
    bet.se <- sqrt(diag(V.bet))
    V.tau <- cov(tau.boots)*(boots$nboot-1)/boots$nboot
    tau.se <- sqrt(diag(V.tau))
  }

  # summary information for this current step
  # for bet
  zvalue.bet <- fit$bet/bet.se
  pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
  resbet <- data.frame(Est=fit$bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                       row.names=colnames(X))
  # for tau
  zvalue.tau <- resss$tau/tau.se
  pvalue.tau <- 2*(1-pnorm(abs(zvalue.tau)))
  restau <- data.frame(Est=resss$tau,SE=tau.se,zvalue=zvalue.tau,pvalue=pvalue.tau,
                       row.names=colnames(X))

  ### output: coef matrix + next step's summary information ###
  out <- list(
    resbet = resbet,
    restau = restau,
    bettau.path = resss$bettau.path,
    fit=fit,
    fit.boots=fit.boots
  )
  return(out)

}


#==== Online estimation in Cox proportional hazards model with dynamic covariate effects (one batch, Renewable, ALasso, no boot) ====#
#' @title Online estimation in Cox proportional hazards model with dynamic covariate effects (without bootstrap rocedure)
#'
#' @description Conduct Online estimation in Cox proportional hazards model with dynamic covariate effects (without bootstrap rocedure).
#'
#' @aliases PH.Renew.Hetero.ALasso.Batch.fit
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param initial a logical value. The default specification is \code{TRUE}, indicating that the current data batch is the first batch and no historical summary statistics is available.
#' @param prev indicates the historical summary data. The default specification is \code{NULL} (with \code{initial=TRUE}). Otherwise, it is a list with the following five elements:
#'   \code{bet} is the historical estimated result of regression coefficients;
#'   \code{bet.batch} is the historical estimated result of regression coefficients using only the previous dataset;
#'   \code{InfoM} is the historical estimated result of the information matrix;
#'   \code{N} is the sample of historical raw data.
#' @param tunpara.num the number of tuning parameter values that will be fitted at.
#' @param tunpara.min.ratio Smallest value for tuning parameter, as a fraction of tunpara.max, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero).
#' @param tunpara A user supplied \code{tunpara} sequence. Typical usage is to have the program compute its own \code{tunpara} sequence based on \code{tunpara.num} and \code{tunpara.min.ratio}. Supplying a value of \code{tunpara} overrides this.
#' @param evatype the information criterion that will be used in the selection of optimal tuning parameter.
#' @param learnrate the learning rate in solving the optimization problem.
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#' @param threshold a threshold value. The deafult value is \code{1e-5}. If the absolute value of the each component of coefficients is smaller than \code{threshold}, it will be converted into zero.
#' @param trace a logical value. The default specification is \code{TRUE}, indicating that several useful information about the current process of fitting will be displayed.
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{resbet}.
#'
#' @export PH.Renew.Hetero.ALasso.Batch.fit
PH.Renew.Hetero.ALasso.Batch.fit <- function(
    yobs,delta,X,initial=TRUE,prev=NULL,
    tunpara.num=50,
    tunpara.min.ratio=0.0001,
    tunpara=NULL,
    evatype=c("AIC","BIC")[2],
    learnrate=c(0.4),
    maxit=1e3,eps=1e-6,
    threshold=5e-5,
    trace=TRUE
){

  N <- length(yobs)

  ### calculate the updated estimator
  if(initial==TRUE){

    ### fit the model using classical method
    fit.init <- PH.fit(yobs,delta,X,maxit=maxit,eps=eps)
    # for bet
    bet <- fit.init$res[,1]
    bet.batch <- bet
    InfoM <- fit.init$InfoM
    # for tau
    tau <- rep(0,length(bet))
    tunpara.opt <- bettau.path <- NULL

  }else{

    ### do repeation (for better InfoM)
    pbet <- ncol(X)
    N <- length(yobs)
    bet.prev   <- prev$bet
    bet.batch.prev <- prev$bet.batch
    InfoM.prev <- prev$InfoM
    N.prev     <- prev$N

    ### other elements
    fitb <- PH.fit(yobs,delta,X,maxit=maxit,eps=eps)
    resbet.batch <- fitb$res
    bet.batch <- resbet.batch[,1]
    tau.nopen <- bet.batch - bet.batch.prev # bet.batch - bet.batch.prev #
    weights  <- abs(1/tau.nopen)

    ## prepare candidate set of tuning parameters
    if(is.null(tunpara)){
      Score.cc.tauu <- (N.prev/(N+N.prev))*InfoM.prev%*%(bet.batch-bet.prev-tau.nopen)
      InfoM.cc.tauu <- (N.prev/(N.prev+N))*InfoM.prev
      tunpara.max <- max(abs((diag(InfoM.cc.tauu)*tau.nopen+Score.cc.tauu)*tau.nopen))
      tunpara <- exp(seq(log(tunpara.max*2), log(tunpara.max * tunpara.min.ratio),
                         length.out = tunpara.num))
      bet.init <- rep(0,pbet)
      tau.init <- rep(0,pbet)
    }else{
      tunpara.num <- length(tunpara)
      bet.init <- bet.batch
      tau.init <- tau.nopen
    }

    ## prepare boxes for puting results
    bettau.path <- array(0,dim=c(2+2*pbet,tunpara.num),dimnames=list(
      c("IsConverge","EvaC",paste("bet",1:pbet,sep=""),paste("tau",1:pbet,sep="")),
      paste("Tuning",1:tunpara.num,sep="")))

    ## iteratively solving ALasso for each tuning parameter
    if(trace==TRUE){cat("Fit the penalized estimator: \n")}
    for(itunpara in 1:tunpara.num){ # itunpara <- 4
      if(trace==TRUE){cat(" - Current tuning parameter",itunpara,"/",length(tunpara),"(",tunpara[itunpara],") ... ")}

      tunpara.weights <- tunpara[itunpara]*weights

      # fit the model based on the current tuning parameter
      numit <- 1; bet.old <- bet.init; tau.old <- tau.init
      repeat{

        ScoreInfoM <- PH.ScoreInfoM(bet=bet.old,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
        InfoM.prev.bettau <- InfoM.prev%*%(bet.prev-bet.old+tau.old)
        Score.cc <- c(
          (N.prev/(N+N.prev))*InfoM.prev.bettau+(N/(N+N.prev))*ScoreInfoM$Score,
          -(N.prev/(N+N.prev))*InfoM.prev.bettau
        )
        InfoM.cc <- array(0,dim=c(2*pbet,2*pbet))
        InfoM.cc[1:pbet,1:pbet] <- (N/(N.prev+N))*ScoreInfoM$InfoM + (N.prev/(N.prev+N))*InfoM.prev
        InfoM.cc[(pbet+1):(2*pbet),1:pbet] <- InfoM.cc[1:pbet,(pbet+1):(2*pbet)] <- -(N.prev/(N.prev+N))*InfoM.prev
        InfoM.cc[(pbet+1):(2*pbet),(pbet+1):(2*pbet)] <- (N.prev/(N.prev+N))*InfoM.prev
        XXXtYYY <- InfoM.cc%*%c(bet.old,tau.old)+Score.cc

        numit.inner <- 1
        bet.inner.old <- bet.old; tau.inner.old <- tau.old
        repeat{
          # use Proximal Gradient Descent Here
          dev <- as.vector(learnrate * ( InfoM.cc%*%c(bet.inner.old,tau.inner.old) - XXXtYYY ))
          bet.inner <- bet.inner.old - dev[1:pbet]
          tau.inner.temp <- tau.inner.old - dev[-(1:pbet)]
          tau.inner <- Soft.Threshold(tau.inner.temp,learnrate*tunpara.weights)

          # judge the convergence (inner)
          if( max(abs(c(bet.inner,tau.inner)-c(bet.inner.old,tau.inner.old)))>eps & numit.inner<maxit ){
            bet.inner.old <- bet.inner
            tau.inner.old <- tau.inner
            numit.inner <- numit.inner + 1
          }else{
            bet.c <- bet.inner
            tau.c <- tau.inner
            break
          }
        }

        # judge the convergence (inner)
        if( max(abs(c(bet.c,tau.c)-c(bet.old,tau.old)))>eps & numit<maxit ){
          bet.old <- bet.c
          tau.old <- tau.c
          numit <- numit + 1
        }else{
          break
        }
      }

      # calculate evaluating criterion
      EvaC.c <- as.numeric(
        -PH.Renew.Hetero.betTau.Loglik(bet.c,tau.c,yobs=yobs,delta=delta,X=X,bet.prev,InfoM.prev,N.prev)+
          0.5*ifelse(evatype=="AIC",2,log(N+N.prev))*sum(abs(tau.c)>=threshold) )/(N+N.prev)
      if(trace==TRUE){ cat(EvaC.c,"\n") }
      # cat(EvaC.c,"\n")

      # combine and update the estimates
      bettau.path[,itunpara] <- c(numit<maxit,EvaC.c,bet.c,tau.c)
      bet.init <- bet.c
      tau.init <- tau.c
    }

    ## choose the best result according the chosen EvaC
    EvaC.best.idx <- which.min(bettau.path['EvaC',])
    bettau <- bettau.path[-c(1,2),EvaC.best.idx]; tunpara.opt <- tunpara[EvaC.best.idx]
    EvaC <- bettau.path[2,EvaC.best.idx]; convergence <- bettau.path[1,EvaC.best.idx]==1
    bet <- bettau[1:pbet]; tau <- bettau[-c(1:pbet)]; tau[abs(tau)<threshold] <- 0
    round(bet,5); round(tau,5)

    # current elements
    InfoMScore.curr <- PH.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,IsScore=TRUE,IsInfoM=TRUE)
    Score.curr      <- InfoMScore.curr$Score
    InfoM.curr      <- InfoMScore.curr$InfoM
    # update various elements
    InfoM      <- (N.prev/(N.prev+N))*InfoM.prev + (N/(N.prev+N))*InfoM.curr

  }

  ### output: coef matrix + next step's summary information ###
  out <- list(
    fit=list(
      bet = bet,
      bet.batch = bet.batch,
      InfoM  = InfoM,
      N = ifelse(initial==TRUE,N,N+N.prev)
    ),
    tunpara.opt = tunpara.opt,
    tau = tau,
    bettau.path = bettau.path
  )
  return(out)

}


#==== Loss function in Cox proportional hazards model (Online, Heterogeneity) ====#
#' @title Loss function in Cox proportional hazards model (Online, Heterogeneity)
#'
#' @description Calculate the loss function in Cox proportional hazards model (Online, Heterogeneity).
#'
#' @aliases PH.Renew.Hetero.betTau.Loglik
#'
#' @param bet unknown parameters corresponding to the model.
#' @param tau unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates
#' @param bet.prev the historical summary data of regression coefficients.
#' @param InfoM.prev the historical summary data of the information matrix (inverse variance-covariance matrix).
#' @param N.prev the sample of historical raw data
#'
#' @export PH.Renew.Hetero.betTau.Loglik
PH.Renew.Hetero.betTau.Loglik <- function(
    bet,tau,yobs,delta,X,bet.prev,InfoM.prev,N.prev
){

  # calculate the loss function for beta (log form)
  val.log <- PH.Beta.LogLik(bet=bet,yobs=yobs,delta=delta,X=X)
  val.log.online <-
    val.log -
    as.vector(t(bet-bet.prev-tau)%*%InfoM.prev%*%(bet-bet.prev-tau)) * N.prev / 2

  # output
  return(val.log.online)

}




#==== The soft-thresholding function ====#
#' @title The soft-thresholding function
#'
#' @description The soft-thresholding function.
#'
#' @aliases Soft.Threshold
#'
#' @param theta the parameter that we will calculate at.
#' @param lam the inputted tuning parameter
#'
#' @export Soft.Threshold
Soft.Threshold <- function(theta,lam){
  # the solution to LASSO penalty under the simplest situation: 2^{-1}(z-theta)^2+penalty
  # this is the equation (2.6) in Fan and Li (2001)
  res <- sign(theta)*pmax(abs(theta)-lam,0)
  return(res)
}






