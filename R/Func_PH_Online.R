

#========================================================#
# Fit the Cox proportional hazards model with streaming dataset ####
# -- a Renewable approach
#========================================================#

#==== Function for fitting online Cox proportional hazards model (one batch, Renewable approach) ====#
#' @title Online estimation in Cox proportional hazards model
#'
#' @description Fit online Cox proportional hazards model.
#' In other words, this function can analyze streaming survival data.
#'
#' @aliases PH.Renew.Batch
#'
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param initial a logical value. The default specification is \code{TRUE}, indicating that the current data batch is the first batch and no historical summary statistics is available.
#' @param prev indicates the historical summary data. The default specification is \code{NULL} (with \code{initial=TRUE}). Otherwise, it is a list with the following five elements:
#'   \code{bet} is the historical estimated result of regression coefficients;
#'   \code{InfoM} is the historical estimated result of the information matrix;
#'   \code{N} is the sample of historical raw data;
#' @param maxit specifies the maximum iteration number. If the convergence criterion is not met, the iteration will be stopped after emmax iterations and the estimates will be based on the last maximum likelihood iteration. The default \code{maxit = 5e3}.
#' @param eps tolerance for convergence. The default is \code{eps = 1e-6}. Iteration stops once the relative change in deviance is less than \code{eps}.
#'
#' @return The fitted results are returned (a list). The estimates of regression coefficients is \code{res}.
#'
#' @examples
#' ## ---- generate the dataset (full) ---- ##
#' N <- 10000
#' sdata.full  <- sdata.PH(N=N)
#' yobs  <- sdata.full$sdata$yobs
#' delta <- sdata.full$sdata$delta
#' X     <- as.matrix(sdata.full$sdata$X)
#' rm(sdata.full)
#'
#' ## ---- fit the Cox proportional hazards model in an online manner ---- ##
#' # prepare basic elements
#' B <- 40
#' Nb <- N/B
#' batch <- ceiling((1:N) / Nb)
#' Res.Online <- array(0,dim=c(B,3,2),
#'                     dimnames=list(paste("batch",1:B,sep=""),paste("bet",1:3,sep=""),c("EST","SE")))
#' # the online procedure
#' for(b in 1:B){  # b <- 1
#'   cat("Batch",b,"...\n")
#'
#'   # preparation: the batch idx / the previous elements
#'   idxb <- batch==b
#'   if(b==1){ prevb <- NULL }else{
#'     prevb <- c(
#'       fitb$fit
#'     )
#'   }
#'
#'   # fit the current data batch (with current data the historical statistics)
#'   fitb <- PH.Renew.Batch(
#'     yobs=yobs[idxb],delta=delta[idxb],X=X[idxb,,drop=F],
#'     initial=(b==1),prev=prevb)
#'   Res.Online[b,,1] <- fitb$res[,1]
#'   Res.Online[b,,2] <- fitb$res[,2]
#'
#' }
#' # present the fitted results
#' print(Res.Online)
#'
#' @export PH.Renew.Batch
PH.Renew.Batch <- function(yobs,delta,X,initial=TRUE,prev=NULL,maxit=1e3,eps=1e-6){

  ### calculate the updated estimator
  if(initial==TRUE){

    ### fit the model using classical method
    fit.init <- PH.fit(yobs,delta,X,eps=eps,maxit=maxit) # internal dataset's full estiamtor
    # for bet
    bet <- fit.init$res[,1]
    bet.se <- fit.init$res[,2]
    InfoM <- fit.init$InfoM
    Score <- fit.init$Score
    N.prev <- 0
    N    <- length(yobs)


  }else{

    ### do repeation (for better InfoM)
    pbet <- ncol(X)
    N    <- length(yobs)
    bet.prev       <- prev$bet
    InfoM.prev     <- prev$InfoM
    N.prev         <- prev$N

    # updated estimator based on newton method: similar to Luo and Song (2020)
    solve.InfoM <- solve(
      N*PH.ScoreInfoM(bet=bet.prev,yobs=yobs,delta=delta,X=X,IsScore=FALSE,IsInfoM=TRUE)$InfoM +
        N.prev*InfoM.prev)
    numit <- 1; bet.old <- bet.prev
    repeat{
      Score <- PH.ScoreInfoM(bet=bet.old,yob=yobs,delta=delta,X=X,IsScore=TRUE)$Score
      dev <- solve.InfoM%*%(N.prev*InfoM.prev%*%(bet.prev-bet.old)+N*Score)
      bet <- bet.old + as.vector(dev)
      if( max(abs(bet-bet.old))>eps & numit<maxit ){
        bet.old <- bet
        numit <- numit + 1
      }else{
        break
      }
    }

    ### calculate SEs for bet using explicit formula !!!
    # current elements
    InfoM.curr      <- PH.ScoreInfoM(bet=bet,yobs=yobs,delta=delta,X=X,IsScore=FALSE,IsInfoM=TRUE)$InfoM
    # update various elements
    InfoM      <- (N.prev/(N.prev+N))*InfoM.prev      + (N/(N.prev+N))*InfoM.curr
    # calculate SEs for bet
    bet.se <- sqrt(diag(solve(InfoM))/(N.prev+N))

  }

  # summary information for this current step
  # for bet
  zvalue.bet <- bet/bet.se
  pvalue.bet <- 2*(1-pnorm(abs(zvalue.bet)))
  res <- data.frame(Est=bet,SE=bet.se,zvalue=zvalue.bet,pvalue=pvalue.bet,
                    row.names=colnames(X))

  ### output: coef matrix + next step's summary information ###
  out <- list(
    res = res,
    fit=list(
      bet=bet,
      InfoM = InfoM,
      N = ifelse(initial==TRUE,N,N+N.prev)
    )
  )
  return(out)

}



#==== Loss function for bet (Profiled online log-likelihood function) ====#
#' @title Loss function for bet in Cox proportional hazards model (Online version)
#'
#' @description  Loss function for the Cox proportional hazards model (Profiled online log-likelihood function in NPMLE).
#' Note that it highly depends on the function \code{\link{PH.Beta.LogLik}} in offline scenario.
#'
#' @aliases PH.Renew.Beta.LogLik
#'
#' @param bet unknown parameters corresponding to the model.
#' @param yobs time to event of interest.
#' @param delta the censoring indicator, normally 1 = event of interest happens, and 0 = censoring.
#' @param X a matrix of covariates.
#' @param bet.prev the historical summary data of regression coefficients.
#' @param ScoreAppro.prev the historical summary data of a term derived from the online variable selection procedure (the recovering term).
#' @param InfoM.prev the historical summary data of the information matrix (inverse variance-covariance matrix).
#' @param N.prev the sample of historical raw data
#'
#' @export PH.Renew.Beta.LogLik
PH.Renew.Beta.LogLik <- function(bet,yobs,delta,X,bet.prev,InfoM.prev,ScoreAppro.prev,N.prev){

  # calculate the loss function for beta (log form)
  val.log <- PH.Beta.LogLik(bet=bet,yobs=yobs,delta=delta,X=X)
  val.log.online <- val.log -
    as.vector(t(bet-bet.prev)%*%(InfoM.prev%*%(bet-bet.prev)-2*ScoreAppro.prev)) * N.prev / 2

  # output
  return(val.log.online)

}


