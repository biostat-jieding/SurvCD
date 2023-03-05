#------------------------------------------------------------------------#
# data generating function ####
#------------------------------------------------------------------------#
#' @title The data generating function under the promotion time cure model (high-dimensional)
#'
#' @description Generate a simulated dataset from the promotion time cure model with relatively high dimensional covariates.
#'
#' @aliases sdata.PH
#'
#' @param N the sample size.
#' @param bet true value of coefficients. Certain default value is set.
#' @param rho the degree of correlation. Certain default value is set.
#' @param cvalue a value used to control the censoring rate. Certain default value is set.
#'
#' @export sdata.PH
sdata.PH <-  function(N,bet=c(1,-1,1),rho=0.5,cvalue=5){

  ### generate covariates
  pbet <- length(bet)
  Sigma <- rho^(abs(row(diag(pbet-1))-col(diag(pbet-1)))) # AR(1)
  X <- array(NA, dim=c(N,3))
  X[,1] <- rbinom(N,1,0.5)
  X[,2:pbet] <- mvtnorm::rmvnorm(N,mean=rep(0,nrow(Sigma)),sigma=Sigma)

  # set survival function (1-S(t))
  Ftx <- function(t,x,bet){ # population distribution function
    Lam0t <- t
    FF <- 1 - exp(-Lam0t*exp(sum(x*bet)) )
    return(FF)
  }

  # general process for generating survival time
  get.stime <- function(t,u,x,bet){ u - Ftx(t,x,bet) }
  stime <- rep(NA, N)
  for(i in 1:N){
    stime[i] <- uniroot(get.stime,c(0,1e6),u=runif(1,0,1),x=X[i,],bet=bet)$root
  }

  ### generate censoring time
  ctime <- rexp(N, 1/cvalue)
  delta <- as.numeric(stime<=ctime)
  censoring.rate <- 1 - mean(delta)

  # observed failure time
  yobs <- pmin(stime,ctime)

  # output
  out <- list(
    sdata=list(yobs=yobs,delta=delta,X=data.frame(X)),
    cenRt=censoring.rate
  ); out$cenRt


  return(out)

}
