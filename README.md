# TransHS
This is an R package for **Trans**ferring **H**eterogeneous **S**ources in survival targets based on Cox proportional hazards models with dymamic parameters.
- We note here the original name of this package was *StreamPHDE*  (see also the github website of Stat-WangXG) as was originally mentioned in the paper titled "*Renewable risk assessment of heterogeneous streaming time-to-event cohorts*", which has been accepted by *Statistics in Medicine* with DOI: 10.1002/sim.10146. In that paper, we focus on fitting Cox proportional hazards model with dynamic covariate effects in the heterogeneous streaming data environment.
- Now its name is *TransHS*, which emphasize that it can play a role in the topic about transfer leaning as well.

## Two main R function are included (with different scenarios for application):
- Analyze streaming survival data using online Cox proportional hazards model (assume homogeneity directly)
  - using *PH.Renew.Batch()*
- Analyze streaming survival data using online Cox proportional hazards model (with heterogeneity exists)
  - using *PH.Renew.Hetero.ALasso.Batch()*
  
## The corresponding example using simulated dataset is shown below:

- *An example for fitting the online Cox proportional hazards model (Heterogeneity exists)*
```R
## ---- basic setup of the simulated dataset ---- ##
pbet <- 3
B.oT <- 5 
Nb.oT <- 400
bettau.jump <- matrix(
  c(c( 1,-2, 0),  # bet0
    c(-1, 0, 2),  # bet1
    c( 1, 0, 0)), # bet2
  ncol=pbet)
Tnum <- nrow(bettau.jump)
bettau.all  <- matrix(0,nrow=B.oT*Tnum,ncol=pbet)
bettau.all[(0:(Tnum-1))*B.oT+1,]  <- bettau.jump
bet.all <- apply(bettau.all,2,cumsum) 

## ---- fit the  Cox proportional hazards model in an online manner ---- ##
# prepare basic elements
B  <- nrow(bet.all)
Nb <- Nb.oT
N <- B*Nb
Res.Online <- array(0,dim=c(B,pbet,2),
              dimnames=list(paste("batch",1:B,sep=""),paste("bet",1:pbet,sep=""),c("EST","SE")))
Res.Online_Hetero <- list(
  res = array(0,dim=c(B,pbet,2),
              dimnames=list(paste("batch",1:B,sep=""),paste("bet",1:pbet,sep=""),c("EST","SE"))),
  restau = array(0,dim=c(B-1,pbet,2),
                 dimnames=list(paste("batch",2:B,sep=""),paste("taub",1:pbet,sep=""),c("EST","SE")))
)
# Online estimation procedures
for(b in 1:B){ # b <- 1
  cat("-------------------------\nCurrent batch number: ",b,"\n")
  
  ## Generate the current data batch
  sdata.batch  <- sdata.PH(N=Nb,bet=bet.all[b,])
  yobs  <- sdata.batch$sdata$yobs
  delta <- sdata.batch$sdata$delta
  X     <- as.matrix(sdata.batch$sdata$X)
  
  ### ____ Online PH (homogeneous) ____ ###
  ## preparation: the previous elements
  if(b==1){ prevb <- NULL }else{
    prevb <- c(
      fitb$fit
    )
  }
  ## updating
  fitb <- PH.Renew.Batch(
    yobs=yobs,delta=delta,X=X,
    initial=(b==1),prev=prevb)
  Res.Online[b,,1] <- fitb$res[,1]
  Res.Online[b,,2] <- fitb$res[,2]
  
  ### ____ Online PH (heterogeneous) ____ ###
  ## preparation: the previous elements
  if(b==1){ prevb.Hetero <- NULL }else{
    prevb.Hetero <- c(
      fitb.Hetero$fit,
      fit.boots=list(fitb.Hetero$fit.boots)
    )
  }
  ## updating
  fitb.Hetero <- PH.Renew.Hetero.ALasso.Batch(
    yobs=yobs,delta=delta,X=X,
    initial=(b==1),prev=prevb.Hetero,
    boots=list(do=TRUE,nboot=100,trace=TRUE))
  Res.Online_Hetero$res[b,,1] <- fitb.Hetero$resbet[,1] # Res.Online_Hetero$res[,,1]
  Res.Online_Hetero$res[b,,2] <- fitb.Hetero$resbet[,2]
  if(b>1){
    Res.Online_Hetero$restau[b-1,,1] <- fitb.Hetero$restau[,1] # Res.Online_Hetero$restau[,,1]
    Res.Online_Hetero$restau[b-1,,2] <- fitb.Hetero$restau[,2]
  }

}
# present the fitted results
print(Res.Online)
print(Res.Online_Hetero$res)
print(Res.Online_Hetero$restau)
```


