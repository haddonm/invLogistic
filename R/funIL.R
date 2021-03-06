
## return the IL coefficients from an IL object


#' @title coef.IL return the estimated IL parameters from an IL object
#'
#' @param object an IL object from the fitIL function
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return a vector of the estimated IL parameters from an IL object
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=TRUE,sitename="Middle Ground")
#' coef(ans)
coef.IL <- function(object,...) {
  ans <- object$model$estimate
  names(ans) <- c("MaxDL","L50","L95","MaxSig")
  return(ans)
}

#' @title dobootIL conducts a bootstrap analysis on an IL object
#'
#' @description dobootIL conducts a bootstrap analysis on an IL object
#'     so as to characterize the uncertainty in the analysis. In the
#'     process it produces a 'bootIL' object. This contains 'reps'
#'     replicate sets of parameter estimates which can be plotted
#'     using the S3 plot method for bootIL objects.
#'
#' @param x an IL object such as produced by fitIL
#' @param reps how many boostraps to produce, 1000 should be regarded
#'     as a minimum but that will usually take a litle time to run
#'
#' @return a class bootIL object containing the bootstrap replicate
#'     parameter estimates and the median plus lower and upp 95% bounds
#'     of each parameter.
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=TRUE,sitename="Middle Ground")
#' boot <- dobootIL(ans,reps=10)
#' str(boot)
dobootIL <- function(x,reps=1000) {
  columns <- c("MaxDL","MDLL95%","MDL50%","MDLU95%","L50","L50L95%",
               "L5050%","L50U95%","L95","L95L95%","L9550%","L95U95%",
               "MaxSig","MSL95%","MS50%","MSU95%","Obs","Errors")
  bootans <- matrix(0,nrow=1,ncol=18,dimnames=list(x$siteid,columns))
  mainpar <- c(1,5,9,13)
  outboot <- matrix(0,nrow=reps,ncol=4,dimnames=list(seq(1,reps,1),
                                   c("MaxDL","L50","L95","MaxSig")))
  bootans[mainpar] <- x$model$estimate
  xLt <- x$Lt
  yDL <- x$DL
  nb <- x$Nobs
  bootans[17] <- nb
  pick <- seq(1,nb,1)
  bootans[18] <- 0
  doboot <- function(inpick,inLt,inDL) { # internal function
    boots <- sample(pick,replace=T)
    boots <- boots[order(boots)]
    model <- fitIL(xLt[boots],yDL[boots])
    return(model$model$estimate)
  }
  res <- lapply(1:reps, function(i) try(doboot(pick,xLt,yDL), TRUE))
  for (bootnum in 1:reps) {
    if (is.numeric(res[[bootnum]])) {
      outboot[bootnum,] <- res[[bootnum]]
    } else {
      bootans[18] <- bootans[18] + 1
      outboot[bootnum,] <- c(NA,NA,NA,NA)
    }
  }
  initcol <- -2
  for (index in 1:4) {
    CIs <- quantile(outboot[,index],probs=c(0.025,0.5,0.975),na.rm=T)
    initcol <- initcol + 4
    bootans[initcol:(initcol+2)] <- CIs
  }
  ans <- list(bootans,outboot)
  names(ans) <- c("Percentiles","Replicates")
  class(ans) <- "bootIL"
  return(ans)
} # end of dobootIL

#' @title fitIL S3 generic method for fitting an IL curve to data
#'
#' @description fitIL is an overloaded front-end function used to fit
#'     an IL curve to different data types containing at least two
#'     vectors of data required: x=Lt (the initial lengths) and y=DL
#'     (the corresponding growth increments). These data types include
#'     a matrix of x,y, a list of x, y, and the default, which is two
#'     vectors x=Lt, and y=DL. This is an S3 generic method, whose
#'     only task is to direct the analysis to the correct method by
#'     class of data input or go  off to the default.
#'
#' @param x a data type containing, first, the Initial lengths, Lt,
#'     and, second, the growth increments, DL. These can be passed
#'     as separate vectors, as a matrix, a list, or a data.frame.
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return returns a fitted IL curve of class IL
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=FALSE,)
#' str(ans)
fitIL <- function(x, ...) {
  if(is.null(class(x))) class(x) <- data.class(x)
  UseMethod("fitIL", x)
}

## Alternative fitIL Methods

#' @title fitIL.matrix S3 method to apply fitIL to a matrix
#'
#' @param x a matrix containing Lt and DL in the first two columns
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return if x is a matrix this fits an IL curve
#' @export
#'
#' @examples
#' data(midg)
#' fitIL(cbind(midg$Lt,midg$DL))
fitIL.matrix <- function(x, ...) {
  fitIL(x[,1],x[,2], ...)
}

#' @title fitIL.list S3 method to apply fitIL to a matrix
#'
#' @param x a list containing Lt and DL in the first two components
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return if x is a list this fits an IL curve
#' @export
#'
#' @examples
#' data(midg)
#' fitIL(list(midg$Lt,midg$DL))
fitIL.list <- function(x, ...) {
  fitIL(x[[1]], x[[2]], ...)
}

#' @title fitIL.data.frame S3 method to apply fitIL to a matrix
#'
#' @param x a data.frame containing Lt and DL. The names must be exact
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return if x is a data.frame this fits an IL curve
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,siteid=1,outliers=TRUE,sitename="midground")
#' str(ans)
fitIL.data.frame <- function(x, ...) {
  fitIL(x=x[,"Lt"], y=x[,"DL"], ...)
}

## Fit Inverse Logistic to x=Lt and y=DL data, make an object of class IL
##   x=midg$Lt; y = midg$DL; siteid=0; outliers=FALSE; sitename="Middle Ground"

#' @title fitIL.default the default S3 method to apply fitIL
#'
#' @param x a vector containing Lt the names must be exact
#' @param y a vector containing DL the names must be exact
#' @param siteid a numeric site identifier, default = 0
#' @param outliers default = FALSE, should outliers be identified and
#'     then omitted?
#' @param sitename an alphanumeric identifier for the site
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return returns an IL curve fit as a class IL object
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg$Lt,midg$DL,sitename="Middle Ground")
#' str(ans)
fitIL.default <- function(x, y, siteid=0,outliers=FALSE,sitename="",...) {
  negLIL <- function(parsin) {
    expDL <- invlog(p=parsin,x)
    expSD <- invlog(p=c(parsin[4],parsin[3],parsin[3]/0.95),x)
    neglogl <- -sum(dnorm(y,expDL,expSD,log=T))
    return(neglogl)
  }
  parsin <- initpars(x,y)
  best <- optim(parsin,negLIL,method="Nelder-Mead",
                hessian=FALSE,
                control=list(trace=0, maxit=1000))
  parsin <- best$par
  mod <- nlm(negLIL,parsin,hessian=T, gradtol = 1e-7)
  parsin <- mod$estimate
  MaxDL <- mod$estimate[1]
  L50 <- mod$estimate[2]
  L95 <- mod$estimate[3]
  MaxSig <- mod$estimate[4]
  xout <- NULL  # will contain the list of outliers if one exists
  yout <- NULL
  L50out <- NULL
  L95out <- NULL
  MaxDLout <- NULL
  MaxSigout <- NULL
  if (outliers) {
    L50out <- L50
    L95out <- L95
    MaxDLout <- MaxDL
    MaxSigout <- MaxSig
    expDL <-  invlog(c(MaxDL,L50,L95),x)
    resids <- abs(y - expDL)
    expSD <- invlog(c(MaxSig,L95,(L95/0.95)),x)
    outers <- resids - 2.576*expSD   #99% confidence limits
    pick <- which(outers > 0)
    if ((length(pick) > 0)==TRUE) {
      xout <- x[pick]
      yout <- y[pick]
      x <- x[-pick]
      y <- y[-pick]
    }
    best <- optim(parsin,negLIL,method="Nelder-Mead",
                  hessian=FALSE,
                  control=list(trace=0, maxit=2000))
    parsin <- best$par
    mod <- nlm(negLIL,parsin,hessian=T, gradtol = 1e-7)
    MaxDL <- mod$estimate[1]
    L50 <- mod$estimate[2]
    L95 <- mod$estimate[3]
    MaxSig <- mod$estimate[4]
  }
  Ltrg <- range(x,na.rm=T)
  xmin <- min(Ltrg[1],50)
  xmax <- max(Ltrg[2],180)
  predLt <- seq(xmin,xmax,1)
  predDL <- invlog(c(MaxDL,L50,L95),predLt)
  Nobs <- length(x)
  ans <- list(mod,MaxDL,L50,L95,MaxSig,predLt,predDL,Nobs,x,y,xout,
              yout,L50out,L95out,MaxDLout,MaxSigout,siteid,sitename)
  names(ans) <- c("model","MaxDL","L50","L95","MaxSig","PredLt",
                  "PredDL","Nobs","Lt","DL","OutLt","OutDL","L50out",
                  "L95out","MaxDLout","MaxSigout","siteid","sitename")
  class(ans) <- "IL"
  return(ans)
} # end of fitIL.default


#' @title getdyn estimates age-structure given tagging growth fit
#'
#' @description getdyn uses the tagging growth fit, the implied length
#'     at age 0+, and the maximum age, to estimate the length-at-age
#'     and the growth increment at average length at age. It enables
#'     an estimation of the age at maturity and the length of animals
#'     at age maxage-1
#'
#' @param x the IL object from fitIL
#' @param initL the length of post-larval anaimals, default = 2.0mm
#' @param maxage the maximum age+1 for the calculations
#'
#' @return a list containing the length at maximum age, the age at
#'     maturity, the length-at-age, and gradient between the
#'     interquartile growth increments
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=TRUE,sitename="Middle Ground")
#' getdyn(ans,maxage=30)
getdyn <- function(x,initL=2.0,maxage=24) {  # x=ans; initL=2.0; maxage=24
  LaA <- matrix(0,nrow=maxage,ncol=2,dimnames=list(seq(1,maxage,1)-1,
                                                   c("DLarAge","LatAge")))
  LaA[1,2] <- initL
  for (count in 1:(maxage-1)) { #count = 1
    Lt <- LaA[count,2]
    DL <- invlog(c(x$MaxDL,x$L50,x$L95),Lt)
    LaA[count,1] <- DL
    Lt1 <- Lt + DL
    if (Lt1 <= x$L50) {
      SaM <- count + 1 # to account for age 0
    }
    LaA[count+1,2] <- LaA[count,2] + LaA[count,1]
  }
  down <- LaA[SaM,2]
  up <- LaA[SaM+1,2]
  dif <- up - down
  prop <- x$L50 - down
  AgeM <- SaM + prop/dif
  Len23 <- LaA[maxage,2]
  LIQ <- 0.75*x$MaxDL  # InterQuartile Difference
  UIQ <- 0.25*x$MaxDL
  pick <- which(x$PredDL < LIQ)
  low <- x$PredLt[pick[1]-1]
  pick <- which(x$PredDL < UIQ)
  high <- x$PredLt[pick[1]-1]
  iqg <- (LIQ-UIQ)/(high - low)
  out <- list(Len23,AgeM,LaA,iqg)
  names(out) <- c("LenMaxA","AgeM","LaA","InterQG")
  return(out)
} # end of getdyn


#' @title ILfuns prints the syntax of the most used functions in invLogistic
#'
#' @description ILfuns prints the syntax of the most used functions in
#'     invLogistic. This is designed to assist with using the syntax.
#'
#' @return nothing, but it does print some lines to the console
#' @export
#'
#' @examples
#' ILfuns()
ILfuns <- function() {
  print("Lt = vector of initial lengths, DL = vector of growth increments",quote=FALSE)
  print("invlog(p,Lt)",quote=FALSE)
  print("fitIL(x, siteid=0,outliers=FALSE,sitename='')",quote=FALSE)
  print("or", quote=FALSE)
  print("fitIL(Lt, DL, siteid=0,outliers=FALSE,sitename='')",quote=FALSE)
  print("dobootIL(output from fitIL,reps=1000)",quote=FALSE)
  print("summary(output from fitIL)",quote=FALSE)
  print("plot(output from fitIL)",quote=FALSE)
  print("plot(output from dobootIL,col=0,font=7)",quote=FALSE)
  print("print - or just type the name of the output from fitIL",quote=FALSE)
  print("initpar(Lt,DL)",quote=FALSE)
}


#' @title initpars estimates the starting parameters for fitting an IL curve
#'
#' @description initpars estimates the starting parameters for fitting
#'     an IL curve. This uses approximations and rules of thumb.
#'
#' @param Ltin a vector of initial lengths to which the IL is to be fitted
#' @param DLin a vector of growth increments to which the IL is to be
#'     fitted
#'
#' @return a vector of four parameters for the IL (MaxDL, L50, L95, MaxSig)
#' @export
#'
#' @examples
#'   data(midg)
#'   initpars(midg$Lt,midg$DL)
initpars <- function(Ltin,DLin) {
  pars <- numeric(4)
  minL <- min(Ltin,na.rm=T)
  maxL <- max(Ltin,na.rm=T)
  extent <- maxL - minL
  lim <- minL+0.2*extent
  pick <- which(Ltin < lim)
  if (length(pick) > 1) {
    pars[1] <- mean(DLin[pick])
    pars[4] <- sd(DLin[pick])
  } else {
    pars[1] <- 22.0
    pars[4] <- 4.5
  }
  pars[2] <- minL + 0.5*extent
  pars[3] <- minL + 0.9*extent
  names(pars) <- c("MaxDL","L50","L95","SigMax")
  return(pars)
} # end of initpars

## The inverse logistic functions
#invlog <- function(x, ...) {
#  if(is.null(class(x))) class(x) <- data.class(x)
#  UseMethod("invlog",x)
#}

#' @title invlog calculates the mean growth increment for the IL
#'
#' @description invlog is the simplest equation dealing with the
#'     inverse logistic equation
#'
#' @param p a vector of parameters with at least c(MaxDL,L50,L95)
#' @param Lt the initial lengths, Lt, for which to calculate the
#'     predicted growth increment
#'
#' @return a vector of growth increments the same length as x
#' @export
#'
#' @examples
#'  data(midg)
#'  pars <- initpars(midg$Lt,midg$DL)
#'  invlog(pars,sort(midg$Lt[1:10]))
invlog <- function(p,Lt) {
  ans <- p[1]/(1+exp(log(19)*(Lt-p[2])/(p[3]-p[2])))
  return(ans)
}

#invlog.vector <- function(x, ...) invlog(x[1],x[2],x[3],...)

#' @title plot.bootIL an S3 method to plot the bootstrap results
#'
#' @description plot.bootIL an S3 method for plotting the results of
#'     bootsrapping a inverse logistic model fit to a set of tagging\
#'     data. A 2 x 2 plot of histograms is generated illustrating the
#'     distribution of the bootstrap replicate parameter estimates.
#'     This has no net effect on the par defnition.
#'
#' @param x the bootIL object from the dobootIL function
#' @param col the col of the histogram cells, default = 0 = empty
#' @param font default=7, bold serif, 1 is sans-serif
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return nothing but oit does generate a plot
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=TRUE,sitename="Middle Ground")
#' boot <- dobootIL(ans,reps=10)
#' plot(boot)
plot.bootIL <- function(x,col=0,font=7,...) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  collab = c("MaxDL","L50","L95","MaxSig")
  initcol <- -2
  sitenum <- rownames(x$Percentiles)
  # if (length(grep("font",ls())) == 0) font <- 7
  # if (length(grep("col",ls())) == 0) col <- 0
  par(mfrow = c(2,2))
  par(mai=c(0.5,0.5,0.1,0.1), oma=c(0,0,2,0),cex=0.8)
  par(mgp=c(1.35,0.35,0),font.axis=font,font=font,font.lab=font)
  for (index in 1:4) {
    initcol <- initcol + 4
    CIs <- x$Percentiles[initcol:(initcol+2)]
    hist(x$Replicates[,index],xlab="",main="",col=col)
    title(xlab=list(collab[index], cex=1.0,font=font))
    abline(v=CIs[1],col=2,lty=2)
    abline(v=CIs[2],col=4,lty=2)
    abline(v=CIs[3],col=2,lty=2)
  }
  mtext(paste("Site ",sitenum,sep=""),side=3,line=0.5,outer=T,font=7,
        cex=1.25)
} # end of S3 plot.bootIL


#' @title plot.IL an S3 method to plot an IL object
#'
#' @description plot.IL is an S3 method for plotting the results of
#'     fitting an inverse logistic curve to a set of tagging data.
#'     it plot the data and the model fit, with curves depicted the
#'     90th percentile spread and, if outliers was set to TRUE in the
#'     call to fitIL, it identifies outliers as red dots, and removes
#'     those points from the calculation. It plots the residuals, the
#'     rate of change in DL, to find the turnover length, and the
#'     density of data by initial length
#'
#' @param x an IL object from the fitIL function
#' @param outliers identify outliers on the plot. Default = FALSE
#' @param maxx allows one to set a particular maximum growth increment.
#'     default=0 so that ymax is set by the data
#' @param miny allows one to set the lower bound of the growth plots. Default=
#'     -3 as befits abalone
#' @param minx allows one to set the lower bound of the x-axis on growth plots.
#'     default=0
#' @param nbreaks number of bins in the histogram of data density, default=25
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return nothing but it does generate a plot
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=TRUE,sitename="Middle Ground")
#' plot(ans)
plot.IL <- function(x,outliers=FALSE,maxx=0,miny=-3,minx=0,nbreaks=25,...) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
  plotmodelIL(x,outliers=outliers,maxx=maxx,ymin=miny,xmin=minx)
  expDL <-  invlog(c(x$MaxDL,x$L50,x$L95),x$Lt)
  resids <- x$DL - expDL
  expSD <- invlog(c(x$MaxSig,x$L95,210),x$PredLt)
  outer99 <- 2.5760 * expSD
  outer90 <- 1.965 * expSD
  par(mfrow = c(2,2))
  par(mai=c(0.2,0.4,0.1,0.1), oma=c(2,0,2,0))
  par(cex=0.8, mgp=c(1.35,0.35,0),font.axis=7,font=7,font.lab=7)
  # Plot the basic fit with outliers if any
  ymax <- max(x$DL,x$OutDL)*1.025
  xmax <- max(max(x$Lt,x$OutLt)*1.025,maxx)
  xmin <- min(min(x$Lt,x$OutLt) - 1,minx)
  plot(x$Lt,x$DL,type="p",pch=20,xlab="",ylab="",xaxs="r",yaxs="r",
       xlim=c(xmin,xmax),ylim=c(miny,ymax))
  lines(x$PredLt,x$PredDL,col=2,lwd=2)
  lines(x$PredLt,x$PredDL+outer99,col=2,lty=2)
  lines(x$PredLt,x$PredDL-outer99,col=2,lty=2)
  lines(x$PredLt,x$PredDL+outer90,col=4,lty=2)
  lines(x$PredLt,x$PredDL-outer90,col=4,lty=2)

  abline(h=0,col="grey")
  abline(h=miny,col="grey")
  if (length(x$OutLt)>0) {
    points(x$OutLt,x$OutDL,col=2,pch=20)
  }
  title(ylab=list("Growth Increment DL", cex=1.1, col=1, font=7))
  text(170,0.95*ymax,round(x$MaxDL,3),cex=0.8,font=7)
  text(170,0.9*ymax,round(x$L50,3),cex=0.8,font=7)
  text(170,0.85*ymax,round(x$L95,3),cex=0.8,font=7)
  text(170,0.8*ymax,round(x$MaxSig,3),cex=0.8,font=7)
  # Plot the residuals
  plot(x$Lt,resids,type="p",pch=20,xlab="",ylab="",xaxs="r",yaxs="r",
       xlim=c(xmin,xmax))
  lines(x$PredLt,outer99,col=2,lty=2)
  lines(x$PredLt,-outer99,col=2,lty=2)
  lines(x$PredLt,outer90,col=4,lty=2)
  lines(x$PredLt,-outer90,col=4,lty=2)
  abline(h=0,col="grey")
  title(ylab=list("Residuals mm", cex=1.1, col=1, font=7))
  # Plot the rate of change in DL
  N <- length(x$PredDL)
  diffDL <- numeric(N-1)
  for (index in 2:N) {
    diffDL[index-1] <- x$PredDL[index] - x$PredDL[index-1]
  }
  plot(x$PredLt[1:(N-1)],diffDL,type="l",xlim=c(xmin,xmax),xlab="",
       ylab="")
  abline(v=x$L50,col=2)
  title(ylab=list("Rate of Change of DL", cex=1.1, col=1, font=7))

  bins <- seq(xmin,xmax,length=nbreaks)
  hist(x$Lt,breaks=bins,xlab="",ylab="",main="")
  abline(v=x$L50,col="grey")
  title(ylab=list("Density of Data Points", cex=1.1, col=1, font=7))
  mtext("Initial Length Lt",side=1,line=0.5,outer=T,font=7,cex=1.25)
  if (nchar(x$sitename) > 0) label=x$sitename
     else label=paste0("Site ",x$siteid)
  mtext(label,side=3,line=0.5,outer=T,font=7,cex=1.25)
} # end of plot.IL

#' @title plotmodelIL plots out the model, data, and fit
#'
#' @description plotmodelIL plots the data, the model fit, and the
#'     90th and 99th percentile bounds on the expected distribution
#'     of observations.
#'
#' @param x an IL object as produced by fitIL
#' @param outliers identify outliers on the plot. Default = FALSE
#' @param maxx The maximum value for the x axis. default=180
#' @param ymin allows one to set the lower bound of the growth plots. Default=
#'     -3 as befits abalone
#' @param xmin allows one to set the lower bound of the x-axis on growth plots.
#'     default=0
#' @param defpar should plot par values be used, default = FALSE
#'
#' @return nothing but does plot a graph
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=TRUE,sitename="Middle Ground")
#' plotmodelIL(ans,outliers=TRUE,defpar=TRUE)
plotmodelIL <- function(x,outliers=FALSE,maxx=0,ymin=-3,xmin=0,defpar=FALSE) {
  expDL <-  invlog(c(x$MaxDL,x$L50,x$L95),x$Lt)
  expSD <- invlog(c(x$MaxSig,x$L95,210),x$PredLt)
  outer99 <- 2.5760 * expSD
  outer90 <- 1.965 * expSD
  if (defpar) {
    par(mfrow = c(1,1))
    par(mai=c(0.4,0.4,0.1,0.1), oma=c(0,0,0,0))
    par(cex=0.8, mgp=c(1.35,0.35,0), font.axis=7)
  }
  # Plot the basic fit with outliers if any
  ymax <- getmax(c(x$DL,x$OutDL),mult=1.025)
  xmax <- getmax(c(x$Lt,x$OutLt,maxx),mult=1.025)
  xmin <- getmin(c(c(c(x$Lt,x$OutLt) - 1),xmin),mult=1.025)
  plot(x$Lt,x$DL,type="p",pch=20,xlab="",ylab="",xaxs="r",yaxs="r",
       xlim<- c(xmin,xmax),ylim=c(ymin,ymax),panel.first=grid())
  lines(x$PredLt,x$PredDL,col=2,lwd=2)
  lines(x$PredLt,x$PredDL+outer99,col=2,lty=2)
  lines(x$PredLt,x$PredDL-outer99,col=2,lty=2)
  lines(x$PredLt,x$PredDL+outer90,col=4,lty=2)
  lines(x$PredLt,x$PredDL-outer90,col=4,lty=2)
  abline(h=0,col="grey")
  abline(h=ymin,col="grey")
  if ((length(x$OutLt)>0) & (outliers == TRUE)) {
    points(x$OutLt,x$OutDL,col=2,pch=20)
  }
  title(ylab=list("Growth Increment DL", cex=1.0, col=1, font=7),
        xlab=list("Shell Length (mm)", cex=1.0, col=1, font=7))
  xloc <- 0.9*xmax
  text(xloc,0.95*ymax,round(x$MaxDL,3),cex=0.8,font=7,pos=4)
  text(xloc,0.875*ymax,round(x$L50,3),cex=0.8,font=7,pos=4)
  text(xloc,0.8*ymax,round(x$L95,3),cex=0.8,font=7,pos=4)
  text(xloc,0.725*ymax,round(x$MaxSig,3),cex=0.8,font=7,pos=4)
} # end of plotmodelIL

##' @title print.IL an S3 method ofr IL objects
#'
#' @description print.IL is an S3 method for IL objects, as generated
#'     by the fitIL function
#'
#' @param x an IL objects, as generated by fitIL
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return nothing but it prints a few details to the console
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,sitename="midgrd")
#' print(ans)
print.IL <- function(x, ...) {
  cat("\n Site Id  : ",x$siteid)
  cat("\n sitename : ",x$sitename)
  cat("\n MaxDL    : ",x$MaxDL)
  cat("\n L50      : ",x$L50)
  cat("\n L95      : ",x$L95)
  cat("\n MaxSig   : ",x$MaxSig)
  cat("\n ")
} # end of print.IL

#' @title summary.IL an S3 method for the summary generic for IL objects
#'
#' @description summary.IL is an S3 method for the summary generic for
#'     IL objects, as generated by the fitIL function.
#'
#' @param object an IL object generated by fitIL
#' @param console print the summary to console, default=TRUE
#'
#' @return a vector of character strings. and can prints summary to the console
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,sitename="middleground")
#' out <- summary(ans,console=TRUE)
summary.IL <- function(object,console=TRUE) {
  txt <- vector("character",16)
  Ltrge <- range(object$Lt,na.rm=TRUE)
  DLrge <- range(object$DL,na.rm=TRUE)
  outs <- FALSE
  if (length(object$MaxDLout) > 0) { outs <- TRUE }
  txt[1] <- paste0("siteid  : ",object$siteid)
  txt[2] <- paste0("sitename: ",object$sitename)
  if (outs){ txt[3] <- paste0("MaxDL   : ",round(object$MaxDL,digits=4),"   ",
                              round(object$MaxDLout,digits=4))
  } else { txt[3] <- paste0("MaxDL   : ",round(object$MaxDL,digits=4))
  }
  if (outs){ txt[4] <- paste0("L50     : ",round(object$L50,digits=4),"  ",
                              round(object$L50out,digits=4))
  } else { txt[4] <- paste0("L50     : ",round(object$L50,digits=4))
  }
  if (outs){ txt[5] <- paste0("L95     : ",round(object$L95,digits=4)," ",
                              round(object$L95out,digits=4))
  } else { txt[5] <- paste0("L95     : ",round(object$L95,digits=4))
  }
  if (outs){ txt[6] <- paste0("MaxSig  : ",round(object$MaxSig,digits=4),"   ",
                              round(object$MaxSigout,digits=4))
  } else { txt[6] <- paste0("MaxSig  : ",round(object$MaxSig,digits=4))
  }
  txt[7] <- paste0("N       : ",object$Nobs)
  txt[8] <- paste0("Outliers: ",length(object$OutLt))
  txt[9] <- paste0("Range Lt: ",Ltrge[1],"    ",Ltrge[2])
  txt[10] <- paste0("Range DL: ",DLrge[1],"    ",DLrge[2])
  txt[11] <- paste0("-ve LL  : ",round(object$model$minimum,5))
  txt[12] <- paste0("Other Components")
  txt[13] <- paste0("$Lt and $DL are the input data minus any outliers")
  txt[14] <- paste0("$model  contains the nlm fit")
  txt[15] <- paste0("$PredLt and PredDL = fitted line")
  txt[16] <- paste0("$OutLt and $OutDL = outlier values")
  if (console) for (i in 1:16) cat(txt[i],"\n")
  return(txt)
} # end of summary.IL
