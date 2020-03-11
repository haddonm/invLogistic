
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
#' ans <- fitIL(midg)
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
#' ans <- fitIL(midg)
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
    expSD <- invlog(p=c(parsin[4],parsin[3],210),x)
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
    expSD <- invlog(c(MaxSig,L95,210),x)
    outers <- resids - 2.576*expSD   #99% confidence limits
    pick <- which(outers > 0)
    if ((length(pick) >0)==TRUE) {
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
  print("invlog(MaxDL,L50,L95,Lt)",quote=FALSE)
  print("fitIL(x, siteid=0,outliers=FALSE,sitename='')",quote=FALSE)
  print("or", quote=FALSE)
  print("fitIL(Lt, DL, siteid=0,outliers=FALSE,sitename='')",quote=FALSE)
  print("dobootIL(output from fitIL)",quote=FALSE)
  print("summary(output from fitIL)",quote=FALSE)
  print("plot(output from fitIL)",quote=FALSE)
  print("plot(output from dobootIL)",quote=FALSE)
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
#' @param x the initial lengths, Lt, for which to calculate the
#'     predicted growth increment
#'
#' @return a vector of growth increments the same length as x
#' @export
#'
#' @examples
#'  data(midg)
#'  pars <- initpars(midg$Lt,midg$DL)
#'  invlog(pars,sort(midg$Lt[1:10]))
invlog <- function(p,x) {
  ans <- p[1]/(1+exp(log(19)*(x-p[2])/(p[3]-p[2])))
  return(ans)
}

#invlog.vector <- function(x, ...) invlog(x[1],x[2],x[3],...)

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
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return nothing but it does generatea plot
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg,outliers=TRUE,sitename="Middle Ground")
#' plot(ans)
plot.IL <- function(x, ...) {
  opar <- par(no.readonly=TRUE)
  on.exit(par(opar))
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
  xmax <- max(max(x$Lt,x$OutLt)*1.025,180)
  xmin <- min(min(x$Lt,x$OutLt) - 1,50)
  plot(x$Lt,x$DL,type="p",pch=20,xlab="",ylab="",xaxs="r",yaxs="r",
       xlim=c(xmin,xmax),ylim=c(-3,ymax))
  lines(x$PredLt,x$PredDL,col=2,lwd=2)
  lines(x$PredLt,x$PredDL+outer99,col=2,lty=2)
  lines(x$PredLt,x$PredDL-outer99,col=2,lty=2)
  lines(x$PredLt,x$PredDL+outer90,col=4,lty=2)
  lines(x$PredLt,x$PredDL-outer90,col=4,lty=2)

  abline(h=0,col="grey")
  abline(h=-3,col="grey")
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

  bins <- seq(xmin,xmax,5)
  hist(x$Lt,breaks=bins,xlab="",ylab="",main="")
  abline(v=x$L50,col="grey")
  title(ylab=list("Density of Data Points", cex=1.1, col=1, font=7))
  mtext("Initial Length Lt",side=1,line=0.5,outer=T,font=7,cex=1.25)
  if (nchar(x$sitename) > 0) label=x$sitename
     else label=paste0("Site ",x$siteid)
  mtext(label,side=3,line=0.5,outer=T,font=7,cex=1.25)
} # end of plot.IL


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
#' ans <- fitIL(midg)
#' print(ans)
print.IL <- function(x, ...) {
  cat("\n Site Id  : ",x$siteid)
  cat("\n sitename : ",x$sitename)
  cat("\n MaxDL    : ",x$MaxDL)
  cat("\n L50      : ",x$L50)
  cat("\n L95      : ",x$L95)
  cat("\n MaxSig   : ",x$MaxSig)
  cat("\n ")
}


#' @title summary.IL an S3 method for the summary generic for IL objects
#'
#' @description summary.IL is an S3 method for the summary generic for
#'     IL objects, as generated by the fitIL function.
#'
#' @param object an IL object generated by fitIL
#' @param ... the ellipsis is for any remaining parameters
#'
#' @return nothing but prints a summary ot the console
#' @export
#'
#' @examples
#' data(midg)
#' ans <- fitIL(midg)
#' summary(ans)
summary.IL <- function(object, ...) {
  Ltrge <- range(object$Lt,na.rm=TRUE)
  DLrge <- range(object$DL,na.rm=TRUE)
  outs <- FALSE
  if (length(object$MaxDLout) > 0) { outs <- TRUE }
  cat("\n siteid : ",object$siteid)
  cat("\n sitename: ",object$sitename)
  if (outs){ cat("\n MaxDL   : ",round(object$MaxDL,digits=4),"   ",
                 round(object$MaxDLout,digits=4)) }
  else  { cat("\n MaxDL   : ",round(object$MaxDL,digits=4)) }
  if (outs){ cat("\n L50     : ",round(object$L50,digits=4),"  ",
                 round(object$L50out,digits=4)) }
  else  { cat("\n L50     : ",round(object$L50,digits=4)) }
  if (outs){ cat("\n L95     : ",round(object$L95,digits=4)," ",
                 round(object$L95out,digits=4)) }
  else  { cat("\n L95     : ",round(object$L95,digits=4)) }
  if (outs){ cat("\n MaxSig  : ",round(object$MaxSig,digits=4),"   ",
                 round(object$MaxSigout,digits=4)) }
  else  { cat("\n MaxSig  : ",round(object$MaxSig,digits=4)) }
  cat("\n N       : ",object$Nobs)
  cat("\n Outliers: ",length(object$OutLt))
  cat("\n Range Lt: ",Ltrge)
  cat("\n Range DL: ",DLrge)
  cat("\n -ve LL  : ",object$model$minimum)
  cat("\n Other Components")
  cat("\n $Lt and $DL are the input data minus any outliers")
  cat("\n $model  contains the nlm fit")
  cat("\n $PredLt and PredDL = fitted line")
  cat("\n $OutLt and $OutDL = outlier values")
  cat("\n ")
}
