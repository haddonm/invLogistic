




dobootIL <- function(x,reps=100) {
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
   doboot <- function(inpick,inLt,inDL) {
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
}

plot.bootIL <- function(x) {
   opar <- par(no.readonly=TRUE)
   collab = c("MaxDL","L50","L95","MaxSig")
   initcol <- -2
   sitenum <- rownames(x$Percentiles)
   par(mfrow = c(2,2))
   par(mai=c(0.5,0.5,0.1,0.1), oma=c(0,0,2,0))
   par(cex=0.8, mgp=c(1.35,0.35,0), font.axis=7)
   for (index in 1:4) {
      initcol <- initcol + 4
      CIs <- x$Percentiles[initcol:(initcol+2)]
      hist(x$Replicates[,index],xlab="",main="")
      title(xlab=list(collab[index], cex=1.0, col=1, font=7))
      abline(v=CIs[1],col=2,lty=2)
      abline(v=CIs[2],col=4,lty=2)
      abline(v=CIs[3],col=2,lty=2)
   }
   mtext(paste("Site ",sitenum,sep=""),side=3,line=0.5,outer=T,font=7,cex=1.25)
   par(opar)
}



## calculate the mean
estGR <- function(x,tMaxDL,tL50,tL95) {
   pick <- which(x$Lt < x$L50)
   xLt <- x$Lt[pick]
   yDL <- x$DL[pick]
   PredDL <- invlog(tMaxDL,tL50,tL95,xLt)
   GR <- mean(yDL-PredDL)
   return(GR)
}


## x <- model
plotsingleIL <- function(x,removeoutlier=T) {
   expDL <-  invlog(c(x$MaxDL,x$L50,x$L95),x$Lt)
   expSD <- invlog(c(x$MaxSig,x$L95,210),x$PredLt)
   outer99 <- 2.5760 * expSD
   outer90 <- 1.965 * expSD
   par(mfrow = c(1,1))
   par(mai=c(0.4,0.4,0.1,0.1), oma=c(0,0,0,0))
   par(cex=0.8, mgp=c(1.35,0.35,0), font.axis=7)
 # Plot the basic fit with outliers if any
   ymax <- max(x$DL,x$OutDL)*1.025
   xmax <- max(max(x$Lt,x$OutLt)*1.025,180)
   xmin <- min(min(x$Lt,x$OutLt) - 1,50)
   plot(x$Lt,x$DL,type="p",pch=20,xlab="",ylab="",xaxs="r",yaxs="r",
        xlim<- c(xmin,xmax),ylim=c(-3,ymax))
   lines(x$PredLt,x$PredDL,col=2,lwd=2)
 #  lines(x$PredLt,x$PredDL+outer99,col=2,lty=2)
 #  lines(x$PredLt,x$PredDL-outer99,col=2,lty=2)
 #  lines(x$PredLt,x$PredDL+outer90,col=4,lty=2)
 #  lines(x$PredLt,x$PredDL-outer90,col=4,lty=2)
   abline(h=0,col="grey")
   abline(h=-3,col="grey")
   if ((length(x$OutLt)>0) & (removeoutlier == F)) {
       points(x$OutLt,x$OutDL,col=2,pch=20)
   }
   title(ylab=list("Growth Increment DL", cex=1.0, col=1, font=7),
         xlab=list("Shell Length (mm)", cex=1.0, col=1, font=7))
}

