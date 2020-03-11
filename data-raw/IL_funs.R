

## calculate the mean
estGR <- function(x,tMaxDL,tL50,tL95) {
   pick <- which(x$Lt < x$L50)
   xLt <- x$Lt[pick]
   yDL <- x$DL[pick]
   PredDL <- invlog(c(tMaxDL,tL50,tL95),xLt)
   GR <- mean(yDL-PredDL)
   return(GR)
}


