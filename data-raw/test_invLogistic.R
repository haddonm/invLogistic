

library(rutilsMH)
library(invLogistic)

options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)



data(midg)

ans <- fitIL(midg$Lt,midg$DL,outliers=TRUE,sitename="Middle Ground")

print(ans)

plotprep(width=7,height=6,newdev=FALSE)
plot(ans)

dyn <- getdyn(ans,maxage=30)

summary(ans)










