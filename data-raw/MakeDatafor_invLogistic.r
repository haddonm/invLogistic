
library(rutilsMH)
library(invLogistic)
options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)

# abdat tagging data -------------------------------------------------
# From Middle Ground in the Acteaons, block 13, recaptured in 2003,
# Site 478, Latitude -43.54 long 146.99, 350 observations. All
# DeltaT = 1 year

datdir <- "C:/Users/Malcolm/Dropbox/rcode2/invLogistic/data-raw/"


midg <- read.csv(paste0(datdir,"midgrd.csv"),header=TRUE)
head(midg,20)

pick <- which(midg$DL < -4)
if (length(pick) > 0) midg[pick,]
midg <- midg[-pick,]

plotprep(width=7, height=5, newdev = FALSE)
plot(midg$Lt,midg$DL,type="p",pch=16,cex=1,xlim=c(5,180))
abline(h=0,col=1)

save(midg,file="data-raw/midg.RData")

# blackisland abalone tagging data ------------------------------

infile <- paste0(datdir,"blackisland.csv")
bi <- read.csv(infile,header=TRUE)
bi <- droplevels(bi[,-c(1,2,3,5,7)])
colnames(bi) <- c("dt","l1","l2","dl")

pick <- which(bi$dl < 0)
bi[pick,]
bi[pick[1],"l2"] <- 155; bi[pick[1],"dl"] <- 0
bi[pick[2],"l2"] <- 169; bi[pick[2],"dl"] <- 0
bi[pick,]

filename <- paste0(datdir,"bi.RData")
save(bi,file=filename)



# tasab maturity data -------------------------------------------

infile <- paste0(datdir,"tasw.csv")
tasab <- read.csv(infile,header=TRUE)
head(tasab)

tasab <- tasab[order(tasab$site,tasab$length),]

colnames(tasab) <- c("site","sex","length","mature")

filename <- paste0(datdir,"tasab.RData")
save(tasab,file=filename)



# check and re-store data ---------------------------------------
#
ddir <- "C:/Users/Malcolm/Dropbox/rcode2/invLogistic/data-raw"
tools::checkRdaFiles(paths=ddir)
tools::resaveRdaFiles(paths=ddir,compress="auto")
tools::checkRdaFiles(paths=ddir)













