

library(codeutils)
library(hplot)
library(invLogistic)

options("show.signif.stars"=FALSE,"stringsAsFactors"=FALSE,
        "max.print"=50000,"width"=240)

ddir <- getDBdir()

rundir <- pathtopath(ddir,"/rcode2/invLogistic/data-raw/")


urch <- read.csv(file=pathtopath(rundir,"JS_data_clean.csv"),header=TRUE)

pick <- which(urch$Days_free > 0)
length(pick)

urch1 <- urch[pick,]

table(urch1$Site_Code)

pickS <- which(urch1$Site_Code == "3E")
label <- "E3"
used <- urch1[pickS,]
table(used$Days_free)


ans <- fitIL(used$Test_Diameter,used$Jaw_Growth_Increment,outliers=TRUE,sitename=label)

print(ans)

plotprep(width=7,height=6,newdev=FALSE)
plot(ans,miny=-0.5)


summary(ans)


boot <- dobootIL(ans,reps=100)

plotprep(width=7,height=6,newdev=FALSE)
plot(boot)




