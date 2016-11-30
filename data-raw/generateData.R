#script to generate data for package

set.seed(45)
N <- 200

#create data frame
sibylData <- data.frame(ID=1:N)

#Add covariates

sibylData$age <- rbinom(N,100,0.6)
sibylData$race <- factor(c("black","hispanic", "white", "other")[1+rbinom(N,3,0.5)])


#Add treatment group
sibylData$grp <- factor(c("patchOnly","combination")[1+rbinom(N,1,0.5)])

#Add subgroups
sibylData$sub.isMale <- as.logical(rbinom(N,1,0.4))
sibylData$sub.isHeavySmoker <- as.logical(rbinom(N,1,0.7))

#Add 2 endpoints
my.scale <- 1 + 3*as.numeric(sibylData$grp)+ sibylData$age/10
sibylData$ttr <- vapply(my.scale,function(x){rweibull(1,1,my.scale)},FUN.VALUE = numeric(1))
sibylData$ttr.cens <- rbinom(N,1,0.6)
#rands <- 0.7+0.3*runif(N)
#sibylData$ttr <- ifelse(sibylData$ttr.cens, sibylData$ttr*rands, sibylData$ttr)


my.scale <- 1 + 0.6*as.numeric(sibylData$grp)+ 0.1*as.numeric(sibylData$race)
sibylData$end.2 <- 40*vapply(my.scale,function(x){rweibull(1,1,my.scale)},FUN.VALUE = numeric(1))
sibylData$cens.2 <- rbinom(N,1,0.6)
#rands <- 0.7+0.3*runif(N)
#sibylData$end.2 <- ifelse(sibylData$cens.2, sibylData$end.2*rands, sibylData$end.2)


#To save data
devtools::use_data(sibylData, overwrite=TRUE)

