#' ---
#' title: Prototyping R-package for Sibyl
#' author: Daniel Dalevi
#' date: 2016-09-14
#' ---


library( ggplot2 )
library( eventPrediction )
library( flexsurv )
library( plyr )
rm( list=ls() )
source( "R/survivalData.R")
source( "R/survivalModels.R")
source( "R/averageModelResult.R")

## PREPARE A DATASET FROM ASAUR PACKAGE 
library( asaur )
data( "pharmacoSmoking" )
pharmacoSmoking$ttr <- as.numeric( pharmacoSmoking$ttr )
pharmacoSmoking$censor <- as.integer( pharmacoSmoking$relapse==0 )
pharmacoSmoking$sub.gender <- pharmacoSmoking$gender=="Male"


## CREATE A SURVIVAL DATA OBJECT WHICH WILL BE THE BASIS FOR ALL ANALYSIS 
my.data <- SurvivalData(  pharmacoSmoking, 
                          arm="grp",
                          covariates=c( "age", "race" ),
                          subgroups=c( "sub.gender" ),
                          ctrl.arm="patchOnly",
                          active.arm="combination",
                          subject="id",
                          censor="censor",
                          time="ttr" )

## LOGLOG SURVIVAL CURVES FOR THE SUBGROUPS
plot( my.data )


## MAKE SUMMARY TABLE
summary( my.data )

# SURVIVAL FORMULA WITH COVARIATES 
f.surv <- Surv( time, has.event) ~ arm + age

# FIT STANDARD MODELS
models <- fitStandardModels( f.surv, my.data )

# ADD GENERALIZED GAMMA MODEL 
models <- c( models, fitGeneralizedGammaModel( f.surv, my.data ) )

# PRINT AIC TABLE
printICtable( models )


# PLOT AVERAGE SURVIVAL CURVE FOR GENERALIZED GAMMA MODEL
times <- seq( 0, 6 * 365.25 / 12, length.out=100 )
result <- averageModel( models$gengamma, my.data, times, 100 )
plot( result, main="Generalized Gamma Model" )



# ALSO PLOT THE AVERAGE SURVIVAL CURVE FOR GOMPERTZ
times <- seq( 0, 6 * 365.25 / 12, length.out=100 )
result <- averageModel( models$gompertz, my.data, times, 100 )
plot( result, main="Gompertz Model" )

