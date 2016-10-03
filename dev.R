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
source( "R/covariates.R")
source( "R/survivalData.R")
source( "R/survivalModels.R")
source( "R/averageModelResult.R")
source( "R/common.R")

## PREPARE A DATASET FROM ASAUR PACKAGE 
library( asaur )
data( "pharmacoSmoking" )
pharmacoSmoking$ttr <- as.numeric( pharmacoSmoking$ttr )
pharmacoSmoking$censor <- as.integer( pharmacoSmoking$relapse==0 )
pharmacoSmoking$sub.gender <- pharmacoSmoking$gender=="Male"

def.covariates = list(
  Covariate( name = "age",
           type = "numeric",
           unit= "years" ),
  Covariate( name = "race",
             type = "categorical", 
             categories = c( "black", "hispanic", "other", "white" ),
             unit= "" ) )


## CREATE A SURVIVAL DATA OBJECT WHICH WILL BE THE BASIS FOR ALL ANALYSIS 
my.data <- SurvivalData(  pharmacoSmoking, 
                          arm="grp",
                          covariates=c( "age", "race" ),
                          covdef = def.covariates,
                          subgroups=c( "sub.gender" ),
                          ctrl.arm="patchOnly",
                          active.arm="combination",
                          subject="id",
                          censor="censor",
                          time="ttr" )

## LOGLOG SURVIVAL CURVES FOR THE SUBGROUPS
plot( my.data, separate.plots=TRUE  )


## MAKE SUMMARY TABLE
summary( my.data )

# SURVIVAL FORMULA WITH COVARIATE AGE AND ARM AS A FACTOR IN MODEL
f.witharm <- survivalFormula( my.data, armAsFactor = TRUE, covariates="age" )
print( f.witharm )
# SURVIVAL FORMULA WITH NO COVARIATES, FIT SEPARATE MODELS FOR EACH ARM ("BY ARM")
f.byarm   <- survivalFormula( my.data, armAsFactor = FALSE )  
print( f.byarm )

# FIT STANDARD MODELS
models.witharm <- fitStandardModels( f.witharm, my.data, subgroups=c( "sub.gender" ) )
models.byarm   <- fitStandardModels( f.byarm, my.data, by="arm", subgroups=c( "sub.gender" ) )

# ADD AND FIT GENERALIZED GAMMA MODEL 
models.witharm <- addModel( models.witharm,  "gengamma" )
models.byarm   <- addModel( models.byarm,    "gengamma" )

# PRINT AIC/BIC TABLE
printICtable( models.witharm )
printICtable( models.byarm )

# PRINT COEFFICIENTS WITH ARM AS INTERCEPT IN MODEL 
printCoefficientEstimates( models.witharm )
# PRINT COEFFICIENTS WITH SEPARATE MODELS PER ARM 
printCoefficientEstimates( models.byarm )


# PLOT MODELS WITH ARM AS INTERCEPT IN MODEL 
plot( models.witharm )
# PLOT MODELS WITH SEPARATE MODELS PER ARM 
plot( models.byarm  )




# PLOT AVERAGE SURVIVAL CURVE FOR GENERALIZED GAMMA MODEL 
times <- seq( 0, 6 * 365.25 / 12, length.out=100 )
result <- averageModel( models.witharm@models$gengamma, my.data, times, 100 )
plot( result, main="Generalized Gamma Model" )


# ALSO PLOT THE AVERAGE SURVIVAL CURVE FOR GOMPERTZ
times <- seq( 0, 6 * 365.25 / 12, length.out=100 )
result <- averageModel( models.witharm@models$gompertz, my.data, times, 100 )
plot( result, main="Gompertz Model" )

