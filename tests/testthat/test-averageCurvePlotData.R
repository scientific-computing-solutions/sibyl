source("setupFunctions.R")
context("AvCurvePlotData")

test_that("non_positive_finite_maxTime_gives_error",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse")
  expect_error(createAvCurvePlotData(fit,maxTime=-5))
  expect_error(createAvCurvePlotData(fit,maxTime=Inf))
  expect_error(createAvCurvePlotData(fit,maxTime=c(5,6)))
  expect_error(createAvCurvePlotData(fit,maxTime="notpositivenumber"))
})

test_that("non_positive_integer_Nsim_gives_error",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse")
  expect_error(createAvCurvePlotData(fit,Nsim=100.7))
  expect_error(createAvCurvePlotData(fit,Nsim=-5))
  expect_error(createAvCurvePlotData(fit,Nsim="notanumber"))
})

test_that("values_not_integers_greater_than_one_as_Npoints_argument_gives_error",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse")
  expect_error(createAvCurvePlotData(fit,Npoints=45.5))
  expect_error(createAvCurvePlotData(fit,Npoints=-5))
  expect_error(createAvCurvePlotData(fit,Npoints=0))
  expect_error(createAvCurvePlotData(fit,Npoints=1))
  expect_error(createAvCurvePlotData(fit,Npoints=c(4,5)))
})

test_that("error_if_requested_model_not_in_survival_model",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm"))
  expect_error(createAvCurvePlotData(fit,models=c("exponential")))
  
})
  
test_that("by_default_all_fitted_models_are_used_to_generate_output_data",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm"))
  avCurvePlotData <- createAvCurvePlotData(fit,Npoints=5,Nsim=20)
  
  expect_true(all(c("weibull","lnorm") %in% avCurvePlotData@plotdata$model))
  expect_true(all(avCurvePlotData@plotdata$model %in% c("weibull","lnorm","KM")))
  
})  

test_that("only_models_parameter_(if_not_NULL)_model_is_used_to_generate_averaged_data",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm","gompertz"))
  avCurvePlotData <- createAvCurvePlotData(fit,Npoints=5,Nsim=20,models=c("weibull","lnorm"))
  expect_true(all(c("weibull","lnorm") %in% avCurvePlotData@plotdata$model))
  expect_true(all(avCurvePlotData@plotdata$model %in% c("weibull","lnorm","KM")))
})

test_that("number_of_rows_per_modelArm_matches_Npoints_argument",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm","gompertz"))
  avCurvePlotData <- createAvCurvePlotData(fit,Npoints=5,Nsim=20,models=c("weibull","lnorm"))
  
  lnormData <- avCurvePlotData@plotdata[avCurvePlotData@plotdata$model=="lnorm",]
  
  expect_equal(nrow(lnormData),2*5) #2 arms, 5 points each
})

test_that("values_of_time_are_calculated_correctly_from_Npoints_and_maxTime_argument",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm","gompertz"))
  avCurvePlotData <- createAvCurvePlotData(fit,Npoints=5,Nsim=20,models=c("weibull","lnorm"),maxTime=100)
  
  lnormData <- avCurvePlotData@plotdata[avCurvePlotData@plotdata$model=="lnorm" & 
                                        avCurvePlotData@plotdata$Arm=="patchOnly",]
  
  expect_equal(lnormData$t, c(0,25,50,75,100))
})

test_that("maxTime_less_than_last_KM_time_uses_KMtime_as_maxTime",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model="weibull")
  avCurvePlotData <- createAvCurvePlotData(fit,Npoints=5,Nsim=20,maxTime=10)
  
  actualMaxTime <- max(avCurvePlotData@plotdata$t[avCurvePlotData@plotdata$model!="KM"])
  
  km <- survfit(Surv(ttr,!ttr.cens) ~ arm,survivalData@subject.data)
  
  expect_equal(actualMaxTime,max(km$time))
  expect_true(actualMaxTime>10)
  
})

test_that("null_maxTime_uses_last_KM_time_as_maxTime",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model="weibull")
  avCurvePlotData <- createAvCurvePlotData(fit,Npoints=5,Nsim=20)
  
  actualMaxTime <- max(avCurvePlotData@plotdata$t[avCurvePlotData@plotdata$model!="KM"])
  
  km <- survfit(Surv(ttr,!ttr.cens) ~ arm,survivalData@subject.data)
  
  expect_equal(actualMaxTime,max(km$time))
})

test_that("KM_part_of_avcurveplot_matches_outputfromsurvfit_when_covariates_are_used_to_fit_the_model",{
  
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model="weibull",covariates="race")
  avCurvePlotData <- createAvCurvePlotData(fit,Npoints=5,Nsim=20,maxTime=10)
  
  #No covariates here
  km <- survfit(Surv(ttr,!ttr.cens) ~ arm,survivalData@subject.data)

  
  kmPatchOnlyData <- avCurvePlotData@plotdata[avCurvePlotData@plotdata$model=="KM" &
                                              avCurvePlotData@plotdata$Arm== "patchOnly",]
  
  #Need a t=0, S=1 added to top of KM if it is not there 
  expectedTime <- km$time[1:km$strata[1]]
  expectedS <- km$surv[1:km$strata[1]]
  if(expectedTime[1] != 0){
    expectedTime <- c(0,expectedTime)
    expectedS <- c(1,expectedS)
  }
  
  expect_equal(kmPatchOnlyData$t,expectedTime) 
  expect_equal(kmPatchOnlyData$S,expectedS) 
  
})


