source("setupFunctions.R")
context("LifeTable")

test_that("getLastPoint__extracts_the_data_from_the_largest_timepoint_less_than_input_parameter",{
  
  km <- data.frame(t=c(0,1,4,7,10),
                   S=c(1,0.8,0.6,0.4,0.3))
  
  expect_equal(getLastPoint(4,km),0.6)
  expect_equal(getLastPoint(0,km),1)
  expect_equal(getLastPoint(8,km),0.4)
})

test_that("getLastPoint_returns_NA_if_time_is_greater_than_largest_timepoint",{
  km <- data.frame(t=c(0,1,4,7,10),
                   S=c(1,0.8,0.6,0.4,0.3))
  expect_equal(getLastPoint(25,km),as.numeric(NA))
})

test_that("error_is_thrown_if_negative_time",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse")
  
  expect_error(createLifeTable(fit,times=c(4,8,-5)))
  expect_error(createLifeTable(fit,times=-5))
  
})

test_that("error_is_thrown_if_non_numeric_time",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse")
  
  expect_error(createLifeTable(fit,times=c(5,"text")))
})

test_that("error_is_thrown_if_model_requested_has_not_been_fit",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse")
  
  expect_error(createLifeTable(fit,times=c(5,10),model="invalidmodel"))
  expect_error(createLifeTable(fit,times=c(5,10),model="spline"))
})

test_that("getBestModel_uses_sum_of_AICs_of_models_when_armAsFactor_is_FALSE",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",armAsFactor=FALSE,model=c("weibull","gompertz"))
  
  icTable <- internalGetICTable(fit)
  
  if(sum(icTable$AIC[icTable$Model=="weibull"]) < sum(icTable$AIC[icTable$Model=="gompertz"])){
    bestModel <- "weibull"
  }else{
    bestModel <- "gompertz"
  }
  
  expect_equal(getBestModel(fit), bestModel)
  
})

test_that("getBestModel_uses_AIC_of_models_when_armAsFactor_is_TRUE",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",armAsFactor=TRUE,model=c("weibull","gompertz"))
  
  icTable <- internalGetICTable(fit)
  if(icTable$AIC[icTable$Model=="weibull"] < icTable$AIC[icTable$Model=="gompertz"]){
    bestModel <- "weibull"
  }else{
    bestModel <- "gompertz"
  }
  
  expect_equal(getBestModel(fit), bestModel)
  
})

test_that("nocovariate_and_armAsFactor_match_summary.flexsurvreg",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",armAsFactor=TRUE,model=c("weibull","gompertz"))
  lT <- createLifeTable(fit, times = seq(0, 20, 4), class="data.frame", model="gompertz")
  lifeTablesByArm <- split(lT,lT$Arm)
  expect_equal(lifeTablesByArm$patchOnly$t, seq(0,20,4))
  
  summ <- summary(fit@models$gompertz[[1]],t=seq(0,20,4),newdata=data.frame(arm="patchOnly"))
  expect_equal(lifeTablesByArm$patchOnly$gompertz,summ[[1]]$est)

})

test_that("nocovariate_and_no_armAsFactor_match_summary.flexsurvreg",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",armAsFactor=FALSE,model=c("weibull","gompertz"))
  lT <- createLifeTable(fit, times = seq(0, 20, 4), class="data.frame", model="gompertz")
  lifeTablesByArm <- split(lT,lT$Arm)
  expect_equal(lifeTablesByArm$patchOnly$t, seq(0,20,4))
  
  summ <- summary(fit@models$gompertz[[1]],t=seq(0,20,4))
  expect_equal(lifeTablesByArm$patchOnly$gompertz,summ[[1]]$est)
})
