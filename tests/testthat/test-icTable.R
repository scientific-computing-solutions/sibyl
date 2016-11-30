source("setupFunctions.R")
context("icTable")

test_that("internalGetICTable_produces_data_frame_with_appropriate_rows_and_columns",{
  
  sData <- createSurvivalDataObject()
  models <- fitModels(sData,endPoint="relapse",armAsFactor=FALSE)
  
  results <- internalGetICTable(models)
  
  expect_equal(class(results),"data.frame")
  expect_equal(colnames(results),c("Model","Arm","AIC","BIC"))
  
  #one model per arm when armAsFactor is false
  expect_equal(nrow(results),2*length(models@models))
  
  #two models per arm when armAsFactor is true
  models <- fitModels(sData,endPoint="relapse",armAsFactor=TRUE)
  expect_equal(nrow(internalGetICTable(models)),length(models@models))
})

test_that("internalGetICTable_extracts_AIC_and_BIC_from_flexsurvreg",{
  sData <- createSurvivalDataObject()
  models <- fitModels(sData,endPoint="relapse",armAsFactor=TRUE)
  
  results <- internalGetICTable(models) 
  gompertzAIC <- vapply( models@models$gompertz,function(x){akaike_ic(x)}, FUN.VALUE = numeric(1) )
  lnormBIC <- vapply( models@models$lnorm,function(x){bayes_ic(x)}, FUN.VALUE = numeric(1) )
  
  names(gompertzAIC) <- NULL
  names(lnormBIC) <- NULL
  
  expect_equal(results$AIC[results$Model=="gompertz"],gompertzAIC) 
  expect_equal(results$BIC[results$Model=="lnorm"],lnormBIC) 
  
})

test_that("internalGetICTable_extracts_AIC_from_survreg",{
  sData <- createSurvivalDataObject()
  models <- fitModels(sData,endPoint="relapse",armAsFactor=FALSE,model=c("weibull", "lognormal"),preferredPackage="survival")

  results <- internalGetICTable(models) 
  expect_equal(results$BIC,rep(NA,4))
  
  weibullAIC <- vapply( models@models$weibull,function(x){extractAIC(x)[2]}, FUN.VALUE = numeric(1) )
  
  names(weibullAIC) <- NULL
  expect_equal(results$AIC[results$Model=="weibull"],weibullAIC)  
  
})


test_that("error_if_invalid_summaryFn_argument",{
  sData <- createSurvivalDataObject()
  models <- fitModels(sData,endPoint="relapse",armAsFactor=FALSE)
  expect_error(createIcTable(models,summaryFn="product"))
})

test_that("armAsFactor_as_true_displays_single_AIC_BIC_per_model_sorted_byAIC",{
  sData <- createSurvivalDataObject()
  models <- fitModels(sData,endPoint="relapse",armAsFactor=TRUE,model=c("exponential","weibull"))
  
  results <-  createIcTable(models,class="data.frame")
  expect_equal(nrow(results),2)
  expect_equal(ncol(results),4) #Model, Arm=All, AIC, BIC
  
  #sorted
  theAICs <- results[1:2,3]
  expect_equal(theAICs, sort(theAICs))
  
})

test_that("armAsFactor_false_summaryFn_as_sum_displays_single_AIC_BIC_per_model_sorted_byAIC",{
  sData <- createSurvivalDataObject()
  models <- fitModels(sData,endPoint="relapse",armAsFactor=FALSE)
  
  results <-  createIcTable(models,summaryFn="sum",class="data.frame")
  expect_equal(nrow(results),5)
  expect_equal(ncol(results),4)
  
  #sorted
  theAICs <- results[1:5,3]
  expect_equal(theAICs, sort(theAICs))
})

test_that("armAsFactor_false_summaryFn_as_identity_displays__AIC_BIC_per_arm",{
  sData <- createSurvivalDataObject()
  models <- fitModels(sData,endPoint="relapse",armAsFactor=FALSE)
  
  results <-  createIcTable(models,summaryFn="identity",class="data.frame")
  expect_equal(nrow(results),10) #2 for each arm
  expect_equal(ncol(results),4) #model, Arm, AIC, BIC
  
})