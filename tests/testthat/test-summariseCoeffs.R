context("summariseCoeffs")

generateListOfModels <- function(){
  data("sibylData")
  
  weibull <- list(flexsurvreg(Surv(ttr,!ttr.cens)~grp,dist="weibull",data=sibylData))
  lnorm <- list(flexsurvreg(Surv(ttr,!ttr.cens)~grp,dist="lnorm",data=sibylData))
  gompertz <- list(
    patchOnly=flexsurvreg(Surv(ttr,!ttr.cens)~race,dist="gompertz",data=sibylData[sibylData$grp=="patchOnly",]),
    combination=flexsurvreg(Surv(ttr,!ttr.cens)~age,dist="gompertz",data=sibylData[sibylData$grp=="combination",])
  )
  
  models <- list(weibull=weibull,
                 gompertz=gompertz,
                 lnorm=lnorm)
  
}

test_that("output_is_list_of_lists_mirroring_the_shape_of_the_input_list",{
  models <- generateListOfModels() 

  dummyFn <- function(fittedModel){class(fittedModel)}
  
  results <- calcModelTables(models,summaryFn = function(modelName,modelClass){dummyFn},class="matrix")
  
  expect_equal(names(results),c("weibull","gompertz","lnorm"))
  expect_equal(length(results$weibull),1)
  expect_equal(length(results$weibull),1)
  expect_equal(length(results$gompertz),2)
  
})

test_that("vcov_displays_covariances",{
  models <- generateListOfModels() 
  
  #create dummy SurvivalModel object
  vcovFn <- function(fittedModel){vcov(fittedModel)}
  results <- calcModelTables(models,summaryFn = function(modelName,modelClass){vcovFn}, class="matrix")
  expect_equal(results$gompertz$patchOnly, vcov(models$gompertz$patchOnly) )
  
})

test_that("getParamSummariser_gives_error_for_unsupported_survreg_models",{
  expect_error(getParamSummariser("lnorm","survreg"))
})

test_that("getParamSummariser_outputs_parameter_summaries",{
  data("sibylData")
  
  weibull <- flexsurvreg(Surv(ttr,!ttr.cens)~grp,dist="weibull",data=sibylData)
  expect_equal(getParamSummariser("weibull","flexsurvreg")(weibull), weibull$res[, 1:3])
})