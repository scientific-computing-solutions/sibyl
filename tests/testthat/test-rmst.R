context("parametricRmst")

test_that("error_thrown_if_armAsFactor_TRUE",{
  
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=TRUE, endPoint="relapse")
  expect_error(calcModelRmst(fit, model="weibull", times=50))
  
})

test_that("error_thrown_if_fit_contains_covariates",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=FALSE, endPoint="relapse", covariate="age")
  expect_error(calcModelRmst(fit, model="weibull", times=50))
  
})

test_that("error_thrown_if_invalid_model",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=FALSE, endPoint="relapse",model="weibull")
  expect_error(calcModelRmst(fit, model="exponential", times=50))
})

test_that("error_thrown_if_invalid_time",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=FALSE, endPoint="relapse",model="weibull")
  expect_error(calcModelRmst(fit, model="weibull", times=c("dsfwe",50)))
})

test_that("one_column_is_generated_per_time_requested",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=FALSE, endPoint="relapse",model="weibull")
  tab <- calcModelRmst(fit, model="weibull", times=c(67,Inf, 10))
  expect_equal(tab$numcol,4) #1 for each time + 1 for arm names
})

test_that("one_row_is_generated_per_arm",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=FALSE, endPoint="relapse",model="weibull")
  tab <- calcModelRmst(fit, model="weibull", times=c(67,Inf, 10))
  expect_equal(tab$numrow, 4) # 1 for difference, 1 for "time" + 2 for arm
})

test_that("time_inf_gives_mean_of_distribution",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=FALSE, endPoint="relapse")
  tab <- calcModelRmst(fit, model="exponential", times=Inf, class="data.frame")
  expect_equal(tab[[1]][2], 1/fit@models$exponential[[1]]$res[1,"est"])
  expect_equal(tab[[1]][3], 1/fit@models$exponential[[2]]$res[1,"est"])
})

test_that("difference_row_is_difference_of_values_of_each_arm",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, armAsFactor=FALSE, endPoint="relapse")
  tab <- calcModelRmst(fit, model="llogis", times=c(5,10,20), class="data.frame") 
  expect_equal(as.numeric(tab[4,]),as.numeric(tab[3,])-as.numeric(tab[2,]))
})

