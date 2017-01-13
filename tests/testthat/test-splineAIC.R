context("splineAIC")

test_that("error_produced_if_invalid_k_argument",{
  
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, endPoint="relapse", model="weibull")
  expect_error(createSplineAICTable(fit, k="s"))
  expect_error(createSplineAICTable(fit, k="NA"))
})


test_that("error_produced_if_invalid_scale_argument",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, endPoint="relapse", model="weibull")
  expect_error(createSplineAICTable(fit, scale=c("odds","invalid")))
})


test_that("table_has_correct_number_of_rows_and_columns",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, endPoint="relapse", model="weibull")
  table <- createSplineAICTable(fit,k=2:4,scale=c("odds"))
  expect_equal(table$numrow, 3)
  expect_equal(table$numcol, 2)
})

test_that("AIC_value_matches_individual_fit_when_armAsFactor_TRUE",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, endPoint="relapse", model="spline", armAsFactor=TRUE,
                   modelOptions=list(spline=list(k=2,scale="odds")))
  table <- createSplineAICTable(fit,k=2,scale=c("odds"), class="data.frame")
  expect_equal(table$AIC,createIcTable(fit,class="data.frame")$AIC)
})


test_that("AIC_value_matches_when_armAsFactor_FALSE",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, endPoint="relapse", model="spline", armAsFactor=FALSE,
                   modelOptions=list(spline=list(k=3,scale="odds")))
  table <- createSplineAICTable(fit, k=3, scale=c("odds"), class="data.frame")
  expect_equal(table$AIC,createIcTable(fit,class="data.frame")$AIC)
  
})
