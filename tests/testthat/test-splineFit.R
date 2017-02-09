context("fitSplines")

test_that("error_produced_if_invalid_k_argument",{
  
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, endPoint="relapse", model="weibull")
  expect_error(fitSplines(fit, k="s"))
  expect_error(fitSplines(fit, k="NA"))
})


test_that("error_produced_if_invalid_scale_argument",{
  sD <- createSurvivalDataObject()
  fit <- fitModels(sD, endPoint="relapse", model="weibull")
  expect_error(fitSplines(fit, scale=c("odds","hazard")))
})

