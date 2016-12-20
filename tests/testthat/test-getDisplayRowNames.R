context("getDisplayRowNames")

test_that("no_covariate_leaves_rownames_unchanged",{
  
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1))
  
  
  model <- flexsurvreg(Surv(time,!cens) ~ 1,data=df, dist="exponential")
  def <- list()
  
  result <- getDisplayRowNames(model,def)
  expect_equal(result,"rate")
})


test_that("numeric_covariate_output_as_displayName",{
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   cov=c(1,4.5,1,4,3,2.3,2,2.5))
  
  
  model <- flexsurvreg(Surv(time,!cens) ~ cov,data=df, dist="gompertz")
  def <- list(ColumnDef(columnName = "cov",displayName = "my name",type = "numeric"))
  
  result <- getDisplayRowNames(model,def)
  expect_equal(result, c("shape","rate","my name"))
})


test_that("logical_covariate_outputs_as_displayName:TRUE",{
  
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   cov=c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE))
  
  
  model <- flexsurvreg(Surv(time,!cens) ~ cov,data=df, dist="weibull")
  def <- list(ColumnDef(columnName = "cov",displayName = "my name",type = "logical"))
  
  result <- getDisplayRowNames(model,def)
  expect_equal(result, c("shape","scale","my name:TRUE"))
})


test_that("multiple_covariates_output_in_unchanged_order",{
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   cov1=c(TRUE,FALSE,TRUE,FALSE,TRUE,FALSE,FALSE,TRUE),
                   cov2=c(1,4.5,1,4,3,2.3,2,2.5))
  
  model <- flexsurvreg(Surv(time,!cens) ~ cov1 + cov2,data=df, dist="weibull")
  def <- list(ColumnDef(columnName = "cov1",displayName = "my name",type = "logical"),
              ColumnDef(columnName = "cov2",displayName = "my name2",type = "numeric"))
  
  result <- getDisplayRowNames(model,def)
  
  expect_equal(result,c("shape","scale","my name:TRUE","my name2"))
  
})

test_that("factor_variable_outputs_all_but_first_level",{
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   cov=c("A","B","A","A","C","A","C","C"))
  
  model <- flexsurvreg(Surv(time,!cens) ~ cov, data=df, dist="weibull")
  def <- list(ColumnDef(columnName = "cov",displayName = "the display",type = "categorical",
                        categories = factor(c("A","B","C"))))
                   
  result <- getDisplayRowNames(model,def)
  
  expect_equal(result, c("shape","scale","the display:B","the display:C"))
})


test_that("factor_variable_with_missing_level_does_not_output_this_level",{
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   cov=c("A","B","A","A","C","A","C","C"))
  
  model <- flexsurvreg(Surv(time,!cens) ~ cov, data=df, dist="weibull")
  def <- list(ColumnDef(columnName = "cov",displayName = "the display",type = "categorical",
                        categories = factor(c("A","B","C","D"))))
  
  result <- getDisplayRowNames(model,def)
  
  expect_equal(result, c("shape","scale","the display:B","the display:C"))
})

test_that("factor_variable_with_missing_first_level_does_not_output_this_level",{
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   cov=c("D","B","D","D","C","D","C","C"))
  
  model <- flexsurvreg(Surv(time,!cens) ~ cov, data=df, dist="weibull")
  def <- list(ColumnDef(columnName = "cov",displayName = "the display",type = "categorical",
                        categories = factor(c("A","B","C","D"))))
  
  result <- getDisplayRowNames(model,def)
  
  expect_equal(result, c("shape","scale","the display:C","the display:D"))
})


test_that("variable_called_arm_does_not_output_arm:value",{
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   arm=c("D","B","D","D","C","D","C","C"))
  
  model <- flexsurvreg(Surv(time,!cens) ~ arm, data=df, dist="weibull")
  def <- list(ColumnDef(columnName = "arm",displayName = "arm",type = "categorical",
                        categories = factor(c("A","B","C","D"))))
  
  result <- getDisplayRowNames(model,def)
  
  expect_equal(result, c("shape","scale","C","D"))
  
})

test_that("variable_called_arm_does_outputs_arm:value_if_keeparmDisplayName_is_true",{
  df <- data.frame(time=c(5,10,34,33,22,28,79,100),
                   cens=c(1,1,0,0,0,1,0,1),
                   arm=c("D","B","D","D","C","D","C","C"))
  
  model <- flexsurvreg(Surv(time,!cens) ~ arm, data=df, dist="weibull")
  def <- list(ColumnDef(columnName = "arm",displayName = "arm",type = "categorical",
                        categories = factor(c("A","B","C","D"))))
  
  result <- getDisplayRowNames(model,def, keeparmDisplayName=TRUE)
  
  expect_equal(result, c("shape","scale","arm:C","arm:D"))
  
})