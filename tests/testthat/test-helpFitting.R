context("survivalFormula")

test_that("Surv_object_takes_time_and_negated_censor_column_as_arguments",{
  
  sF <- survivalFormula(armAsFactor=FALSE,
                        timeCol="my.time",
                        censorCol="my.cens")
  
  expect_equal(formula("Surv(my.time, !my.cens) ~ 1"),sF)
})

test_that("armAsFactor_true_adds_arm_to_rhs_of_formula",{
  sF <- survivalFormula(armAsFactor=TRUE,
                        timeCol="my.time",
                        censorCol="my.cens")
  
  expect_equal(formula("Surv(my.time, !my.cens) ~ arm"),sF)
  
})

test_that("covariates__added_to_rhs_of_formula",{
  sF <- survivalFormula(armAsFactor=TRUE,
                        timeCol="my.time",
                        censorCol="my.cens",
                        covariates="cov1")
  
  expect_equal(formula("Surv(my.time, !my.cens) ~ arm + cov1"),sF)
  
  sF <- survivalFormula(armAsFactor=FALSE,
                        timeCol="my.time",
                        censorCol="my.cens",
                        covariates=c("cov1","cov2"))
  
  expect_equal(formula("Surv(my.time, !my.cens) ~ cov1 + cov2"),sF)
  
})

test_that("strata_are_added_to_rhs_of_formula",{
  sF <- survivalFormula(armAsFactor=TRUE,
                        timeCol="my.time",
                        censorCol="my.cens",
                        strata="strata1")
  
  expect_equal(formula("Surv(my.time, !my.cens) ~ arm + strata(strata1)"),sF)
  
  sF <- survivalFormula(armAsFactor=FALSE,
                        timeCol="the.time",
                        censorCol="the.cens",
                        strata=c("strata1","strata2"),
                        covariates=c("cov1","cov2"))
  
  expect_equal(formula("Surv(the.time, !the.cens) ~ cov1 + cov2 + strata(strata1)+strata(strata2)"),sF)
  
})

context("helpFitting")

test_that("removeNonPositiveTimes_removes_non_positive_times",{
  
  df <- data.frame(time=c(-6,0,4),
                  cens=c(FALSE,FALSE,TRUE))
  
  expect_warning(outputData <- removeNonPositiveTimes(df,"time"))
  
  
  expecteddf <- data.frame(time=4,cens=TRUE)
  rownames(expecteddf) <- NULL
  rownames(outputData) <- NULL
  
  expect_equal(outputData,expecteddf)
  
})

test_that("removeNonPositiveTimes_keeps_data_frame_unchanged_if_all_times_positive",{
  df <- data.frame(time=c(-6,0,4),
                   time2=c(5,6,7),
                   cens=c(FALSE,FALSE,TRUE))
  
  expect_equal(removeNonPositiveTimes(df,"time2"),df)
  
})

test_that("extract_subgroup_extracts_given_subgroup",{
  
  df <- data.frame(subgroup=c(TRUE,TRUE,FALSE,FALSE,TRUE),
                   ID=c(4,7,1,3,6))
  
  
  extractedDf <- extractSubgroup(df,"subgroup")
  
  expect_equal(extractedDf$ID,c(4,7,6))
  expect_equal(extractedDf$subgroup,c(TRUE,TRUE,TRUE))
  
})

test_that("extract_subgroup_leaves_data_unchanged_if_no_subgroup",{
  df <- data.frame(subgroup=c(TRUE,TRUE,FALSE,FALSE,TRUE),
                   ID=c(4,7,1,3,6))
  
  
  extractedDf <- extractSubgroup(df,subgroup=NA)
  expect_equal(df,extractedDf)
  
})