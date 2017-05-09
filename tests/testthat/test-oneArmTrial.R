context("oneArmTrials")

createOneArmSurvivalData <- function(){
  data("sibylData")
  df <- sibylData[sibylData$grp == "patchOnly",] 
  
  armDef <- ColumnDef(columnName = "grp",
                      type = "categorical",
                      categories = factor(c("patchOnly"),
                                          levels=c("patchOnly")))
  
  covariateDef <- list(
    ColumnDef(columnName = "age",
              type = "numeric",
              unit= "years"),
    ColumnDef(columnName = "sub.isHeavySmoker",
              type = "logical",
              displayName = "Heavy smoker"),
    ColumnDef(columnName = "race",
              type = "categorical",
              categories = factor(c("black", "hispanic", "other", "white"),
                                  levels=c("black", "white", "hispanic", "other"))))
  
  # Define columns corresponding to subgroups
  subgroupDef <- list(
    ColumnDef(columnName = "sub.isMale",
              type = "logical",
              displayName = "Male"))

  SurvivalData(data = df,
               armDef = armDef,
               covDef = covariateDef,
               subgroupDef = subgroupDef,
               subjectCol = "ID",
               endPointNames = c("relapse", "newEndpoint"),
               censorCol = c("ttr.cens", "cens.2"),
               timeCol = c("ttr", "end.2"))   
}


test_that("isSingleArm_is_TRUE_for_SurvivalData_objects_with_one_arm",{
  survivalData <- createOneArmSurvivalData()
  expect_true(isSingleArm(survivalData))
})

test_that("error_thrown_if_try_to_create_single_arm_object_if_data_set_contains_more_than_one_arm",{
  
  data("sibylData")
  
  #armDef has only one arm
  armDef <- ColumnDef(columnName = "grp",
                      type = "categorical",
                      categories = factor(c("patchOnly"),
                                          levels=c("patchOnly")))
  
  #but SibylData contains 2 arms
  expect_error(SurvivalData(data = sibylData,
               armDef = armDef,
               subjectCol = "ID",
               endPointNames = c("relapse", "newEndpoint"),
               censorCol = c("ttr.cens", "cens.2"),
               timeCol = c("ttr", "end.2")))
  
})

test_that("isSingleArm_is_TRUE_for_SemiParametricModel_objects_with_one_arm",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  expect_true(isSingleArm(sP))
})

test_that("km_is_calculate_correctly_for_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse", conf.type="plain")
  
  km <- survfit(Surv(ttr, !ttr.cens) ~ 1, data=sP@survData@subject.data, conf.type="plain")
  
  #change call
  km$call <- NULL
  sP@km$call <- NULL
  
  expect_equal(km, sP@km)
  
})

test_that("cox_slot_is_not_included_for_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  
  expect_true(is.null(sP@cox))
})

test_that("coxWithStrata_is_not_null_if_covariates_used",{
  
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse",
                          covariate="age")
  
  expect_false(is.null(sP@coxWithStrata))
})

test_that("error_thrown_if_calculating_logrank_test_for_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  
  expect_error(coxphLogRankTest(sP))
  
})

test_that("median_TTE_contains_only_one_arm_for_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  
  expect_equal(ncol(summary(sP, class="data.frame", type="medianTTE")), 1)
  
})

test_that("error_thrown_when_use.facet_is_true_for_diagnostic_plots_for_one_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  expect_error(plot(sP, type="CumHaz", use.facet=TRUE))
})


test_that("error_thrown_if_calculating_coxSnell_if_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  expect_error(coxSnellRes(sP))
})

test_that("error_thrown_if_calculating_schoenfeld_residuals_if_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  expect_error(schoenfeldResiduals(sP))
})

test_that("error_thrown_if_calculating_semi_parametric_restricted_mean_if_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  sP <- fitSemiParametric(survivalData, endPoint="relapse")
  expect_error(calcRmst(sP, alpha=0.05, truc=40, class="rmst"))
})


test_that("error_thrown_if_armAsFactor_TRUE_for_one_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  expect_error(fitModels(survivalData, armAsFactor=TRUE, endPoint="relapse"))
})

test_that("isSingleArm_is_true_for_SurvivalModel_for_single_arm_trial",{
  survivalData <- createOneArmSurvivalData()
  fits <- fitModels(survivalData, armAsFactor=FALSE, endPoint="relapse")
  expect_true(isSingleArm(fits))
})

test_that("restricted_mean_table_has_one_row_of_data",{
  survivalData <- createOneArmSurvivalData()
  fits <- fitModels(survivalData, armAsFactor=FALSE, endPoint="relapse")
  modelRmsts <- calcModelRmst(fits, model="exponential", times=c(0,4,Inf),
                              class="data.frame")
  
  expect_equal(nrow(modelRmsts), 2) #1 row for data, one row for heading
  
})

