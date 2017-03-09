source("setupFunctions.R")
context("semiParametricFitting")

test_that("using_endpoint_not_in_SurvivalData_object_gives_error",{
  survivalData <- createSurvivalDataObject()

  expect_error(fitSemiParametric(survivalData,endPoint="nonsense"))

  # The defined endpoints are not vector-valued
  expect_error(fitSemiParametric(survivalData,endPoint=c("relapse","relapse")))
})

test_that("using_subgroup_not_in_SurvivalData_object_gives_error",{
  survivalData <- createSurvivalDataObject()
  expect_error(fitSemiParametric(survivalData,endPoint="relapse",subgroup="mysubgroup"))
})

test_that("error_if_an_arm_contains_no_data", {

  data("sibylData")

  for (s in c("patchOnly", "combination")){
    # Create subgroup that is just an indicator for arm membership, so that
    # subsetting by it results in no data in any other arm
    sibylData$sub.isMale <- sibylData$grp == s

    inputs <- survivalDataConstuctorTestSetUp()

    survivalData <- SurvivalData(data = sibylData,
                                 armDef = inputs[["arm"]],
                                 covDef = inputs[["cov"]],
                                 subgroupDef = inputs[["sub"]],
                                 subjectCol = "ID",
                                 endPointNames = c("relapse", "newEndpoint"),
                                 censorCol = c("ttr.cens", "cens.2"),
                                 timeCol = c("ttr", "end.2"))

    expect_error(fitSemiParametric(survivalData, endPoint="relapse", subgroup = "sub.isMale"))
  }
})

test_that("error_if_arm_has_no_events", {

  data("sibylData")

  for (a in c("patchOnly", "combination")){
    # Censor all subjects on one arm
    sibylData$ttr.cens <- sibylData$grp == a

    inputs <- survivalDataConstuctorTestSetUp()

    survivalData <- SurvivalData(data = sibylData,
                                 armDef = inputs[["arm"]],
                                 covDef = inputs[["cov"]],
                                 subgroupDef = inputs[["sub"]],
                                 subjectCol = "ID",
                                 endPointNames = c("relapse", "newEndpoint"),
                                 censorCol = c("ttr.cens", "cens.2"),
                                 timeCol = c("ttr", "end.2"))

    for (s in list(as.character(NA), "sub.isMale")){
      expect_error(fitSemiParametric(survivalData, endPoint="relapse", subgroup = s))
    }
  }
})

test_that("using_covariate_or_strata_not_in_SurvivalData_gives_error",{
  survivalData <- createSurvivalDataObject()
  expect_error(fitSemiParametric(survivalData,endPoint="relapse",covariates=c("age","otherCovar")))
  expect_error(fitSemiParametric(survivalData,endPoint="relapse",strata="otherCovar"))
})

test_that("invalid_conf.type_throws_error",{
  survivalData <- createSurvivalDataObject()
  
  expect_error(fitSemiParametric(survivalData,endPoint="relapse", conf.type="invalid"))
  
})

test_that("SemiParametricModelObjects_can_be_created_with_KM_and_Cox_fitted_approrpriately",{
  survivalData <- createSurvivalDataObject()
  
  sP <- fitSemiParametric(survivalData,endPoint="relapse")
  expect_equal(class(sP)[1],"SemiParametricModel")
  
  #km:
  km <- survfit(Surv(ttr,!ttr.cens) ~ arm, data=survivalData@subject.data)
  #set calls to be the same
  km$call <- sP@km$call
  
  expect_equal(sP@km, km)
  
  #Cox:
  cox <- coxph(Surv(ttr,!ttr.cens) ~ arm, data=survivalData@subject.data, ties="breslow", model=TRUE)
  cox$call <- sP@cox$call
  expect_equal(sP@cox, cox)
  
})

test_that("conf.type_argument_is_passed_to_survfit",{
  
  survivalData <- createSurvivalDataObject()
  
  sP <- fitSemiParametric(survivalData,endPoint="relapse", conf.type="log-log")
  #km:
  km <- survfit(Surv(ttr,!ttr.cens) ~ arm, data=survivalData@subject.data, conf.type="log-log")
  
  expect_equal(quantile(sP@km, prob=0.5, conf.int=TRUE), quantile(km, prob=0.5, conf.int=TRUE) )
  
  
})

test_that("SemiParametricModelObjects_can_be_created_with_covariates",{
  survivalData <- createSurvivalDataObject()
  
  sP <- fitSemiParametric(survivalData,endPoint="relapse",covariates=c("age","race"))

  #km:
  km <- survfit(Surv(ttr,!ttr.cens) ~ arm, data=survivalData@subject.data)
  #set calls to be the same
  km$call <- sP@km$call
  
  expect_equal(sP@km, km)
  
  #Cox:
  cox <- coxph(Surv(ttr,!ttr.cens) ~ arm+age+race, data=survivalData@subject.data, ties="breslow", model=TRUE)
  cox$call <- sP@coxWithStrata$call
  expect_equal(sP@coxWithStrata, cox)
  
})

test_that("SemiParametricModelObjects_can_be_created_with_subgroups_and_strata",{
  survivalData <- createSurvivalDataObject()
  
  sP <- fitSemiParametric(survivalData,endPoint="relapse",strata="race",subgroup="sub.isMale")
  
  df <- survivalData@subject.data[survivalData@subject.data$sub.isMale,]
  
  #km:
  km <- survfit(Surv(ttr,!ttr.cens) ~ arm, data=df)
  #set calls to be the same
  km$call <- sP@km$call
  
  expect_equal(sP@km, km)
  
  #Cox:
  cox <- coxph(Surv(ttr,!ttr.cens) ~ arm+strata(race), data=df, ties="breslow", model=TRUE)
  cox$call <- sP@coxWithStrata$call
  expect_equal(sP@coxWithStrata, cox)
})

test_that("only_appropriate_subgroup_data_is_added_to_survdata_slot",{
  survivalData <- createSurvivalDataObject()
  
  sP <- fitSemiParametric(survivalData,endPoint="relapse",strata="race",subgroup="sub.isMale")
  
  expect_true(all(sP@survData@subject.data$sub.isMale))
  expect_equal(nrow(sP@survData@subject.data),
               nrow(survivalData@subject.data[survivalData@subject.data$sub.isMale,]))
  
})

test_that("subjects_with_missing_endpoint_data_are_not_added_to_survdata_slot",{
  survivalData <- createSurvivalDataObject()
  survivalData@subject.data$ttr[1] <- NA
  survivalData@subject.data$ttr.cens[1] <- NA
  
  sP <- fitSemiParametric(survivalData,endPoint="relapse",strata="race")
  expect_equal(sP@survData@subject.data, survivalData@subject.data[2:nrow(survivalData@subject.data),])
  
})

test_that("all_data_is_added_to_survdata_slot_if_no_subgroups",{
  survivalData <- createSurvivalDataObject()
  
  sP <- fitSemiParametric(survivalData,endPoint="relapse")
  expect_equal(sP@survData,survivalData)
})

context("semiParametricFittingOutput")

test_that("logrank_test_matches_independentCoxFit_with_strata",{
  survivalData <- createSurvivalDataObject()
  sP <- fitSemiParametric(survivalData,endPoint="relapse",strata="race")
  logrankOutput <- coxphLogRankTest(sP)
  
  coxWithStrata <- coxph(Surv(ttr,!ttr.cens)~ arm + strata(race),
                         data=survivalData@subject.data)
  summStrata <- summary(coxWithStrata)[9]$logtest
  names(summStrata) <- NULL
  
  expect_equal(logrankOutput[2,1],summStrata[1])
  expect_equal(logrankOutput[2,2],summStrata[2])
  expect_equal(logrankOutput[2,3],summStrata[3])

})

test_that("logrank_test_with_no_strata_matches_even_strata_also_used",{
  survivalData <- createSurvivalDataObject()
  sP <- fitSemiParametric(survivalData,endPoint="relapse")
  logrankOutput <- coxphLogRankTest(sP)
  
  cox <- coxph(Surv(ttr,!ttr.cens)~ arm ,
                         data=survivalData@subject.data)
  summ <- summary(cox)[9]$logtest
  names(summ) <- NULL
  
  expect_equal(logrankOutput[1,1],summ[1])
  expect_equal(logrankOutput[1,2],summ[2])
  expect_equal(logrankOutput[1,3],summ[3])

})


test_that("number_of_events_is_correctly_calculated",{
  survivalData <- createSurvivalDataObject()
  sP <- fitSemiParametric(survivalData,endPoint="relapse",subgroup="sub.isMale")
  summarysP <- summary(sP, class="data.frame")
  
  subgroupData <- survivalData@subject.data[survivalData@subject.data$sub.isMale,]
  
  numEvents <- c(patchOnly=nrow(subgroupData[subgroupData$arm=="patchOnly" & !subgroupData$ttr.cens,]),
                 combination=nrow(subgroupData[subgroupData$arm=="combination" & !subgroupData$ttr.cens,]))
  
  expect_equal(summarysP[1,2:1],numEvents)
  
})


context("extractCumHazData")

test_that("outputs_one_dataframe_per_arm",{
  data("sibylData")
  
  km <- survfit(Surv(ttr,!ttr.cens)~grp, data=sibylData)  
  
  results <- extractCumHazData(km,armnames=c("B","A"))
  
  expect_equal(length(results),2)
  expect_true(is.data.frame(results[[1]]))
})

test_that("adds_given_armnames_to_output_dataframe",{
  data("sibylData")
  
  km <- survfit(Surv(ttr,!ttr.cens)~grp, data=sibylData)  
  
  results <- extractCumHazData(km,armnames=c("B","A"))
  
  expect_true(all(results[[1]]$Arm=="B"))
  expect_true(all(results[[2]]$Arm=="A"))
})


test_that("outputs_confidence_intervals_when_requested",{
  
  data("sibylData")
  km <- survfit(Surv(ttr,!ttr.cens)~grp, data=sibylData)  
  results <- extractCumHazData(km,armnames=c("B","A"),outputCI = TRUE)
  expect_equal(colnames(results[[1]]),c("t","S","Arm","lower","upper"))
  
})

test_that("t0_S1_row_is added_to_dataframes",{
  data("sibylData")
  km <- survfit(Surv(ttr,!ttr.cens)~grp, data=sibylData)  
  results <- extractCumHazData(km,armnames=c("B","A"))
  
  expect_equal(results[[1]][1,1],0) #t
  expect_equal(results[[1]][1,2],1) #S
})

