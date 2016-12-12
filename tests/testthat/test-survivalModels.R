source("setupFunctions.R")

context("fitModels")

test_that("invalid_endPoint_throws_error",{
  survivalData <- createSurvivalDataObject()
  expect_error(fitModels(survivalData,endPoint="invalid"))
  expect_error(fitModels(survivalData,endPoint=c("relapse","newEndpoint")))
})

test_that("invalid_subgroups_throws_error",{
  survivalData <- createSurvivalDataObject()
  expect_error(fitModels(survivalData,endPoint="relapse",subgroup="Male"))
})

test_that("invalid_covariates_throws_error",{
  survivalData <- createSurvivalDataObject()
  expect_error(fitModels(survivalData,endPoint="relapse",covariates=c("nonsense","race")))  
})

test_that("invalid_model_argument_throws_error",{
  survivalData <- createSurvivalDataObject()
  expect_error(fitModels(survivalData,endPoint="relapse",model=c("wrong","exp")))   
})

test_that("survival_model_not_recognised_if_preferredPacakage_is_flexSurv",{
  survivalData <- createSurvivalDataObject()
  expect_error(fitModels(survivalData,endPoint="relapse",model=c("loglogistic"))) 
})

test_that("warning_if_flexsurv_model_used_when_preferredPacakage_is_survival",{
  survivalData <- createSurvivalDataObject()
  expect_warning(fitModels(survivalData,endPoint="relapse",model=c("llogis"),preferredPackage="survival")) 
})

test_that("no_subgroup_argument_fits_models_to_whole_dataset",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull"))
  
  #survivalData object is unchanged
  expect_equal(fit@survData,survivalData)
  expect_true(is.na(fit@subgroup))
  
  #Weibull model fit to whole data set
  expect_equal(fit@models$weibull[[1]]$N,nrow(survivalData@subject.data))
  
})

test_that("subgroup_argument_fits_models_to_just_subgroup",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",subgroup="sub.isMale",model="weibull")
  
  expect_equal(fit@subgroup,"sub.isMale")
  expect_equal(fit@models$weibull[[1]]$N, sum(survivalData@subject.data$sub.isMale))
  expect_equal(nrow(fit@survData@subject.data), sum(survivalData@subject.data$sub.isMale))
  
  
})

test_that("subjects_with_missing_endpoint_data_are_removed_when_fitting_models",{
  
  survivalData <- createSurvivalDataObject()
  survivalData@subject.data$ttr[1] <- NA
  survivalData@subject.data$ttr.cens[1] <- NA
  
  fit <- fitModels(survivalData,endPoint="relapse",model="weibull")
  expect_equal(fit@survData@subject.data, survivalData@subject.data[2:nrow(survivalData@subject.data),])
  
})


test_that("requested_models_are_fitted_by_arm_when_armAsFactor_is_FALSE",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm"),armAsFactor=FALSE)
  
  expect_equal(names(fit@models),c("weibull","lnorm"))
  #two separate models
  expect_equal(names(fit@models$lnorm),c("patchOnly","combination"))
  
})

test_that("requested_models_are_fitted_on_while_dataSet_when_armAsFactor_it_TRUE",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",covariates="race",model=c("weibull","gompertz"),armAsFactor=TRUE)
  
  expect_equal(names(fit@models),c("weibull","gompertz"))
  #one models
  expect_equal(length(fit@models$gompertz),1)

  expectedModel <- flexsurvreg(Surv(ttr,!ttr.cens)~arm+race,data=survivalData@subject.data,dist="gompertz")
  expectedModel$call <- fit@models$gompertz[[1]]$call  
  
  
  expect_equal(fit@models$gompertz[[1]] , expectedModel)
})

test_that("slots_of_SurvivalModel_object_match_given_inputs_when_no_subgroup",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",covariates="race",model=c("weibull","gompertz"),armAsFactor=TRUE)
  
  expect_equal(fit@covariates,"race")
  expect_equal(fit@armAsFactor,TRUE)
  expect_equal(fit@survData,survivalData) #No subgroup so survivalData is not subsetted
  expect_equal(fit@endPointDef,list(timeCol="ttr",censorCol="ttr.cens"))
  expect_true(is.na(fit@subgroup))
  
})

test_that("slots_of_SurvivalModel_object_match_given_inputs_when_using_subgroup",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",covariates=c("race","age"),
                   model=c("weibull","gompertz"),armAsFactor=FALSE,subgroup="sub.isMale")
  
  expect_equal(fit@covariates,c("race","age"))
  expect_equal(fit@armAsFactor,FALSE)
  
  #keep only subjects in subgroup in survivalData object
  survivalData@subject.data <- survivalData@subject.data[survivalData@subject.data$sub.isMale,]
  expect_equal(fit@survData,survivalData)
  expect_equal(fit@endPointDef,list(timeCol="ttr",censorCol="ttr.cens"))
  expect_equal(fit@subgroup,"sub.isMale")
  
})


test_that("error_if_all_times_zero",{
  data("sibylData")
  
  sibylData$ttr <- 0
  input <- survivalDataConstuctorTestSetUp()
  
  survivalData <- SurvivalData(data = sibylData,
               armDef = input$arm,
               subjectCol = "ID",
               endPointNames = c("relapse", "newEndpoint"),
               censorCol = c("ttr.cens", "cens.2"),
               timeCol = c("ttr", "end.2"))
  
  expect_error(fitModels(survivalData,endPoint="relapse"))
  
})

test_that("error_if_all_times_in_one_arm_are_zero",{
  data("sibylData")
  
  sibylData$ttr[sibylData$grp=="patchOnly"] <- 0
  input <- survivalDataConstuctorTestSetUp()
  
  survivalData <- SurvivalData(data = sibylData,
                               armDef = input$arm,
                               subjectCol = "ID",
                               endPointNames = c("relapse", "newEndpoint"),
                               censorCol = c("ttr.cens", "cens.2"),
                               timeCol = c("ttr", "end.2"))
  
  expect_warning(expect_error(fitModels(survivalData,endPoint="relapse",armAsFactor=FALSE)))
  expect_warning(expect_error(fitModels(survivalData,endPoint="relapse",armAsFactor=TRUE)))
})

test_that("error_if_no_events",{
  data("sibylData")
  
  sibylData$ttr.cens <- TRUE
  input <- survivalDataConstuctorTestSetUp()
  
  survivalData <- SurvivalData(data = sibylData,
                               armDef = input$arm,
                               subjectCol = "ID",
                               endPointNames = c("relapse", "newEndpoint"),
                               censorCol = c("ttr.cens", "cens.2"),
                               timeCol = c("ttr", "end.2"))
  
  expect_error(fitModels(survivalData,endPoint="relapse",armAsFactor=FALSE))
  expect_error(fitModels(survivalData,endPoint="relapse",armAsFactor=TRUE))
  
})

test_that("error_if_no_events_on_one_arm",{
  data("sibylData")
  
  sibylData$ttr.cens[sibylData$grp=="patchOnly"] <- TRUE
  input <- survivalDataConstuctorTestSetUp()
  
  survivalData <- SurvivalData(data = sibylData,
                               armDef = input$arm,
                               subjectCol = "ID",
                               endPointNames = c("relapse", "newEndpoint"),
                               censorCol = c("ttr.cens", "cens.2"),
                               timeCol = c("ttr", "end.2"))
  
  expect_error(fitModels(survivalData,endPoint="relapse",armAsFactor=FALSE))
  expect_error(fitModels(survivalData,endPoint="relapse",armAsFactor=TRUE))
  
})

test_that("adding_a_model_already_included_throws_a_warning",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm"),armAsFactor=FALSE)
  
  expect_warning(addModel(fit,c("weibull")))
})


test_that("adding_a_model_already_included_overwrites_old_model",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","spline"),
                   modelOptions=list(spline=list(k=3)),armAsFactor=FALSE)
  
  #check spline has k=3 
  expect_equal(fit@models$spline$patchOnly$k, 3)
  
  expect_warning(fit <- addModel(fit,"spline",list(spline=list(k=4))))
  
  #check spline now has k=4 
  expect_equal(fit@models$spline$patchOnly$k, 4)
  
})



test_that("adding_a_new_model_adds_it_into_the_SurvivalModel_object",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm"),armAsFactor=FALSE)
  
  addedFit <- addModel(fit,c("gengamma","exponential"))
  expect_equal(names(addedFit@models),c("weibull","lnorm","gengamma","exponential"))
})

test_that("removing_an_existing_model_removes_it_from_the_fit",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm","gompertz"),armAsFactor=FALSE)
  
  removedFit <- removeModel(fit,c("lnorm","gompertz"))
  expect_equal(names(removedFit@models),c("weibull"))
})

test_that("attempting_to_remove_a_model_which_does_not_exist_throws_a_warning",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm","gompertz"),armAsFactor=FALSE)
  
  expect_warning(removeModel(fit,"gengamma"))
  
})

test_that("attempting_to_remove_all_models_throws_error",{
  survivalData <- createSurvivalDataObject()
  fit <- fitModels(survivalData,endPoint="relapse",model=c("weibull","lnorm"))
  expect_error(removeModel(fir,c("weibull","lnorm")))
})