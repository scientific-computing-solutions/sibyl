source("setupFunctions.R")

# Test helper function - invalidate a column of the data and check for an error
errorIfColumnInvalid <- function(df, inputs, colName, invalidValues, condition=expect_error,
                                 colName2=NULL,invalidValues2=NULL ){

  for (x in seq_along(invalidValues)){
    testDf <- df
    testDf[, colName] <- invalidValues[[x]]
    
    if(!is.null(colName2)){
      testDf[,colName2] <- invalidValues2[[x]]
    }
    
    condition(SurvivalData(data=testDf,
                           armDef=inputs$arm,
                           subjectCol="ID",
                           covDef = inputs$cov,
                           subgroupDef = inputs$sub,
                           endPointNames="relapse",
                           censorCol="ttr.cens",
                           timeCol="ttr"))
  }
}


context("SurvivalData_constructor")

test_that("can_create_SurvivalData_objects_with_no_covariates",{
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()

  survivalData <- SurvivalData(data=sibylData,
               armDef=inputs$arm,
               subjectCol="ID",
               subgroupDef=inputs$sub,
               endPointNames="relapse",
               censorCol="ttr.cens",
               timeCol="ttr")

  expect_equal(survivalData@covDef,list())

})

test_that("can_create_SurvivalData_objects_with_no_subgroups",{
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()

  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               covDef=inputs$cov,
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr")

  expect_equal(survivalData@subgroupDef,list())
})

test_that("can_create_SurvivalData_with_a_subject_with_NAs_for_endpoint_data",{
  data("sibylData")
  
  inputs <- survivalDataConstuctorTestSetUp()
  
  sibylData$ttr[1] <- NA
  sibylData$ttr.cens[1] <- NA
  
  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               covDef=inputs$cov,
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr")
  
  expect_equal(class(survivalData@subject.data$ttr), "numeric")
  expect_equal(class(survivalData@subject.data$ttr.cens), "logical")
  
  expect_true(is.na(survivalData@subject.data$ttr[1]))
  expect_true(is.na(survivalData@subject.data$ttr.cens[1]))
  
})

test_that("can_create_SurvivalData_with_a_subject_with_empty_string_for_endpoint_data",{
  data("sibylData")
  
  inputs <- survivalDataConstuctorTestSetUp()
  
  sibylData$ttr[1] <- ""
  sibylData$ttr.cens[1] <- ""
  
  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               covDef=inputs$cov,
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr")
  
  expect_equal(class(survivalData@subject.data$ttr), "numeric")
  expect_equal(class(survivalData@subject.data$ttr.cens), "logical")
  expect_true(is.na(survivalData@subject.data$ttr[1]))
  expect_true(is.na(survivalData@subject.data$ttr.cens[1]))
  
  
})


context("SurvivalData_invalid_arguments")

test_that("error_if_invalid_endPointUnit",{
  data("sibylData")
  
  inputs <- survivalDataConstuctorTestSetUp()
  expect_error(SurvivalData(data=sibylData,
               armDef=inputs$arm,
               subjectCol="ID",
               endPointNames="relapse",
               censorCol="ttr.cens",
               timeCol="ttr",
               endPointUnit = "madeUpUnits"))
})


test_that("error_if_subject_ID_not_unique",{
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()

  #set non-unique-ID column
  sibylData$nonUniqueID <- 1:nrow(sibylData)
  sibylData$nonUniqueID[1] <- 2

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="nonUniqueID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))


})

test_that("error_if_data_is_not_a_data_frame",{
  inputs <- survivalDataConstuctorTestSetUp()

  expect_error(SurvivalData(data=c(3,4,5),
                            armDef=inputs$arm,
                            subjectCol="nonUniqueID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))

})

test_that("error_if_subgroup_displayNames_not_unique",{
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()
  inputs$sub[[1]]@displayName <- "theName"
  inputs$sub[[2]]@displayName <- "theName"
  
  expect_error(SurvivalData(data = sibylData,
                             armDef = inputs$arm,
                             covDef = inputs$cov,
                             subgroupDef = inputs$sub,
                             subjectCol = "ID",
                             endPointNames = c("relapse", "newEndpoint"),
                             censorCol = c("ttr.cens", "cens.2"),
                             timeCol = c("ttr", "end.2")))
  
})

test_that("error_if_covariate_displayNames_not_unique",{
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()
  inputs$cov[[1]]@displayName <- "theName"
  inputs$cov[[2]]@displayName <- "theName"
  
  expect_error(SurvivalData(data = sibylData,
                            armDef = inputs$arm,
                            covDef = inputs$cov,
                            subgroupDef = inputs$sub,
                            subjectCol = "ID",
                            endPointNames = c("relapse", "newEndpoint"),
                            censorCol = c("ttr.cens", "cens.2"),
                            timeCol = c("ttr", "end.2")))
  
})

test_that("error_if_endPointNames_not_unique",{
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()
  expect_error(SurvivalData(data = sibylData,
                            armDef = inputs$arm,
                            covDef = inputs$cov,
                            subgroupDef = inputs$sub,
                            subjectCol = "ID",
                            endPointNames = c("MyName", "MyName"),
                            censorCol = c("ttr.cens", "cens.2"),
                            timeCol = c("ttr", "end.2")))
  
})

test_that("error_if_subgroup_also_covariate",{
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()
  
  inputs$cov[[3]] <-  ColumnDef(columnName = "sub.isMale",
                                type = "logical",
                                displayName = "the diplay names")
  
  
  expect_error(SurvivalData(data = sibylData,
                            armDef = inputs$arm,
                            covDef = inputs$cov,
                            subgroupDef = inputs$sub,
                            subjectCol = "ID",
                            endPointNames = c("ttr", "endpoint2"),
                            censorCol = c("ttr.cens", "cens.2"),
                            timeCol = c("ttr", "end.2")))
})


context("SurvivalData_columnDef_mismatches")

test_that("no_error_if_1_element_list_of_arm_columns", {
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()
  inputs$arm <- list(inputs$arm)

  SurvivalData(data=sibylData,
               armDef=inputs$arm,
               subjectCol="ID",
               endPointNames="relapse",
               censorCol="ttr.cens",
               timeCol="ttr")
})

test_that("error_if_list_of_more_than_1_arm_column", {
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()
  inputs$arm <- list(ColumnDef(columnName = "grp",
                               type = "categorical",
                               categories = factor(c("patchOnly", "combination"),
                                                   levels=c("patchOnly", "combination"))),
                     ColumnDef(columnName = "race",
                               type = "categorical",
                               categories = factor(c("black", "white", "hispanic", "other"),
                                                   levels=c("black", "white", "hispanic", "other"))))

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))
})

test_that("error_if_arm_column_not_in_data_set",{
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()
  inputs$arm@columnName <- "missingArm"


  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))
})

test_that("error_if_covariate_column_not_in_data_set",{
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()
  inputs$cov[[1]]@columnName <- "missingColName"

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            covDef = inputs$cov,
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))
})

test_that("error_if_endpoint_time_column_not_in_data_set",{
  data("sibylData")

  inputs <- survivalDataConstuctorTestSetUp()

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            subgroupDef=inputs$sub,
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="notInColumn"))
})


context("Invalid_data")

test_that("error_if_times_non-numeric_NA_or_negative", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]

  allInvalid <- list(c("1", "hello", "3"),
                     c(TRUE, FALSE, TRUE),
                     factor(c(2, 3, 4)),
                     c(-1, 1.1, 1.2),
                     c(-1, -2.7, 3.1))

  errorIfColumnInvalid(sibylData, inputs, "ttr", allInvalid)
})

test_that("error_if_censor_values_not_TRUE/FALSE_or_0/1", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]

  allInvalid <- list(c(0, 1, 2),
                     as.factor(c(1, 0, 1)),
                     c(-1, 0, 0),
                     c(1.01, 0, 0),
                     c(1, 1, 0.00001),
                     c(1, 1, -0.00001))

  errorIfColumnInvalid(sibylData, inputs, "ttr.cens", allInvalid)
})


test_that("error_if_timeCol_or_cenosrCol_misssing_but_other_is_not",{
  inputs <- survivalDataConstuctorTestSetUp()
  
  data("sibylData")
  sibylData <- sibylData[1:3, ]
  
  time <- list(c(NA, 5, 10),
               c("4","5",""),
               c(5,8,4),
               c(0, 5, 2),
               c(NA,3,2))
  
  cens <- list(c(TRUE,FALSE,TRUE),
               c(0,1,1),
               c("TRUE","FALSE",""),
               c(NA, 1, 1),
               c(NA, NA, "1"))
  
  
  errorIfColumnInvalid(sibylData, inputs, "ttr", time, 
                       colName2 = "ttr.cens", invalidValues2 = cens)
  
})

test_that("error_if_an_endpoint_has_no_time_or_censor_col",{
  inputs <- survivalDataConstuctorTestSetUp()
  
  data("sibylData")
  sibylData <- sibylData[1:3, ]
  
  time <- list(as.character(c(NA, NA, NA)),
               as.numeric(c(NA,NA,NA)))
  
  cens <- list(as.numeric(c(NA, NA, NA)),
               c(NA, NA, NA))
  
  
  errorIfColumnInvalid(sibylData, inputs, "ttr", time, 
                       colName2 = "ttr.cens", invalidValues2 = cens)
  
})


test_that("error_if_subgroup_values_not_TRUE/FALSE_or_0/1", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]

  allInvalid <- list(c(0, 1, 2),
                     c(-1, 0, 0),
                     c(0, NA, 0),
                     c(0, 1, NA),
                     c(1.01, 0, 0),
                     c(1, 1, 0.00001),
                     c(1, 1, -0.00001))

  errorIfColumnInvalid(sibylData, inputs, "sub.isMale", allInvalid)
})

test_that("error_if_fewer_than_2_levels_specified_for_arm", {
  inputs <- survivalDataConstuctorTestSetUp()
  inputs$arm@categories <- factor(c("patchOnly"), levels="patchOnly")

  data("sibylData")

  # Set all values to agree with defined levels, to avoid getting a different
  # error (values don't match levels)
  sibylData$grp <- "patchOnly"

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))

})

test_that("error_if_no_data_for_a_defined_arm", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")

  # Set all values of arm column to one value. The second arm therefore has no data.
  sibylData$grp <- "patchOnly"

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))
})

test_that("error_if_unexpected_value_in_arm_columns", {
  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")

  # Set an unexpected value in the arm column. Set the levels of the column to
  # include this value, so we don't get a warning from this.
  levels(sibylData$grp) <- c("patchOnly", "combination", "unexpectedValue")
  sibylData$grp[1] <- "unexpectedValue"

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))
})

test_that("error_if_specified_arm_not_in_raw_data", {

  inputs <- survivalDataConstuctorTestSetUp()
  inputs$arm <- ColumnDef(columnName = "NonExistentColumn",
                          type = "categorical",
                          categories = factor(c("patchOnly", "combination"),
                                              levels=c("patchOnly", "combination")))

  data("sibylData")
  expect_error(SurvivalData(data=testData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            endPointNames="relapse",
                            censorCol="ttr.cens",
                            timeCol="ttr"))
})

test_that("error_if_values_of_categorical_covariates_don't_match_levels", {
  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]

  allInvalid <- list(c("Black", "hispanic", "other"),
                     c("black", "HISPANIC", "other"),
                     c("black", "hispanic", "non-existent"),
                     c("white", 0, "other"),
                     c("white", 0, 1),
                     c(0, 1, 2),
                     c(1, 2, 3),
                     c(3, 4, 5),
                     c(-1, 0, 1))

  errorIfColumnInvalid(sibylData, inputs, "race", allInvalid)
})

test_that("error_if_NA_in_arm_column", {
  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]

  allInvalid <- list(c("patchOnly", "combination", NA),
                     c("patchOnly", NA, NA),
                     c(NA, NA, NA))

  errorIfColumnInvalid(sibylData, inputs, "grp", allInvalid)
})

test_that("error_if_NA_in_subgroup_column", {
  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]

  allInvalid <- list(c(TRUE, FALSE, NA))

  errorIfColumnInvalid(sibylData, inputs, "sub.isMale", allInvalid)
})

test_that("no_error_for_NA_as_value_of_categorical_covariate", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]
  sibylData$race <- c("black", "white", NA)

  SurvivalData(data=sibylData,
               armDef=inputs$arm,
               subjectCol="ID",
               covDef = inputs$cov,
               endPointNames="relapse",
               censorCol="ttr.cens",
               timeCol="ttr")

})

test_that("error_if_different_number_of_endpoint_names_and_time_columns", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")
  sibylData <- sibylData[1:3, ]

  # Pairs of inconsistent names and time columns
  allEndPoints <- list(endPoints = list(c(), c("relapse"), c("relapse"), c("relapse", "relapse")),
                       timeCols = list(c("ttr"), c(), c("ttr", "ttr"), c("ttr")))

  for (idx in seq_len(length(allEndPoints[["endPoints"]]))){

    thisNames <- allEndPoints[["endPoints"]][[idx]]
    thisTimes <- allEndPoints[["timeCols"]][[idx]]

    expect_error(SurvivalData(data=sibylData,
                              armDef=inputs$arm,
                              subjectCol="ID",
                              covDef = inputs$cov,
                              subgroupDef = inputs$sub,
                              endPointNames=thisNames,
                              censorCol="ttr.cens",
                              timeCol=thisTimes))
  }
})

test_that("error_if_different_number_of_endpoint_names_and_censor_columns", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")

  # Pairs of inconsistent names and time columns
  allEndPoints <- list(endPoints = list(c(), c("relapse"), c("relapse"), c("relapse", "relapse")),
                       censorCols = list(c("ttr.cens"), c(), c("ttr.cens", "ttr.cens"), c("ttr.cens")))

  for (idx in seq_len(length(allEndPoints[["endPoints"]]))){

    thisNames <- allEndPoints[["endPoints"]][[idx]]
    thisCensor <- allEndPoints[["censorCols"]][[idx]]

    expect_error(SurvivalData(data=sibylData,
                              armDef=inputs$arm,
                              subjectCol="ID",
                              covDef = inputs$cov,
                              subgroupDef = inputs$sub,
                              endPointNames=thisNames,
                              censorCol=thisCensor,
                              timeCol="ttr"))
  }
})

test_that("error_if_repeated_event_censor_columns", {

  inputs <- survivalDataConstuctorTestSetUp()

  data("sibylData")

  expect_error(SurvivalData(data=sibylData,
                            armDef=inputs$arm,
                            subjectCol="ID",
                            covDef = inputs$cov,
                            subgroupDef = inputs$sub,
                            endPointNames=c("relapse", "relapse2"),
                            censorCol=c("ttr.cens", "ttr.cens"),
                            timeCol=c("ttr", "ttr.cens")))

})


context("minOfMaxObserved")

test_that("error_if_object_not_SurvivalData",{
  expect_error(minOfMaxObserved(object="erer",endPointName = "relapse"))
})



test_that("error_if_invalid_endPoint_selected",{
  survivalData <- createSurvivalDataObject()
  expect_error(minOfMaxObserved(survivalData,endPointName = "NotAnEndpoint"))
})

test_that("minimum_is_calculated",{
  survivalData <- createSurvivalDataObject()
  result <- minOfMaxObserved(survivalData,endPointName = "relapse")
  
  arm1Result <- max(survivalData@subject.data$ttr[survivalData@subject.data$arm=="patchOnly"])
  arm2Result <- max(survivalData@subject.data$ttr[survivalData@subject.data$arm=="combination"])
  
  expect_equal(result, min(arm1Result,arm2Result))
})

