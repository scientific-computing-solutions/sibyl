source("setupFunctions.R")

context("summarySurvivalData")

test_that("error_if_invalid_summary_type",{
  survivalData <- createSurvivalDataObject()
  expect_error(summary(survivalData,type="unknownType"))
})

test_that("subgroup_summary_displays_only_number_in_arm_if_no_subgroups",{
  data("sibylData")
  
  inputs <- survivalDataConstuctorTestSetUp()
  
  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr")  

  result <- extractSubgroupTable(survivalData)  
  
  expect_equal(nrow(result), 1)
  
  #number in each arm:
  patchNum <- sum(sibylData$grp=="patchOnly")
  comboNum <-sum(sibylData$grp=="combination") 
  
  expect_equal(as.numeric(result),c(patchNum,comboNum))
  
})

test_that("subgroup_summary_displays_one_row_per_subgroup_and one_column_per_arm",{
  survivalData <- createSurvivalDataObject()
  result <- extractSubgroupTable(survivalData) 
  
  #the 1+ is for the number in each arm row in the data frame
  expect_equal(nrow(result),1+length(survivalData@subgroupDef))
  expect_equal(ncol(result),length(survivalData@armDef@categories))
  
})

test_that("subgroup_summary_uses_subgroup_display_names",{
  survivalData <- createSurvivalDataObject()
  result <- extractSubgroupTable(survivalData)  
  
  #subgroup display names
  displayNames <- vapply(1:2, function(i){survivalData@subgroupDef[[i]]@displayName}, FUN.VALUE=character(1))
  
  expect_equal(rownames(result)[2:3], displayNames)
  
})

test_that("endpoint_summary_displays_results_for_all_and_each_subgroup_split_by_arm",{
  survivalData <- createSurvivalDataObject()
  result <- summary(survivalData,type="endPoints") 
  
  expect_equal(result$numcol,8) #2 (header cols) + 2 arms x (1[all] + 2 [subgroups])
  
})

test_that("endpoint_summary_has_only_'all'_column_if_no_subgroups",{
  data("sibylData")
  
  inputs <- survivalDataConstuctorTestSetUp()
  
  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr")
  result <- summary(survivalData,type="endPoints") 
  expect_equal(result$numcol,2+2) #2 header columns then 2 All columns
})

test_that("endpoint_summaryt_has_two_rows_per_endpoint",{
  data("sibylData")
  
  inputs <- survivalDataConstuctorTestSetUp()
  
  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               endPointNames=c("relapse","endpoint2"),
                               censorCol=c("ttr.cens","cens.2"),
                               timeCol=c("ttr","end.2"))
  result <-  summary(survivalData,type="endPoints") 
  
  expect_equal(result$numrow, 2*2)  #2 for each endpoint
  
})

test_that("extractEndPointOutput_uses_func_as_summary_function",{
  survivalData <- createSurvivalDataObject()
  
  
  result <- extractEndPointOutput(survivalData,
                  function(time,cens){
                    as.character(round(0.567845,3))
                  })
  
  ans <- t(data.frame(relapse=as.character(rep(0.568,6)),
                    newEndpoint=as.character(rep(0.568,6))))
  
  colnames(ans) <- rep(c("patchOnly","combination"),3)
  expect_equal(ans, result)
  
})

test_that("extractEndPointOutput_with_func_as_maturity_calculates_maturity",{
  survivalData <- createSurvivalDataObject()
  
  
  result <- extractEndPointOutput(survivalData,
                                  maturityFunc <- function(time, cens){
                                    as.character((length(cens)-sum(cens))) 
                                  })
  
  #isMale maturity:
  #Number in each arm who had events in subgroup male
  isMaleData <- survivalData@subject.data[survivalData@subject.data$sub.isMale,]
  maturityIsMale <- as.character(c(sum(!isMaleData$ttr.cens & isMaleData$arm == "patchOnly"),
                      sum(!isMaleData$ttr.cens & isMaleData$arm == "combination")))
  
  names(maturityIsMale) <- c("patchOnly","combination")
  
  expect_equal(result[1,3:4],maturityIsMale)
})


test_that("no_covariates_produces_NULL_covariate_summaries",{
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()
  survivalData <-  SurvivalData(data=sibylData,
                                armDef=inputs$arm,
                                subjectCol="ID",
                                subgroupDef=inputs$sub,
                                endPointNames="relapse",
                                censorCol="ttr.cens",
                                timeCol="ttr")
  
  result <- summary(survivalData,type="covariates")
  
  expect_true(is.null(result$numeric))
  expect_true(is.null(result$categorical))
  
})

test_that("logical_covariates_are_treated_as_factors_for_covariate_summaries",{
  
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()
  
  sibylData$covMale <- sibylData$sub.isMale
  
  covariateDef <- list(
    ColumnDef(columnName = "covMale",
              type = "logical",
              displayName = "Male"),
    ColumnDef(columnName = "race",
              type = "categorical",
              categories = factor(c("black", "hispanic", "other", "white"),
                                  levels=c("black", "white", "hispanic", "other"))))
  
  
  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               subgroupDef=inputs$sub,
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr",
                               covDef=covariateDef)
  result <- summary(survivalData,type="covariates")
  expect_equal(result$categorical$numrow,6) #4 for race and 2 for male
})

test_that("extractCovariateOutput_keeps_order_of_categories",{
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()
  
  sibylData$covMale <- sibylData$sub.isMale
  
  covariateDef <- list(
    ColumnDef(columnName = "covMale",
              type = "logical",
              displayName = "male"),
    ColumnDef(columnName = "race",
              displayName = "race",
              type = "categorical",
              categories = factor(c("black", "hispanic", "other", "white"),
                                  levels=c("white", "black", "hispanic", "other"))))
  
  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               subgroupDef=inputs$sub,
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr",
                               covDef=covariateDef)
  
  #dummy summary function outputting the string "1" for all results
  f <- function(covVals, cens){
    #coerce logical to factor for the analysis
    if(is.logical(covVals)){
      covVals <- factor(covVals, levels=c(TRUE,FALSE)) 
    }
    
    #named vector of results one per category
    vapply(levels(covVals),function(x){
      as.character(1)
    },character(1))
  }
  
  #have to change logical covariate to categorical 
  survivalData <- convertMissingFactorsToOwnLevel(survivalData)
  result <- extractCovariateOutput(survivalData, func=f, requiredTypes = c("logical","categorical"),endPoint=NULL)
  
  expect_equal(rownames(result),c("TRUE","FALSE","white","black","hispanic","other"))
  
})


test_that("subgroup_as_a_covariate_sets_percent_as_100_in_table",{
  data("sibylData")
  inputs <- survivalDataConstuctorTestSetUp()

  sibylData$covMale <- sibylData$sub.isMale
  covariateDef <- list(
    ColumnDef(columnName = "covMale",
              type = "logical",
              displayName = "Male"))

  survivalData <- SurvivalData(data=sibylData,
                               armDef=inputs$arm,
                               subjectCol="ID",
                               subgroupDef=inputs$sub,
                               endPointNames="relapse",
                               censorCol="ttr.cens",
                               timeCol="ttr",
                               covDef=covariateDef)
  
  #summary function to get the number and percentage
  f <- getTypeSpecificValues(FALSE, digits=2, endPoint=NULL)$summaryFunc
  
  #have to change logical covariate to categorical 
  survivalData <- convertMissingFactorsToOwnLevel(survivalData)
  
  result <- extractCovariateOutput(survivalData, func=f, requiredTypes = c("logical","categorical"), endPoint=NULL)

  expect_equal(result[1,3], "100 (42)") #male patchOnly = 100%
  expect_equal(result[1,4], "100 (42)") # male combination = 0
  expect_equal(result[2,3], "0 (0)") #not male patchOnly = 0
  expect_equal(result[2,4], "0 (0)") #not male combination = 0
})


test_that("covariate_maturity_summary_function_calculates_maturity",{
  fn <- categoricalMaturityTypeSpecific()$summaryFunc
  
  vals <- c("A","B","A","B","C","A","C","C","A")
  vals <- factor(vals,levels=c("A","B","C","D"))
  cens <- c(TRUE,FALSE,NA,FALSE,NA,FALSE,NA,NA,TRUE)
  
  ans <- c("1/3","2/2","-/0","-/0")
  names(ans) <- c("A","B","C","D")  
  
  expect_equal(fn(vals,cens), ans)
})
