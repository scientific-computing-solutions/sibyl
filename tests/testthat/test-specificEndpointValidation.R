context("SpecificEndpoinValidation")

#see Roxygen commens on specificEndpointRules function
#for more details


#endPoint is a vector of endPoint names which have timeCol name_TIME
#and censorCol name_CENS
simpleSurvivalDataCreation <- function(data, endPoint){
  data$subject <- 1:nrow(data)
  data$arm <- c("A", rep("B",nrow(data)-1))
  
  armDef <- ColumnDef(columnName="arm",
                      type="categorical",
                      categories=factor(c("A","B"),levels=c("A","B")))
  
  SurvivalData(data,armDef=armDef,
               subjectCol="subject",
               endPointNames=endPoint,
               censorCol=paste(endPoint,"CENS",sep="_"),
               timeCol=paste(endPoint,"TIME",sep="_"))
  
  
}

test_that("object_not_SurvivalData_throws_error",{
  expect_error(specificEndpointRules("hello"))
})

test_that("rulesToCheck_negative_throw_error",{
  data <- data.frame(PFS_CENS=c(1,0,1,1,0,1),
                     PFS_TIME=c(45,11,12,23,67,100))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS"))
  expect_error(specificEndpointRules(sD, -1))
})

test_that("rule1_ignored_if_no_PFS2_endpoint",{
  
  data <- data.frame(PFS_CENS=c(1,0,1,1,0,1),
                     PFS_TIME=c(45,11,12,23,67,100))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS"))
  expect_true(is.null(specificEndpointRules(sD, 1)))
})

test_that("rule1_ignored_if_no_PFS_endpoint",{
  
  data <- data.frame(PFS2_CENS=c(1,0,1,1,0,1),
                     PFS2_TIME=c(45,11,12,23,67,100))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS2"))
  expect_true(is.null(specificEndpointRules(sD, 1)))
})

test_that("NULL_returned_if_rule1_not_violated",{
  data <- data.frame(PFS_CENS=c(1,0,1,1,0,1),
                     PFS_TIME=c(45,11,12,23,67,100),
                     PFS2_CENS=c(NA,0,NA,NA,1,NA),
                     PFS2_TIME=c(NA,24,NA,NA,75,NA))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS2","PFS"))
  expect_true(is.null(specificEndpointRules(sD, 1)))
})

test_that("Warning_string_returned_if_rule1_violated",{
  data <- data.frame(PFS_CENS=c(1,0,1,1,0,1),
                     PFS_TIME=c(45,11,12,23,67,100),
                     PFS2_CENS=c(NA,0,NA,NA,1,NA),
                     PFS2_TIME=c(NA,10,NA,NA,75,NA))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS2","PFS"))
  expect_equal(specificEndpointRules(sD, 1),
    "WARNING: Subject 2 has PFS2 time < PFS time")
})

test_that("Warning_string_returned_if_rule1_violated_more_than_once",{
  data <- data.frame(PFS_CENS=c(1,0,1,1,0,1),
                     PFS_TIME=c(45,11,12,23,100,100),
                     PFS2_CENS=c(NA,0,NA,NA,1,NA),
                     PFS2_TIME=c(NA,10,NA,NA,75,NA))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS2","PFS"))
  expect_equal(specificEndpointRules(sD, 1),
               paste("WARNING: Subject 2 has PFS2 time < PFS time",
                     "WARNING: Subject 5 has PFS2 time < PFS time",sep="\n")
  )
})

test_that("rule_ignored_if_not_in_rulesToCheck",{
  data <- data.frame(PFS_CENS=c(1,0,1,1,0,1),
                     PFS_TIME=c(45,11,12,23,100,100),
                     PFS2_CENS=c(NA,0,NA,NA,1,NA),
                     PFS2_TIME=c(NA,10,NA,NA,75,NA))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS2","PFS"))
  expect_true(is.null(specificEndpointRules(sD, 2)))
})

test_that("rule_2_ignored_if_no_OS_endpoint",{
  data <- data.frame(PFS_CENS=c(1,0,1,1,0,1),
                     PFS_TIME=c(45,11,12,23,100,100),
                     PFS2_CENS=c(NA,0,NA,NA,1,NA),
                     PFS2_TIME=c(NA,10,NA,NA,75,NA))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("PFS2","PFS"))
  expect_true(is.null(specificEndpointRules(sD, 2)))
  
})

test_that("rule_2_ignored_if_only_OS_endpoint",{
  data <- data.frame(OS_CENS=c(1,0,1,1,0,1),
                     OS_TIME=c(45,11,12,23,100,100))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("OS"))
  expect_true(is.null(specificEndpointRules(sD, 2)))
})

test_that("rule_2_outputs_warnings_only_for_those_who_had_OS_event",{
  data <- data.frame(OS_CENS=c(1,0,0,0,0,NA),
                     OS_TIME=c(45,11,12,23,100,NA),
                     OTHER1_CENS=c(1,0,NA,1,1,0),
                     OTHER1_TIME=c(34,15,NA,25,56,10),
                     OTHER2_CENS=c(0,0,1,1,0,1),
                     OTHER2_TIME=c(34,15,12,15,56,15))
  
  sD <- simpleSurvivalDataCreation(data,endPoint=c("OS","OTHER1","OTHER2"))
  expect_equal(specificEndpointRules(sD, 2),
  paste("WARNING: Subject 2 has, for an endpoint other than OS, a time > OS event time",
        "WARNING: Subject 4 has, for an endpoint other than OS, a time > OS event time",
        sep="\n"))
})


context("getZeroTimes")
test_that("error_thrown_if_not_SurvivalData_object",{
 expect_error(getZeroTimes(c(4,5,6,7))) 
})


test_that("no_zero_times_returns_NULL",{
  data <- data.frame(OS_CENS=c(1,0,0,0,0,NA),
                     OS_TIME=c(45,11,12,23,100,NA),
                     OTHER_CENS=c(1,0,NA,1,1,0),
                     OTHER_TIME=c(23,56,NA,15,12,56))
  sD <- simpleSurvivalDataCreation(data,endPoint=c("OS","OTHER"))
  expect_true(is.null(getZeroTimes(sD)))
})

test_that("single_zero_time_is_output_as_warning_string",{
  data <- data.frame(OS_CENS=c(1,0,0,0,0,NA),
                     OS_TIME=c(0,11,12,23,100,NA),
                     OTHER_CENS=c(1,0,NA,1,1,0),
                     OTHER_TIME=c(23,1,NA,23,12,56))
  sD <- simpleSurvivalDataCreation(data,endPoint=c("OS","OTHER"))
  expect_equal(getZeroTimes(sD),
               "WARNING: Subject 1 has time=0 for endpoint OS and this subject will not be used when fitting parametric models for this endpoint")
})

test_that("zero_times_in_multiple_subjects_outputs_each_as_a_warning_string",{
  data <- data.frame(OS_CENS=c(1,0,0,0,0,NA),
                     OS_TIME=c(0,11,12,23,100,NA),
                     OTHER_CENS=c(1,0,NA,1,1,0),
                     OTHER_TIME=c(23,0,NA,0,12,56))
  sD <- simpleSurvivalDataCreation(data,endPoint=c("OS","OTHER"))
  expect_equal(getZeroTimes(sD),
               paste("WARNING: Subject 1 has time=0 for endpoint OS and this subject will not be used when fitting parametric models for this endpoint",
                     "WARNING: Subject 2 has time=0 for endpoint OTHER and had an event and this subject will not be used when fitting parametric models for this endpoint",
                     "WARNING: Subject 4 has time=0 for endpoint OTHER and this subject will not be used when fitting parametric models for this endpoint",sep="\n"))
})


test_that("zero_times_for_multiple_endpoints_for_same_subject_outputs_all_warning_strings",{
  data <- data.frame(OS_CENS=c(1,0,0,0,0,NA),
                     OS_TIME=c(0,11,12,23,100,NA),
                     OTHER_CENS=c(1,0,NA,1,1,0),
                     OTHER_TIME=c(0,56,NA,1,12,56))
  sD <- simpleSurvivalDataCreation(data,endPoint=c("OS","OTHER"))
  expect_equal(getZeroTimes(sD),
               paste("WARNING: Subject 1 has time=0 for endpoint OS and this subject will not be used when fitting parametric models for this endpoint",
                     "WARNING: Subject 1 has time=0 for endpoint OTHER and this subject will not be used when fitting parametric models for this endpoint", sep="\n")
  )
})
