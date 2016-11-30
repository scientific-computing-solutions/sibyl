# Test helper function - set up standard inputs to SurvivalData()
survivalDataConstuctorTestSetUp <- function(){
  armDef <- ColumnDef(columnName = "grp",
                      type = "categorical",
                      categories = factor(c("patchOnly", "combination"),
                                          levels=c("patchOnly", "combination")))
  
  covariateDef <- list(
    ColumnDef(columnName = "age",
              type = "numeric",
              unit= "years"),
    ColumnDef(columnName = "race",
              type = "categorical",
              categories = factor(c("black", "hispanic", "other", "white"),
                                  levels=c("black", "white", "hispanic", "other"))))
  
  # Define columns corresponding to subgroups
  subgroupDef <- list(
    ColumnDef(columnName = "sub.isMale",
              type = "logical",
              displayName = "Male"),
    ColumnDef(columnName = "sub.isHeavySmoker",
              type = "logical",
              displayName = "Heavy smoker"))
  
  
  list(arm=armDef,
       cov=covariateDef,
       sub=subgroupDef)
  
}

createSurvivalDataObject <- function(){
  
  inputs <- survivalDataConstuctorTestSetUp()
  
  SurvivalData(data = sibylData,
               armDef = inputs[["arm"]],
               covDef = inputs[["cov"]],
               subgroupDef = inputs[["sub"]],
               subjectCol = "ID",
               endPointNames = c("relapse", "newEndpoint"),
               censorCol = c("ttr.cens", "cens.2"),
               timeCol = c("ttr", "end.2")) 
}