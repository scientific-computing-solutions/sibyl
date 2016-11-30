#function to calcuate standard error
se <- function(x){
  sd(x)/sqrt(length(x))
}
