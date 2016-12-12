#function to calcuate standard error
se <- function(x,na.rm){
  sd(x,na.rm = na.rm)/sqrt(length(x))
}
