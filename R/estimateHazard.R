#this function is an (unstable) estimate of the hazard given the survival function S(t)
#it should be replaced with something more stable before being used
#see http://stackoverflow.com/questions/11081069/calculate-the-derivative-of-a-data-function-in-r
estimateHazard <- function(S, t){
  warning("The estimate for h(t) is unlikely to be stable so proceed",
          "with extreme caution")

  H <- -log(S)
  spl <- smooth.spline(t, H)
  predict(spl, deriv=1, x=t)$y
}