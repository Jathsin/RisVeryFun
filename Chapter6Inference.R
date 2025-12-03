#CHAPTER 6: Statistical Inference 

# Estimators

variance <- function(x, y, estimator) {
  mean <- sum(x*y)/sum(y)
  num <- sum((x-mean)**2)
  N <- sum(y)
  if (estimator) {
    N <- N-1
  }
  num/N
}



# When using class marks
breaks <- seq(0,100, by=10)
mids <- (breaks[-length(breaks)] + breaks[-1])/2 # xi but cool



# if sum(p) != 1 
normalise <- function(p) {
  if (abs(sum(p) - 1) >= 1e-6) {
    p <- p/sum(p)
  }
}


