#CHAPTER 6: Statistical Inference 
variance <- function(x, y, estimator) {
  mean <- sum(x*y)/sum(y)
  num <- sum((x-mean)**2)
  N <- sum(y)
  if (estimator) {
    N <- N-1
  }
  num/N
}