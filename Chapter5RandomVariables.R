# ===================== DISCRETE VARIABLE ==============================

# TO ADD: corrections by continuity


# ====== Uniform ======

# Probability function
pUniform <- function(k){  # k := number of elements in sample space -> equally probable
  sol <- 1/k
  sol
}

# Estimators
uniMean <- function(xi, k){
  sol <- sum(xi*(1/k))
  sol
}
uniVar <- function(mean, k, xi){
  sol <- sum(((xi - mean)^2)*1/k)
  sol
}

# ====== Bernoulli (p+q = 1) ====== 

pBernoulli <- function(x, p) {
  if (x == 1) {
    p
  } else {
    1-p
  }
}

bernMean <- function(p){
  p
}
bernVar <- function(p){
  q <- 1-p
  sol <- p*q
  sol
}

# ====== Binomial ======

# P(X = x)  dbinom(x, n := number of experiments, p := probability success)
#Used to see the success probability, a fixed number of attempts, count on successes
pBinomial <- function(x,n,p) {
  dbinom(x,n, p) # x can be either a vector or a number
}

# P(X <= x) pbinom(x, n, p, lower.tail := X <= x ?, log.p := return probability in log form?) 
# lower.tail = FALSE <-> P(X > x)
dBinomial <- function(x,n,p,lowerTail,logP){
  pbinom(x,n,p,lowerTail,logP) 
}


# Alternative
dBinomial2 <- function(vector, n, p){
  sum(dbinom(vector, n, p))
}


# ====== Geometric :=  independent Bernoulli trials ======
#Number of attepmts until you get the first success

# P(X = x)
pGeometric <- function(k, p){ # k:= number of trials, p := probability of success
  q <- 1-p
  sol <- q^(k-1)*p
  sol
}

# P(X <= x)
dGeometric <- function(vector, p) {
  ps <- pGeometric(vector, p)
  sum(ps)
}

# Estimators
geoMean <- function(p){
  sol <- 1/p
  sol
}
geoVar <- function(p){
  q <- 1-p
  sol <- q/(p^2)
}

# ====== Poisson ======
#Ratio, for example, telephone exchange per minute (lambda = 400)
# P(X = x)
pPoisson <- function(num_successes, lambda) {
  dpois(num_successes, lambda)
}

dPoisson <- function(vector, lambda) {
  sum(dpois(vector, lambda))
}



# mean = var = lambda

  

# ===================== CONTINUOUS VARIABLE ==============================

#NOTE: no need to interpolate! Also, P(X = x) does not make sense with continuous variables
#                                    therefore we will use f(x) to refer to probability densities

# pnorm, ppois... are all left gives left side

# ====== Uniform ======
# P(X <= x)
pUnif <- function(x, a, b) {
  x*1/(b-a);
}

dUnif <- function(x, a, b) {
  punif(x, min=a, max=b)
}

uniMeanCont <- function(a,b) {
  (a+b)/2
} 

uniVarCont <- function(a,b) {
  (b-a)^2/12
}
 

# ====== Normal ======
# NOTE: do not correct by continuity in R!

# f(x)
fNormal <- function(x, mean, std) {
  p <- dnorm(x, mean, std) 
  cat("f(",x,") =", p, "\n")
  p
}
  
# P( X >= x)
dNormal <- function(a, mean, std) { 
  p <- 1 - pnorm(a, mean = mean, sd = std)
  cat("P[X >=",a, "] =", p, "\n")
  p
}
#dNormal(450, 400-0.5, 20)



# ====== Poisson ======

# lowerTail = TRUE <-> P(X <= x), else P(X > x)
# lowerTail = TRUE, logP = FALSE by default
dPoisson <- function(x, lambda, lowerTail, logP) {
  ppois(x, lambda, lowerTail, logP)
}


# ====== Chi-Squared ======
#Does it follow x distribution
#Also used for independency or interaction
# P(X <= q)
dChi2 <- function(p, df, lowerTAIL) { # p <- percentile and df <- degrees of freedom
  pchisq(p, df, lowerTAIL)
}

 
# What is the value q such that P(X <= q) = p
get_chi2 <- function(p, df) {
  qchisq(p, df)
}


# ====== Student-T ======
qt(c(lowerBound, upperBound), df) # lowerBound <- lower percentile
# upperBound <- upperPercentile, df <- degrees of freedom


# ====== Fisher-Snedecor F-Distribution ======
qf(p, df1, df2) # p <- percentile, df1 and df2 <- degrees of freedom

