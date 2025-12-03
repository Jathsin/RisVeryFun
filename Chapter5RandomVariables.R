# ===================== DISCRETE VARIABLE ==============================
# DISTRIBUTION + CUMULATIVE
# ====== Normal ======

# ====== Uniform ======
# Dist Func
uniformDistFunc <- function(k){
  sol <- 1/k
  sol
}
# Mean and Variance
uniMean <- function(xi, k){
  sol <- sum(xi*(1/k))
  sol
}
uniVar <- function(mean, k, xi){
  sol <- sum(((xi - mean)^2)*1/k)
  sol
}

# ====== Bernoulli ======
# p+q = 1
# Mean and Variance
bernMean <- function(p){
  p
}
bernVar <- function(p){
  q <- 1-p
  sol <- p*q
  sol
}

# ====== Binomial ======
#dbinom(x, n(totalSize), p(probability of each outcome))
P <- dbinom(x,n, 0.5) # x can be either a vector or a number

# Cumulative
# probability of obtaining 4 our fewer heads, these are equivalent
p_four_fewer_heads <- sum(dbinom(c(1,2,3,4), n, p_heads))
p_four_fewer_heads_2 <- pbinom(4,10,0.7,lower.tail=TRUE, log.p=FALSE) #prefer this one


# ====== Geometric ======
# probab Func
geoProbFunc <- function(k, p){
  q <- 1-p
  sol <- q^(k-1)*p
  sol
}
# Mean and Variance
geoMean <- function(p){
  sol <- 1/p
  sol
}
geoVar <- function(p){
  q <- 1-p
  sol <- q/(p^2)
}

# ====== Poisson ======
success <- 5
lambda <- 4 # mean
P <- dpois(success, 4)

#Cumulative
P <- ppois(1,4)
# Mean = lambda, Variance = lambda

# ====== Interpolation ====== Use proportionality between triangles



# Generate random vectors following X distribution
data <- rpois(10000,100)

hist(data,
     main = "Histogram of Poisson(λ=4) Samples",
     xlab = "x",
     col = "skyblue",
     border = "white")


# ===================== CONTINUOUS VARIABLE ==============================
# ====== Uniform ======
# cumulative of uniform
p <- punif(5, min=0, max=15)

# ====== Normal ======
# (do not correct by continuity in R!)

x <- 5
p <- pnorm(x, mean = 6, sd = 1.5) # Probability that X is less than x

p_normal_greater <- function(a, mean, std) {
  prob <- 1 - pnorm(a - 1, mean = mean, sd = std)
  # cat("P[X >= ", a, "] =", prob, "\n")
}
a <- p_normal_greater(450, 400-0.5, 20)


# Probability that X is greater equal than a
p_poison_greater<- function(a,lambda) {
  p <- 1 - ppois(a-1, lambda, lower.tail = TRUE)
}
p <- p_poison_greater(450, 400)


# ====== Chi-Squared ======
qchisq(p, df) # p <- percentile and df <- degrees of freedom


# ====== Student-T ======
qt(c(lowerBound, upperBound), df) # lowerBound <- lower percentile
# upperBound <- upperPercentile, df <- degrees of freedom


# ====== Fisher-Snedecor F-Distribution ======
qf(p, df1, df2) # p <- percentile, df1 and df2 <- degrees of freedom

