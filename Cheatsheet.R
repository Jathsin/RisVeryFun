# Cheatsheet
n <- 10
x <- 1:n


# DISCHRETE
# DISTRIBUTION + CUMULATIVE
#Normal

#Uniform

#Bernoulli

#Binomial
P <- dbinom(x,n, 0.5) # x can be either a vector or a number

# Cumulative
# probability of obtaining 4 our fewer heads, these are equivalent
p_four_fewer_heads <- sum(dbinom(c(1,2,3,4), n, p_heads))
p_four_fewer_heads_2 <- pbinom(4,10,0.7,lower.tail=TRUE, log.p=FALSE) #prefer this one



# Poisson
success <- 5
lambda <- 4
P <- dpois(success, 4)

#Cumulative
P <- ppois(1,4)

#INTERPOLATION! Use proportionality between triangles



# Generate random vectors following X distribution
data <- rpois(10000,100)

hist(data,
     main = "Histogram of Poisson(λ=4) Samples",
     xlab = "x",
     col = "skyblue",
     border = "white")


#CONTINUOUS VARIABLE

# cumulative of uniform
p <- punif(5, min=0, max=15)


# NORMAL


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



