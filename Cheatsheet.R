# Cheatsheet
n <- 10
x <- 1:n

# DISTRIBUTION + CUMULATIVE

#Normal

#Uniform

#Bernoulli

#Binomial
# Definition
P <- dbinom(x,n, 0.5) # x can be either a vector or a number

# Cumulative
# probability of obtaining 4 our fewer heads, these are equivalent
p_four_fewer_heads <- sum(dbinom(c(1,2,3,4), n, p_heads))
p_four_fewer_heads_2 <- pbinom(4,10,0.7,lower.tail=TRUE, log.p=FALSE) #prefer this one


