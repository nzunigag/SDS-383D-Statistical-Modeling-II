# Question 2.14 A hierarchical Gaussian Linear Model
rm(list = ls())
library(dplyr)

#------------------------------------------------------
# DATA: DENTAL DISTANCE IN 27 CHILDREN 

dental <- read.csv('dental.csv')
n <- nrow(dental) 

data <- dental %>% 
  mutate(intercept = rep(1, n)) %>% 
  mutate(sex = as.numeric(Sex == "Male")) %>% 
  select(distance, intercept, age, sex)
y <- as.matrix(data$distance)
X <- as.matrix(data[,-1])
p <- ncol(X)

library(MASS)
#------------------------------------------------------
# BAYESIAN LINEAR REGRESSION

# GIBB SAMPLER
# Hyperparameters we need E[omega]=1 and Var[omega]=inf (vague)
K <- diag(rep(0.1, p)) # precision matriz pxp
mu <- matrix(0, p) # prior guess on mu
a <- 1 # prior sample size for the error variance
b <- 1 # prior sum of square errors for the error variance
tau <- 1

# Sampling updated parameters
n.iter <- 5000
beta <- matrix(NA, n.iter, p)
omega <- matrix(NA, n.iter)
lambda <- matrix(NA, n.iter, n)

beta[1,] <- rep(0, p)
lambda[1,] <- rep(1, n)
omega[1] <- 1

for (i in 2:n.iter) {
  Lambda <- diag(lambda[i-1,])
  XtLambdaX <- t(X) %*% Lambda %*% X
  beta.hat <- solve(XtLambdaX) %*% crossprod(X, Lambda) %*% y
  y.hat <- X %*% beta.hat
  K.new <- K + XtLambdaX
  mu.new <- solve(K.new) %*% (XtLambdaX %*% beta.hat + K %*% mu)
  
  # Betas
  beta[i,] <- mvrnorm(1, mu.new, solve(omega[i-1] * K.new))
  
  # omega  
  a.new <- a + (n + 1) / 2
  s <- t(mu) %*% K %*% mu + t(y.hat) %*% Lambda %*% y.hat
  r <- t(mu) %*% (K + XtLambdaX) %*% mu
  b.new <- as.numeric(b + 1 / 2 * (s - r))
  omega[i] <- rgamma(1, a.new, b.new)

  #lambda
  an <- tau + 1/2
  bn <- (1/2) * (omega[i] * (y - X %*% beta[i,])^2 + tau)
  lambda[i,] <- rgamma(1, an, bn)
}

# Fit with burnin of 1000
beta.int <- mean(beta[1001:5000,1])
beta.age <- mean(beta[1001:5000,2])
beta.sex <- mean(beta[1001:5000,3])

# Trace Plots
png(filename=paste('P2_14_Trace.png'),width=15,height=15,units="cm",res=200)
par(mfrow=c(2,2))
# plot(omega[1001:5000,1], type = "l")
plot(beta[1001:5000,1], type = "l")
plot(beta[1001:5000,2], type = "l")
plot(beta[1001:5000,3], type = "l")
dev.off()

res_1 <- abs(y--y_1) 
png(filename=paste('P2_15_Res1.png'),width=15,height=15,units="cm",res=200)
hist(res_1)
dev.off()
res_2 <- abs(y--y_2)
png(filename=paste('P2_15_Res2.png'),width=15,height=15,units="cm",res=200)
hist(res_2)
dev.off()


