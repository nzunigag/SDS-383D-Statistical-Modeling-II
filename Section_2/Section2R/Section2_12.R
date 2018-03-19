
# Question 2.12 A Gaussian Linear Model
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

summary(data)
plot(data$age, data$distance)
plot(data$sex, data$distance)

library(MASS) # Had some confilcts with dplyr 
unloadNamespace("MASS")

#------------------------------------------------------
# FREQUENTIST LINEAR REGRESSION

# Least Squares
mlm <- 
  lm("distance ~ age + sex", data = data)$coefficients 
beta.lm.int <- 15.3856902
beta.lm.age <- 0.6601852
beta.lm.sex <- 2.3210227

# Ridge
mlmridge <- 
  lm.ridge("distance ~ intercept + age + sex + 0", 
           data = data, lambda = 1)
beta.ridge.int <- 13.9440406
beta.ridge.age <- 0.7655923
beta.ridge.sex <- 2.5793120

#------------------------------------------------------
# BAYESIAN LINEAR REGRESSION

# Bayesian Model Description:
# - Lambda = I and K = I 
# - Vague prior for the hyperparameters 
# - (Y|beta,omega) ~ N(X %*% beta, (omega * Lamda)^{-1})
# - (beta|omega) ~ N(mu, (omega * K)^{-1})
# - omega ~ Gamma(a, b)

# Hyperparameters we need E[omega]=1 and Var[omega]=inf (vague)
Lambda <- diag(rep(1, n)) # matrix nxn
K <- diag(rep(0.01, p)) # precision matriz pxp
mu <- matrix(0, p) # prior guess on mu
a <- 0.01 # prior sample size for the error variance
b <- 0.01 # prior sum of square errors for the error variance

# Updated parameters
# (beta|y,omega) ~ N(mu.new, (omega * K.new)^{-1})
# (omega|y) ~ Gamma(a.new, b.new) 

XtLambdaX <- t(X) %*% Lambda %*% X
beta.hat <- solve(XtLambdaX) %*% crossprod(X, Lambda) %*% y
y.hat <- X %*% beta.hat

K.new <- K + XtLambdaX
mu.new <- solve(K.new) %*% (XtLambdaX %*% beta.hat + K %*% mu)
a.new <- a + (n + 1) / 2
s <- t(mu) %*% K %*% mu + t(y.hat) %*% Lambda %*% y.hat
r <- t(mu) %*% (K + XtLambdaX) %*% mu
b.new <- as.numeric(b + 1 / 2 * (s - r))

# Sampling updated parameters
n.iter <- 5000
beta <- matrix(0, n.iter, p)
omega <- rep(0, n.iter)
for (i in 1:n.iter) {
  # omega
  omega.zero <- rgamma(1, a.new, b.new)
  omega[i] <- omega.zero
  # beta
  beta.zero <- mvrnorm(1, mu.new, solve(omega.zero * K.new))
  beta[i,1] <- beta.zero[1]
  beta[i,2] <- beta.zero[2]
  beta[i,3] <- beta.zero[3]
}

# Fit with burnin of 1000
beta.int <- mean(beta[1001:5000,1])
beta.age <- mean(beta[1001:5000,2])
beta.sex <- mean(beta[1001:5000,3])
omega.post <- mean(omega[1001:5000])

# Trace Plots
png(filename=paste('P2_12_Trace.png'),width=15,height=15,units="cm",res=200)
par(mfrow=c(2,2))
plot(omega, type = "l")
plot(beta[,1], type = "l")
plot(beta[,2], type = "l")
plot(beta[,3], type = "l")
dev.off()

#---------------------------------------------------
# RESULTS PLOTS

png(filename=paste('P2_12_Int.png'),width=15,height=10,units="cm",res=200)
plot(density(beta[1001:5000,1]), 
     xlab = 'intercept')
abline(v=beta.lm.int, col="blue", lty=1)
abline(v=beta.ridge.int, col="red", lty=2)
legend('topleft', legend=c("Least Squares", "Ridge"),
       col=c("blue", "red"), lty=1:2, cex=0.8)
dev.off()

png(filename=paste('P2_12_Age.png'),width=15,height=10,units="cm",res=200)
plot(density(beta[1001:5000,2]), 
     xlab = 'age')
abline(v=beta.lm.age, col="blue", lty=1)
abline(v=beta.ridge.age, col="red", lty=2)
legend('topleft', legend=c("Least Squares", "Ridge"),
       col=c("blue", "red"), lty=1:2, cex=0.8)
dev.off()

png(filename=paste('P2_12_Sex.png'),width=15,height=10,units="cm",res=200)
plot(density(beta[1001:5000,3]), 
     xlab = 'sex')
abline(v=beta.lm.sex, col="blue", lty=1)
abline(v=beta.ridge.sex, col="red", lty=2)
legend('topleft', legend=c("Least Squares", "Ridge"),
       col=c("blue", "red"), lty=1:2, cex=0.8)
dev.off()

