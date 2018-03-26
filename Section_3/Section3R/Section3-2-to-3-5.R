# Questions 3.2 to 3.5 Modeling non-Gaussian observations
rm(list = ls())
# Set seed for reproducible results
set.seed(2112)
library(tidyverse)

#---------------------------------------------------------------------
# Data: survival data from the Titanic

titanic <- read_csv('titanic.csv')
summary(titanic)

titanic <- na.omit(titanic) # Age has 447 NAs
y <- as.matrix(ifelse(titanic$Survived == "Yes",1,0))
X <- as.matrix(scale(titanic$Age))
plot(X, y)
n <- nrow(y)
p <- ncol(X)

#---------------------------------------------------------------------

sig <- function(a){
  # Sigmoid Function of a
  1 / (1 + exp(a))
}

n.log.lik <- function(y, X, beta){
  # Negative log likelihood of the posterior
  (beta * beta) / 2 + 
    sum((1 - y) * X * beta + log(1 + exp(-(X * beta)))) 
}

MAP <- optim(0, function(beta) - MAP(y, X, beta), method = "Brent", lower = -1, upper = 1)
MAP$par 
