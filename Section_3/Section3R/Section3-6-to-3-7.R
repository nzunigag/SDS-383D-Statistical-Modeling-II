# Questions 3.6 to 3.7 Modeling non-Gaussian observations
rm(list = ls())
# Set seed for reproducible results
set.seed(2112)
library(tidyverse)

#---------------------------------------------------------------------
# Data: numberof out of school suspensions

data <- read_csv('tea_discipline_oss.csv')
summary(data)

tea <- data %>% filter(ACTIONS >= 0) # Remove dropouts and NAs

summary(tea)
hist(tea$ACTIONS, breaks = 100)

y <- as.matrix(tea$ACTIONS)
n <- nrow(y)
X <- tea %>% mutate(intercept = rep(1, n)) %>% select(intercept, GRADE) 
p <- ncol(X)
plot(X$GRADE, y)

#---------------------------------------------------------------------


