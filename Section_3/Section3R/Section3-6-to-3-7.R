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

#---------------------------------------------------------------------
# Ex 3.6 
summary(tea)

png(filename=paste('P3_6_hist.png'),width=15,height=10,units="cm",res=200)
ggplot(data = tea, aes(tea$ACTIONS)) + geom_histogram(binwidth = 20)
dev.off()

#---------------------------------------------------------------------
# Ex 3.7

y <- as.matrix(tea$ACTIONS)
n <- nrow(y)
X <- tea %>% mutate(intercept = rep(1, n)) %>% select(intercept, GRADE) 
p <- ncol(X)
plot(X$GRADE, y)




