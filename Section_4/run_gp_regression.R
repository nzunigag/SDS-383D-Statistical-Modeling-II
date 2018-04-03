rm(list = ls()) 
library(readr)
library(rstan)

#---------------------------------------------------------------------
# Data: Faithful

data(faithful)
gp_data<-read_rdump('faithful_data.R')

#---------------------------------------------------------------------
# Run Stan

filename = "gp_regression.stan"
stan_code <-readChar(filename, file.info(filename)$size)
fit <- stan(model_code=stan_code,data=gp_data,chains=3,iter=500,warmup=50,thin=10)
params<-extract(fit)

#---------------------------------------------------------------------
# Results

par(mfrow=c(1, 4))

alpha_breaks=10 * (0:50) / 50 - 5
hist(log(params$alpha), main="", xlab="log(alpha)", col="red", yaxt='n')

beta_breaks=10 * (0:50) / 50 - 5
hist(log(params$rho_1), main="", xlab="log(rho)", col="red", yaxt='n')

beta_breaks=10 * (0:50) / 50 - 5
hist(log(params$rho_2), main="", xlab="log(rho)", col="red", yaxt='n')

sigma_breaks=5 * (0:50) / 50
hist(params$sigma, main="", xlab="sigma", col="red", yaxt='n')

readline(prompt="Press [enter] to continue")

probs = c(0.1,0.5,0.9)
cred <- sapply(1:length(gp_data$x_predict), function(n) quantile(params$y_predict[,n], probs=probs))

plot(1, type="n", main=title, xlab="x", ylab="y", xlim=c(0, 120), ylim=c(-10, 10))
polygon(c(gp_data$x_predict, rev(gp_data$x_predict)), c(cred[1,], rev(cred[3,])), border = NA)
lines(gp_data$x_predict, cred[2,], lwd=2)

points(gp_data$x,gp_data$y)

