rm(list = ls()) 
library(readr)
library(rstan)

#---------------------------------------------------------------------
# Data: weather

gp_data_2 <- read_rdump('weather_data.R')

#---------------------------------------------------------------------
# Run Stan

filename = "gp_regression_temp.stan"
stan_code <- readChar(filename, file.info(filename)$size)
fit <- stan(model_code = stan_code, data = gp_data_2, chains=3,iter=500,warmup=50, thin=10)
params <- extract(fit)

#---------------------------------------------------------------------
# Results

par(mfrow=c(1, 2))
hist(log(params$alpha), main="", xlab="log(alpha)", col="red", yaxt='n')
hist(params$sigma, main="", xlab="sigma", col="red", yaxt='n')

par(mfrow=c(1, 2))
hist(log(params$rho_1), main="", xlab="log(rho_long)", col="red", yaxt='n')
hist(log(params$rho_2), main="", xlab="log(rho_lat)", col="red", yaxt='n')

quilt.plot(gp_data_2$x_1_predict, gp_data_2$x_2_predict, 
           colMeans(params$f_predict), xlab = 'Longitude', 
           ylab = 'Latitude', main = 'Heatmap of Temperature')
