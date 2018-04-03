
weather <- read.csv("weather.csv")
N <- nrow(weather)
y <- weather$temperature
x_1 <- weather$lon
x_2 <- weather$lat 

N_predict <- 50 * 50
longrid <- seq(min(weather$lon), max(weather$lon), length.out = 50)
latgrid <- seq(min(weather$lat), max(weather$lat), length.out = 50)
grid <- expand.grid(longrid, latgrid)

x_1_predict <- grid$Var1
x_2_predict <- grid$Var2