library("readxl")
library(tseries)
library(Metrics)
library(sandwich)
library("aTSA")
library(forecast)
library(generics)
library(mixR)
library(copula)
library(MASS)
library(yuima)
library(stats4)
library(expm)
library(Matrix)
data <- read_xlsx("hypothetical_data_set.xlsx", skip = 1)
FTP <- as.data.frame(data[, 2 : 6])
DFTP <- as.data.frame(data[ , 1 ])
ER <- as.data.frame(data[ , 9 : 13])
DER <- as.data.frame(data[ , 8])
er <- ER[, 1]
ftp_2y <- FTP[ , 1]
#--functions writing here-----
vasi_ml <- function(x){
  n <- length(x) - 1
  x_n_1 <- x[1 : n]
  x_n_2 <- x[2 : (n + 1)]
  a <- -12 * log((n * sum(x_n_1 * x_n_2) - sum(x_n_1) * sum(x_n_2)) / (n * sum(x_n_1 * x_n_1) - sum(x_n_1)^2))
  r <- 1 / (n * (1 - exp(-a / 12))) * (sum(x_n_2) - exp(a / 12) * sum(x_n_1))
  v <- 2*a / (n * (1 - exp(-a / 6))) * sum((x_n_2 - x_n_1 * exp(-a / 12) - r * (1 - exp(-a / 12)))^2)
  return(c(a, r, v))
}
vasi_forecast <-  function(par, n, lastdata){
  a <- par[1]
  r <- par[2]
  sim <- rep(0, n)
  sim[1] <- lastdata * exp(-a / 12) + r * (1 - exp(-a / 12))
  for(i in 2 : n){
    sim[i] <- sim[(i - 1)] * exp(-a / 12) + r * (1 - exp(-a / 12))
  }
  return(sim)
}
vasi_sim <- function(par, n, lastdata){
  a <- par[1]
  r <- par[2]
  v <- par[3]
  dt <- 1 / 12
  sim <- rep(0, n)
  for(i in 1 : n){
    m <- lastdata * exp(- a * i * dt) +  r * ( 1- exp(-a * i * dt))
    var <- v / (2 * a) * (1 - exp(-2 * a * i *dt))
    sim[i] <- rnorm(1, mean = m, sd = sqrt(var))
  }
  return(sim)
}
ls_vasi_joint <- function(x, y){
  x1 <- x[1 : (length(x) - 1)]
  x2 <- x[2 : length(x)]
  y1 <- y[1 : (length(y) - 1)]
  y2 <- y[2 : length(y)]
  modx <- lm(x2 ~x1)
  mody <- lm(y2 ~y1)
  ax <- (modx$coefficients[2] - 1) * (-12)
  rx <- modx$coefficients[1] * 12 / ax
  sx <- sqrt(var(modx$residuals)* 12)
  rho <- cor(modx$residuals, mody$residuals)
  ay <- (mody$coefficients[2] - 1) * -12
  ry <- mody$coefficients[1] * 12 / ay
  sy <- sqrt(var(mody$residuals)* 12)
  return(c(ax, rx, sx, ay, ry, sy, rho))
}
vasi_joint_sim <- function(par, n, lastdata1, lastdata2){
  ax <- par[1]
  rx <- par[2]
  vx <- par[3]^2
  ay <- par[4]
  ry <- par[5]
  vy <- par[6]^2
  rho <- par[7]
  dt <- 1 / 12
  sim1 <- rep(0, n)
  sim2 <- rep(0, n)
  for(i in 1 : n){
    m1 <- lastdata1 * exp(- ax* i * dt) +  rx * ( 1- exp(-ax * i * dt))
    m2 <- lastdata2 * exp(- ay* i * dt) +  ry * ( 1- exp(-ay * i * dt))
    var1 <- vx / (2 * ax) * (1 - exp(-2 * ax * i *dt))
    var2 <- vy / (2 * ay) * (1 - exp(-2 * ay * i *dt))
    ex <- rnorm(1, mean = 0, sd = 1)
    ey <- rho * ex + sqrt(1 - rho^2)*rnorm(1, mean = 0, sd = 1)
    sim1[i] <- m1 + sqrt(var1) * ex
    sim2[i] <- m2 + sqrt(var2) *ey
  }
  return(cbind(sim1, sim2))
}

vasi_fitted <- function(data, par){
  dt <- 1/12
  n<- length(data)
  a <- par[1]
  r <- par[2]
  v <- par[3]
  sim <- rep(0, n)
  sim[1] <- data[1]
  for(i in 2 : n){
    mu <- sim[1] * exp(-a * i * dt) + r * ( 1- exp(-a * i * dt))
    var <- v^2 / (2 * a) * (1 - exp(-2 * a * i *dt)) 
    print(var)
    sim[i] <- rnorm(1, mean = mu, sd = sqrt(var))
  }
  return(sim)
}

#--------------
test1 <- 30 
test2 <- 20
test_ftp1 <- ftp_2y[91:120]
test_ftp2 <- ftp_2y[101:120]
#--------------
  #random walk
rw <- rep(ftp_2y[90], test1)
rw2 <- rep(ftp_2y[100], test2)
rmse_rw <- rmse(test_ftp1, rw)
rmse_rw2 <- rmse(test_ftp2, rw2)

  #ARIMA(p, 1, q)
#(0, 0) --> rw
arima00_1 <- Arima(ftp_2y[1:90], order = c(0, 1, 0), include.mean = T)
arima00_2 <- Arima(ftp_2y[1:100], order = c(0, 1, 0), include.mean = T)
arima00_1$bic
arima00_2$bic
sim_arima00_1 <- forecast(arima00_1, h = 30)
sim_arima00_2 <- forecast(arima00_2, h = 20)
rmse_arima00_1 <- rmse(test_ftp1, sim_arima00_1$mean)
rmse_arima00_2 <- rmse(test_ftp2, sim_arima00_2$mean)
#(0, 1) 
arima01_1 <- Arima(ftp_2y[1:90], order = c(0, 1, 1), include.mean = T)
arima01_2 <- Arima(ftp_2y[1:100], order = c(0, 1, 1), include.mean = T)
arima01_1$bic
arima01_2$bic
sim_arima01_1 <- forecast(arima01_1, h = 30)
sim_arima01_2 <- forecast(arima01_2, h = 20)
rmse_arima01_1 <- rmse(test_ftp1, sim_arima01_1$mean)
rmse_arima01_2 <- rmse(test_ftp2, sim_arima01_2$mean)
#(1, 0) 
arima10_1 <- Arima(ftp_2y[1:90], order = c(1, 1, 0), include.mean = T)
arima10_2 <- Arima(ftp_2y[1:100], order = c(1, 1, 0), include.mean = T)
arima10_1$bic
arima10_2$bic
sim_arima10_1 <- forecast(arima10_1, h = 30)
sim_arima10_2 <- forecast(arima10_2, h = 20)
rmse_arima10_1 <- rmse(test_ftp1, sim_arima10_1$mean)
rmse_arima10_2 <- rmse(test_ftp2, sim_arima10_2$mean)
#(1, 1) 
arima11_1 <- Arima(ftp_2y[1:90], order = c(1, 1, 1), include.mean = T)
arima11_2 <- Arima(ftp_2y[1:100], order = c(1, 1, 1), include.mean = T)
arima11_1$bic
arima11_2$bic
sim_arima11_1 <- forecast(arima11_1, h = 30)
sim_arima11_2 <- forecast(arima11_2, h = 20)
rmse_arima11_1 <- rmse(test_ftp1, sim_arima11_1$mean)
rmse_arima11_2 <- rmse(test_ftp2, sim_arima11_2$mean)
#-----AR(1) model
ar1_1 <- Arima(ftp_2y[1:90], order = c(1, 0, 0), include.mean = T)
ar1_2 <- Arima(ftp_2y[1:100], order = c(1, 0, 0), include.mean = T)
ar1 <- Arima(ftp_2y, order = c(1, 0, 0), include.mean = T)
mu <- mean(dif_2y)
sd <- sd(dif_2y)
sim <- rep(0, 119)
sim[1] <- ftp_2y[1]+ rnorm(1, mean = mu, sd = sd)
for(i in 2 : 119 ){
  sim[i] <- sim[(i - 1)] * a + rnorm(1, mean = mu, sd = sd)
}
plot(sim, type = "l")
lines(ftp_2y[2:120], col = "green")
Box.test(ar1_1$residuals)
Box.test(ar1_2$residuals)
shapiro.test(ar1_1$residuals)
shapiro.test(ar1_2$residuals)
ar1_1$coef
ar1_2$coef
acf(ar1_1$residuals)
pacf(ar1_1$residuals)
acf(ar1_2$residuals)
pacf(ar1_2$residuals)
ggqqplot(ar1_1$residuals)
ggqqplot(ar1_2$residuals)
sim_ar1_1 <- forecast(ar1_1, h = 30)
sim_ar1_2 <- forecast(ar1_2, h = 20)
rmse_ar1_1 <- rmse(test_ftp1, sim_ar1_1$mean)
rmse_ar1_2 <- rmse(test_ftp2, sim_ar1_2$mean)
#------vasicek 
rmse_vasi_ls1 <- rmse(test_ftp1, vasi_forecast(ls_vasi_joint(er[1:90], ftp_2y[1:90])[4:6], 30, ftp_2y[90]))
rmse_vasi_ls2 <- rmse(test_ftp2, vasi_forecast(ls_vasi_joint(er[1:100], ftp_2y[1:100])[4:6], 20, ftp_2y[100]))
rmse_vasi_ml1 <- rmse(test_ftp1, vasi_forecast(vasi_ml(ftp_2y[1:90]), 30, ftp_2y[90]))
rmse_vasi_ml2 <- rmse(test_ftp2, vasi_forecast(vasi_ml(ftp_2y[1:100]), length(c(101:120)), ftp_2y[100]))
#-------Joint model 
sim_joint <- vasi_joint_sim(jointpar_all, x, er[length(er)], ftp_2y[length(ftp_2y)])
ftp_sim_joint <- sim_joint[, 2]
euri_sim_joint <- sim_joint[, 1]