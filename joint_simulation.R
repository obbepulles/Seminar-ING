library("readxl")
library(tseries)
library(Metrics)
library(sandwich)
library("aTSA")
library(forecast)
library(generics)
data <- read_xlsx("hypothetical_data_set.xlsx", skip = 1)
FTP <- as.data.frame(data[, 2 : 6])
DFTP <- as.data.frame(data[ , 1 ])
ER <- as.data.frame(data[ , 9 : 13])
DER <- as.data.frame(data[ , 8])
ftp_2y <- FTP[ , 1]
er <- ER[, 1]
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
  v <- par[3]
  sim <- rep(0, n)
  sim[1] <- lastdata * exp(-a / 12) + r * (1 - exp(-a / 12))
  for(i in 2 : n){
    sim[i] <- sim[(i - 1)] * exp(-a / 12) + r * (1 - exp(-a / 12))
    print(sim[i])
  }
  return(sim)
}
vasi_sim <- function(par, n, lastdata){
  a <- par[1]
  r <- par[2]
  v <- par[3]
  sim <- rep(0, n)
  for(i in 1 : n){
    m <- lastdata * exp(- a * t)
  }
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
jointpar_all <- ls_vasi_joint(er, ftp_2y)
x <- 120 ##x is number of months
x_months_ftp <- vasi_forecast(jointpar_all[4:6], x, ftp_2y[length(ftp_2y)])
x_months_er <- vasi_forecast(jointpar_all[1:3], x,  er[length(er)])