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
er <- ER[, 2]
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
  #kappa
  ax <- (modx$coefficients[2] - 1) * (-12)
  #mu
  rx <- modx$coefficients[1] * 12 / ax
  #sigma
  sx <- sqrt(var(modx$residuals)* 12)
  #rho
  rho <- cor(modx$residuals, mody$residuals)
  ay <- (mody$coefficients[2] - 1) * -12
  ry <- mody$coefficients[1] * 12 / ay
  sy <- sqrt(var(mody$residuals)* 12)
  return(c(ax, rx, sx, ay, ry, sy, rho))
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
jointpar_all <- ls_vasi_joint(er, ftp_2y)
jointpar_all_ind <- jointpar_all
jointpar_all_ind[3] <- jointpar_all[3]^2
jointpar_all_ind[6] <- jointpar_all[6]^2
x <- 120 ##x is number of months
  #this is LS estimator for individual 
ftp_sim <- vasi_sim(jointpar_all_ind[4:6], x, ftp_2y[length(ftp_2y)])
euri_sim <- vasi_sim(jointpar_all_ind[1:3], x,  er[length(er)])
  #this is LS estimator for joint simulation
sim_joint <- vasi_joint_sim(jointpar_all, x, er[length(er)], ftp_2y[length(ftp_2y)])
ftp_sim_joint <- sim_joint[, 2]
euri_sim_joint <- sim_joint[, 1]

