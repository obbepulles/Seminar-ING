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
sim_sarima <- function(n, p, I, q, P, Q, s, pv, qv, Pv, Qv, c, model){
  res_0 <- model$residuals
  mu <- mean(res_0)
  sd <- sd(res_0)
  res <- rnorm(n, mean = mu, sd = sd)
  if (I == 0){  
    resa <- c(res_0, res)
    sim <- c(ftp_2y_ts, rep(0, n))
    l <- length(ftp_2y_ts)
    for(i in (l + 1) : (l + n)){
      sP <- Pv * sim[(i - s) : (i - s - P + 1)]
      sQ <- Qv * resa[(i - s) : (i - s - Q + 1)] 
      sP[is.na(sP)] <- 0
      sQ[is.na(sQ)] <- 0
      sim[i] <- sum(pv * sim[(i - 1) : (i - p)]) + sum(qv * sim[(i - 1) : (i - q)]) + c + sum(sP) + sum(sQ)
    }
    return(sim[(l + 1) : (l + n)])
  }
  
  else{
    dif <- diff(ftp_2y, lag = I)
    l <- length(dif)
    sim <- rep(0, n)
    resa <- c(res_0, res)
    sp <- pv * dif[(l) : (l - p + 1)]
    sq <- qv * res_0[(l) : (l - q + 1)]
    sP <- Pv * dif[(l - s + 1) : (l - s - P + 2)]
    sQ <- Qv * resa[(l - s + 1) : (l - s - Q + 2)]
    sp[is.na(sp)] <- 0
    sq[is.na(sq)] <- 0
    sP[is.na(sP)] <- 0
    sQ[is.na(sQ)] <- 0
    delta <- sum(sp) + sum(sq) + c + sum(sP) + sum(sQ)
    sim [1] <- ftp_2y[length(ftp_2y)] + delta
    dif <- c(dif, delta)
    for(i in 2 : n){
      sp <- pv * dif[(l + i) : (l + i- p + 1)]
      sq <- qv * res_0[(l + i) : (l + i- q + 1)]
      sP <- Pv * dif[(l + i - s + 1) : (l + i - s - P + 2)]
      sQ <- Qv * resa[(l + i - s + 1) : (l + i - s - Q + 2)]
      sp[is.na(sp)] <- 0
      sq[is.na(sq)] <- 0
      sP[is.na(sP)] <- 0
      sQ[is.na(sQ)] <- 0
      delta <- sum(sp) + sum(sq) + c + sum(sP) + sum(sQ)
      sim [i] <- sim[(i - 1)] + delta
      dif <- c(dif, delta)
    }
    return(sim)
  }
}
#--------------
#overall data visualization <- still need to fix the x-axis to time 
plot(DFTP, FTP[,1], type = "l")
plot(FTP[, 1], type = "l")
plot(FTP[, 1], ylim = c(-0.007, 0.012), type = "l")
for(i in 2 : 5){
  lines(FTP[, i])
}
hist(ftp_2y)
plot(density(ftp_2y))
#sp_FTP shows that we only need to model one maturity
sp_FTP<- cbind(FTP[ , 2]-FTP[ , 1], 0.5 * (FTP[ , 3] - FTP[ , 2]), 0.5 * (FTP[ , 4] - FTP[ , 3]), (FTP[ , 5] - FTP[ , 4]) / 3)
#visualization of the density and show the two peaks
d <- density(FTP[ , 1])
plot(d)
abline(v =modes(d)$x)
#modes(d)$x is the peak value of FTP_2y maturity
cor(FTP[ , 1], ER[ , 1:3])
print(modes(d)$x)
plot(DFTP, FTP[, 1], type = "l")
abline(h = modes(d)$x , col = "blue", lwd  = 2)
abline(h = mean(FTP[ , 1]), col = "pink", lwd = 2)
sd(FTP[ , 1])
plot( abs(diff(ftp_2y)), type = "l")
plot(diff(ftp_2y), type = "l")
library(ggpubr)
ggqqplot(diff(ftp_2y))
dif_2y <- diff(ftp_2y)
shapiro.test(diff(ftp_2y))
shapiro.test(ftp_2y)## ftp is not normal 
kpss.test(ftp_2y)#non-stationary
adf.test(ftp_2y)#non-stationary
kpss.test(dif_2y)#
adf.test(dif_2y)#this is contradict --> heteroskedasticity
library(vrtest)
vrtest::Auto.Q(ftp_2y)
plot(seq(-0.002, 0.002, by = 0.0001), dnorm(seq(from = -0.002, to =  0.002, by = 0.0001), mean = mean(dif_2y), sd = sd(dif_2y)), type = "l")
lines(density(dif_2y), col = "red")
mean(dif_2y)
cancel_ftp <- rep(0, (length(ftp_2y) - 5))
for(i in 1 : (length(ftp_2y) - 5)){
  cancel_ftp[i] <- mean(ftp_2y[ (1+i) : 120 ] - ftp_2y[1 : (120 - i)])
}
cancel_ftp
plot(c(1 : (length(ftp_2y) - 5)), cancel_ftp, type = "l", ylim = c(-0.02, 0.02))
lines(c(1 : (length(ftp_2y) -5)) * mean(dif_2y), col = "blue")
lines(c(1 : (length(ftp_2y) - 5)) * mean(dif_2y) + 2 *sd(dif_2y) * length(ftp_2y) / 12, col = "red")
rmse(cancel_ftp, c(1 : (length(ftp_2y) -5)) * mean(dif_2y) )
lm(cancel_ftp ~ c(1 : (length(ftp_2y) - 5)))
mean(dif_2y)
cancel_ftp
#------ with consideration of interest rate leave it later can be ignored now 
FTP_i <- cbind((FTP[, 1] - 0.001), FTP[ , 2:3], FTP[ , 5])
ER_i <- ER[ , 2:5]
plot(ftp_2y , type = "l" , ylim = c(-0.01 , 0.02))
lines(ER[ , 1], col = "blue")
plot(diff(ftp_2y), type = "l" , ylim = c(-0.01 , 0.01))
lines(diff(ER[ , 1]) , col = "blue")
#------ use for autocorrelation it need to be garch later on, it shows high autocorrelation
lines(diff(ER[ , 1]), col = "red")
kk <- lm(ftp_2y~ ER[ , 1])
adf.test(diff(ftp_2y))
print(acf(ftp_2y^2))
#------- trying to do normal draws from dif_ftp to predict
plot(density(diff(ftp_2y)))
library(mixR)
plot(pacf(dif_2y))
plot(acf(dif_2y))
plot(dif_2y, type = "l")
abline(h = mean(dif_2y))
x <- lm(dif_2y~diff(ER[,1]))
plot(acf(dif_2y^2))
Box.test(dif_2y, lag = 12, type = "Ljung-Box")
Box.test(dif_2y^2, lag = 12, type = "Ljung-Box")
##Ljung-box shows that it should be iid as it is not significant
train_dif <-  dif_2y[1 : 90]
train <- length(train_dif)
test_dif <- dif_2y[91 : length(dif_2y)]
train_ftp <- ftp_2y[1 : (train + 1)]
test_ftp <- ftp_2y[(train+2) : length(ftp_2y)]
test <- length(test_ftp)
mix_norm1_<- mixfit(train_dif, ncomp = 1, family = "normal")
Box.test(train_dif, lag = 12, type = "Ljung-Box")
n_sim <- 10000
##Ljung-box shows that in the training data it is still iid
#----- normal draws from the distribution of training and back-testing 
norm_predict <- matrix(0 , n_sim , test)
for (i in 1 : n_sim){
  norm_predict[i, 1] <-  rnorm(1, mix_norm1_$mu, sd = mix_norm1_$sd )
  for(j in 2 : test){
    norm_predict[i , j] <- rnorm(1, mix_norm1_$mu, sd = mix_norm1_$sd )
  }
}
norm_predict <- colMeans(norm_predict)
predict_norm_ftp <- rep (0, test)
predict_norm_ftp[1] <- train_ftp[length(train_ftp)] + norm_predict[1]
for(i in 2 : test){
  predict_norm_ftp[i] <- predict_norm_ftp[(i -1)] + norm_predict[i]
}
rmse_normal_draw_ftp <- rmse(test_ftp, predict_norm_ftp)
rmse_normal_dif <- rmse(test_dif, norm_predict)
plot(norm_predict)
points(test_dif, col = "blue")
mean(train_dif)
##-------random walk 
rw <- rep(ftp_2y[train + 1], test)
plot(test_ftp, type = "l")
AR1 <- lm(train_ftp[2 : length(train_ftp)] ~ train_ftp[1 : (length(train_ftp) - 1)])
AR_mean_res <- mean(AR1$residuals)
AR_sd_res <- sd(AR1$residuals)
plot(AR1$residuals)
summary(AR1)
ar1 <- rep(0, test)
ar1[1] <- train_ftp[length(train_ftp)] * AR1$coefficients[2] + AR1$coefficients[1]
for(i in 2: test){
  ar1[i] <- ar1[(i-1)] * AR1$coefficients[2] + AR1$coefficients[1]
}
plot(test_ftp, type = "l")
plot(ar1)
lines(ar1)
summary(AR1)
rmse_ar1 <- rmse(test_ftp, ar1)
pacf(dif_2y, lag = 36)
acf(dif_2y, lag = 36)
lines( rw, col = "blue")
rmse_rw <- rmse(test_ftp, rw)
rmse_normal <- rmse(test_dif, norm_predict)
plot(norm_predict, type = "l")
plot(test_dif, type = "l") ## this is very bad
pacf(dif_2y)
acf(dif_2y)
#----ARIMA(1,1,1)
library(forecast)
arima111 <- Arima(train_ftp_ts, order = c(1, 1, 1), include.mean = T, method = "CSS")
sim_arima111 <- forecast(arima111, h = 29)
plot(sim_arima111)
lines(ftp_2y_ts, col = "blue")
plot(sim_ar11$residuals)
rmse_arima111 <- rmse(test_ftp, sim_arima111$mean)
acf(arima111$residuals)
#AR1 package
ar1p <- Arima(train_ftp_ts, order = c( 1, 0, 0), method = "CSS",  include.constant =T )
ar1p$coef
sim_ar1p <- forecast(ar1p, h = 29)
plot(sim_ar1p)
rmse_ar1p <- rmse(test_ftp, sim_ar1p$mean)
#--- ARIMA( 0 ,1 , 1)
arima011 <- Arima(train_ftp_ts, order = c(0, 1, 1), include.mean = T, method = "CSS")
sim_arima011 <- forecast(arima011, h = 29)
rmse_arima011 <- rmse(test_ftp, sim_arima011$mean)
plot(sim_arima011)
lines(ftp_2y_ts, col = "black")
lines(arima011$fitted, col = "cornflowerblue", lty = 2, lwd = 1.2)
acf(arima011$residuals, lag.max = 48)
pacf(arima011$residuals, lag.max = 48)
#--- ARIMA( 12 ,1 , 0)
arima1210 <- Arima(train_ftp_ts, order = c(12, 1, 0), include.mean = T, method = "CSS")
sim_arima1210 <- forecast(arima1210, h = 29)
rmse_arima1210 <- rmse(test_ftp, sim_arima1210$mean)
plot(sim_arima1210)
lines(ftp_2y_ts, col = "black")
lines(arima1210$fitted, col = "cornflowerblue", lty = 2, lwd = 1.2)
legend("bottomleft", legend = c("fitted value","actual value", "predicted value"), col = c("cornflowerblue","black", "cyan2"), lty = c(2 , 1, 1), lwd = c(1, 1, 2))
#--- ARIMA( 0 ,1 , 12)
arima0112 <- arima(train_ftp_ts, order = c(0, 4, 0),seasonal = c(0,0,1,4))
sim_arima0112 <- forecast(arima0112, h = 29)
rmse_arima0112 <- rmse(test_ftp, sim_arima0112$mean)
plot(sim_arima1210)
lines(ftp_2y_ts, col = "black")
lines(arima1210$fitted, col = "cornflowerblue", lty = 2, lwd = 1.2)
legend("bottomleft", legend = c("fitted value","actual value", "predicted value"), col = c("cornflowerblue","black", "cyan2"), lty = c(2 , 1, 1), lwd = c(1, 1, 2))
#--- SARIMA()
sarima01012 <- arima(train_ftp_ts, order = c(1, 0, 0), seasonal = c(0,0,1,4), method = "ML")
sim_sarima0112 <- forecast(sarima01012, h = 29)
rmse_sarima0112 <- rmse(test_ftp, sim_sarima0112$mean)
plot(sim_sarima0112)
lines(ftp_2y_ts, col = "black")
acf(diff(ftp_2y, lag = 4), lag = 24)
pacf(diff(ftp_2y, lag = 4), lag = 24)
acf(sarima01012$residuals)
plot(sarima01012$residuals)
library("olsrr")
#ols_test_breusch_pagan(vasi) ## test result shows that the residuals of lm model has heteroskedasticity --> vasicek hand wavy ish (?)
#------ ARCH test on diff _ftp to check whether GARCH can be applied
library("FinTS")
ArchTest(dif_2y) ## you cannot model it via Garch
ArchTest(ftp_2y)
#----vasicek paper using MLE method by Axel Gerebenk
va <- vasi_forecast(vasi_ml(train_ftp), test, train_ftp[length(train_ftp)])
plot(ftp_2y[(train + 2): length(ftp_2y)], type = "l")
plot(expected_forecast)
rmse_vasi_paper <- rmse(test_ftp, va)
#----GARCH directly to FTP
# Assuming 'returns' is your financial return time series
library(rugarch)
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)))
fit <- ugarchfit(spec, data = train_ftp_ts)
sim_GARCH <- ugarchforecast(fit, n.ahead = 29)
#fitted(sim_GARCH)[1 : 29]
rmse_GARCH <- rmse(test_ftp, fitted(sim_GARCH))
#plot(sim_GARCH)
#lines(ftp_2y_ts)
#---- labor Hamilton for FTP -- EM algorithm 
  #initial 
set.seed(123)
T <- length(ftp_2y)
p_ij_t <- matrix(0, nrow= 4, ncol = T)
for(i in 1 : T){
  arr <- runif(4)
  p_ij_t[ , i] <- arr/sum(arr)
}
count <- 0
loop <- 1

while(count < 1000){    
  #M-step 
  p11 <- (sum(p_ij_t[1 , ]) - p_ij_t[1 , 1]) / (sum(p_ij_t[1 , ]) - p_ij_t[1, T] + sum(p_ij_t[2, ]) - p_ij_t[2, T])
  p22 <- (sum(p_ij_t[4 , ]) - p_ij_t[4 , 1]) / (sum(p_ij_t[4 , ]) - p_ij_t[4, T] + sum(p_ij_t[3, ]) - p_ij_t[3, T])
  mu_1 <- sum((p_ij_t[1 , ] + p_ij_t[2, ]) * ftp_2y) / (sum(p_ij_t[1 , ]) + sum(p_ij_t[2 , ]))
  mu_2 <- sum((p_ij_t[4 , ] + p_ij_t[3, ]) * ftp_2y) / (sum(p_ij_t[4 , ]) + sum(p_ij_t[3 , ]))
  sig_1 <- sum((p_ij_t[1 , ] + p_ij_t[2, ]) * (ftp_2y - mu_1)^2) /  (sum(p_ij_t[1 , ]) + sum(p_ij_t[2 , ]))
  sig_2 <- sum((p_ij_t[4 , ] + p_ij_t[3, ]) * (ftp_2y - mu_2)^2) /  (sum(p_ij_t[4 , ]) + sum(p_ij_t[3 , ]))
  print(mu_1)
  #E-step
  fs1 <- dnorm(ftp_2y, mean = mu_1, sd = sqrt(sig_1))
  fs2 <- dnorm(ftp_2y, mean = mu_2, sd = sqrt(sig_2))
  fs <- rbind(fs1, fs1, fs2, fs2)
  xi_tt <- matrix(0, nrow = 4, ncol = T)
  xi_pred <- matrix(0, nrow = 4, ncol = (T - 1)) # 2|1 to T|(T-1)
  xi_tt[ , 1] <- c(0.25, 0.25, 0.25, 0.25)
  P <- cbind(c(p11, 0, (1 - p11), 0), c(p11, 0, (1 - p11), 0), c(0, (1 - p22), 0, p22), c(0, (1 - p22), 0, p22))
  xi_tT <- matrix(0, nrow = 4, ncol = T)
  #from t = 2
  for(i in 1 : (T - 1)){
    #prediction
    xi_pred[ , i] <- P %*% xi_tt[ , i]
    #update 
    xi_tt[ ,(i + 1)] <- fs[ , (i + 1)] * xi_pred[ ,i] / (sum(fs[ , (i + 1)] * xi_pred[ ,i]))
  }
  #backward 
  xi_tT[ , T] <- xi_tt[ , T]
  for(i in 1 : (T - 1)){
    xi_tT[ , (T - i)] <- xi_tt[ , (T - i)] * t(P) %*% (xi_tT[ , (T - i + 1)] / xi_pred[ , (T - i)])
  }
  p_ij_t <- xi_tT
  count <- count + 1
}
t_test <- (mu_1 - mu_2) / sqrt(sig_1 / T + sig_2 / T)
t_test
p_value <- 2 * (1 - pt(abs(t_test),  df = ( 2*T  - 2))) # this shows there is no effect 
qt(0.95, df = 2*T -1)
p_value
#--- function of dFTP 
dFTP_func <- function(tau, data, T){
  #T for maturity , tau is cancellation time
  
}
summary(lm(dif_2y~ diff(ER_i[ , 1])))
#------ simulation for different model: # only generate n value tho
#1. SARIMA(1, 0, 0)(0, 0, 1)[4]; 2. GARCH(1, 1); 3. AR(1); 4. ARIMA(0, 1, 1) 5. RW; 6. Vasicek
ftp_forecast <- function(n, type, x){
  if(type == 1){
    sample_sarima <- Arima(x, order = c(1, 0, 0), seasonal = c(0,0,1,4), method = "ML")
    sim_1 <- forecast(sample_sarima, h = n)
    return(sim_1$mean)
  }
  else if(type == 2){
    spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)))
    fit1 <- ugarchfit(spec1, data = x)
    sim_2 <- ugarchforecast(fit1, n.ahead = n)
    return(fitted(sim_2))
  }
  else if(type == 3){
    sample_sarima <- Arima(x, order = c(1, 0, 0), method = "CSS")
    sim_1 <- forecast(sample_sarima, h = n)
    return(sim_1$mean)
  }
  else if(type == 4){
    sample_sarima <- Arima(x, order = c(0, 1, 1), method = "ML")
    sim_1 <- forecast(sample_sarima, h = n)
    return(sim_1$mean)
  }
  else if(type == 5){
    #y <- diff(x, lag  = 1)
    #res <- rnorm(n, mean = mean(y), sd = sd(y))
    sim <- rep(x[length(x)], n)
    return(sim)
  }
  else if (type == 6){
    sim <- vasi_forecast(vasi_ml(x), n, x[length(x)])
    return(sim)
  }
}
m <- cbind(ftp_forecast(29, 1, train_ftp), ftp_forecast(29, 2, train_ftp), ftp_forecast(29, 3, train_ftp), ftp_forecast(29, 4, train_ftp), ftp_forecast(29, 5, train_ftp), ftp_forecast(29, 6, train_ftp))
plot(ts(test_ftp), ylim = c(-0.0025, 0.0005), col = "red3", type = "l", lwd = 3, , ylab = "FTP in test set", xlab= "time index")
lines(ts(m[, 1]), lty = 1, col = "blue", lwd = 2)
lines(ts(m[, 2]), lty = 2, col = "green4", lwd = 2)
lines(ts(m[, 3]), lty = 3, col = "brown4", lwd = 2)
lines(ts(m[, 4]), lty = 4, col = "orange2", lwd = 2)
lines(ts(m[, 5]), lty = 5, col = "cornflowerblue", lwd = 2)
lines(ts(m[, 6]), lty = 6, col = "pink", lwd = 2)
legend("topright", legend = c("Actual", "SARIMA", "GARCH (1, 1)", "AR1", "ARIMA(0, 1, 1)", "RW", "Vasicek"), col = c("red3","blue", "green4", "brown4", "orange2", "cornflowerblue", "pink"), lty = c(1 , 1, 2, 3, 4, 5, 6), lwd = c(3, 2, 2, 2, 2, 2, 2))
#--- empirical cdf of ftp
par(mfrow = c(1, 1))
library(copula)
library(MASS)
emp_ftp <- rank(ftp_2y)/ (length(ftp_2y) + 1)
emp_eur <- rank(ER[, 1])/ (length(ER[, 1]) + 1)
plot(emp_ftp, emp_eur, pch = 16, col = "blue") ## dependence on right tail -> Gumble 
emp_UU <- cbind(emp_ftp, emp_eur)
gumble_obj <- gumbelCopula(dim = 2)
fit_gumble <- fitCopula(gumble_obj, emp_UU, method = "ml")
summary(fit_gumble)
theta <- 1.715
upper_dep <- 2 - 2^(1/ theta)
upper_dep
frank_obj <- frankCopula(dim = 2)
fit_frank <- fitCopula(frank_obj, emp_UU, method = "ml")
clayton_obj <- claytonCopula(dim = 2)
fit_clayton <- fitCopula(clayton_obj, emp_UU, method = "ml")
AIC(fit_gumble)
AIC(fit_frank)
AIC(fit_clayton)
# -------jointly vasicek 
library(mvtnorm)
pdf_mvn <- function(v, mu, cov){
  x <- exp(-0.5 * t(v - mu) %*% solve(cov) %*%(v - mu)) / sqrt((2 * pi)^2 * det(cov))
  return(x[1, 1])
}
log_likelihood <- function(par, data) {
  a_R <- par[1]
  r_R <- par[2] #long tern mean 
  sig_R <- par[3]
  a_ftp <- par[4]
  r_ftp <- par[5] #long tern mean
  sig_ftp <- par[6]
  rho_j <- par[7]
  dt <- 1 / 12 
  obs1_2 <- data[, 1][2 : length(data[ , 1])] 
  obs1_1<- data[ , 1][1 : (length(data[ , 1]) - 1)] #given data in conditional
  obs2_2 <- data[, 2][2 : length(data[ , 2])] 
  obs2_1 <- data[ , 2][1 : (length(data[ , 2]) - 1)] #given data in conditional 
  obs <- cbind(obs1_2, obs2_2)
  ll <- 0
  for(i in 1 : (length(obs2_1) + 1)){
    mu <- c(obs1_1[i] + a_R *(r_R - obs1_1[i]) * dt , obs2_1[i] * + a_ftp *(r_ftp - obs2_1[i]) * dt)
    cov_matrix <- cbind(c(sig_R^2 * dt, (sig_R * sig_ftp * rho_j * dt)), c((sig_R * sig_ftp * rho_j * dt), (sig_ftp^2 * (1 - rho_j^2) * dt)))
    #print(pdf_mvn(obs[i, ], mu, cov_matrix))
    #print(mu)
    #print(cov_matrix)
    ll1 <- pdf_mvn(obs[i, ], mu, cov_matrix)
    print("observation set 1")
    print(obs[i, ])
    print("this is mu")
    print(mu)
    print("this is cov")
    print(cov_matrix)
    print("this is pdf")
    print(ll1)
    
    if(ll1 == 0){
      ll1 <- 0.000000000001
    }
    ll <- ll + log(ll1)
  }
  return(-ll)  # Minimize negative log-likelihood
}
initial_params <- c(vasi_ml(ER[ , 1]), vasi_ml(ftp_2y), cor(ER[, 1], ftp_2y))
joint_data <- cbind(ER[ , 1], ftp_2y)
result <- optim(par = initial_params, fn = log_likelihood, data = joint_data, method = "L-BFGS-B", lower=c(-1, -0.05, 0.0000000001, -1, -0.05, 0.000000001, -0.5), upper = c(1, 0.5, 0.5, 1, 0.5, 0.5, 0.5))
result$par
log_likelihood(initial_params, joint_data)
initial_params
#dmvnorm(c(0.0143, 0.0022), mean = c(0.012, -0.0000000144), sigma = cbind(c(4.054355e-13, 3.011445e-13), c(3.011445e-13, 2.960252e-13)))
matrix(0, nrow = 2, ncol = 2)
(ftp_2y[2])
(ER[, 1][2])     
dmvnorm(c(0, 0), mean = c(0.012, -0.0000000144), sigma = cbind(c(1, 0.6), c(0.6, 1)))
pdf_mvn(c(0, 0), c(0.012,  -0.0000000144), cbind(c(1, 0.6), c(0.6, 1)))
#---- other method: ls calibration for test 
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
er <- ER[, 1]
#--forecast ability in test set
er_test <- er[(length(train_ftp) + 1) : length(er)]
train_er <-er[1 : length(train_ftp)]
jointpar <- ls_vasi_joint(train_er, train_ftp)
plot(ts(test_ftp), ylim = c(-0.006, 0.01), col = "red3", type = "l", lwd = 3, , ylab = "FTP in test set", xlab= "time index")
lines(vasi_forecast(jointpar[4:6], 29, train_ftp[length(train_ftp)]), col = "deeppink", lwd = 2, lty = 2)
lines(er_test, col = "seagreen", lwd = 3)
lines(vasi_forecast(jointpar[1:3], 29, train_er[length(train_er)]), col = "mediumspringgreen", lwd = 2, lty = 3)
#---individually forecast 
plot(ts(test_ftp), ylim = c(-0.006, 0.01), col = "red3", type = "l", lwd = 3, , ylab = "FTP in test set", xlab= "time index")
lines(vasi_forecast(vasi_ml(train_ftp), 29, train_ftp[length(train_ftp)]), col = "deeppink", lwd = 2, lty = 2)
lines(er_test, col = "seagreen", lwd = 3)
lines(vasi_forecast(vasi_ml(train_er), 29, train_er[length(train_er)]), col = "mediumspringgreen", lwd = 2, lty = 3)
#----individually fitting the whole sample
plot(ftp_2y, type = "l", ylim = c(-0.006, 0.01), col = "seagreen", lwd = 2, xlab = "rate", ylab = "time index")
lines(vasi_forecast(vasi_ml(ftp_2y), 120, ftp_2y[1]), col = "mediumspringgreen", lwd = 3, lty = 3)
lines(er, col  = "purple", lwd = 2)
lines(vasi_forecast(vasi_ml(er), 120, er[1]), col = "deeppink", lwd = 3, lty = 3)
legend("topright", legend = c("Actual Euribor 6M", "Actual FTP", "Vasicek Euribor", "Vasicek FTP"), col = c("purple", "seagreen", "deeppink", "mediumspringgreen"), lwd = c(2, 2, 3, 3), lty = c(1, 1, 3, 3))
# ---jointly fitting the whole sample
jointpar_all <- ls_vasi_joint(er, ftp_2y)
plot(ftp_2y, type = "l", ylim = c(-0.006, 0.01), col = "seagreen", lwd = 2, xlab = "rate", ylab = "time index")
lines(vasi_forecast(jointpar_all[4:6], 120, ftp_2y[1]), col = "mediumspringgreen", lwd = 3, lty = 3)
lines(er, col  = "purple", lwd = 2)
lines(vasi_forecast(jointpar_all[1:3], 120,  er[1]), col = "deeppink", lwd = 3, lty = 3)
legend("topright", legend = c("Actual Euribor 6M", "Actual FTP", "Vasicek Euribor", "Vasicek FTP"), col = c("purple", "seagreen", "deeppink", "mediumspringgreen"), lwd = c(2, 2, 3, 3), lty = c(1, 1, 3, 3))
par(mfrow = c(1, 2))
acf(dif_2y)
pacf(dif_2y)
par(mfrow = c(1, 1))
# ---- DF?
df6m <- 1 / (1 + ER[, 1])
df1y <- 1 / (1 + ER[, 2])
df3y <- 1 / (1 + ER[ , 3])
df5y <- 1 / (1 + ER[ , 4])
df10y <- 1 / (1 + ER[ , 5])
plot(df6m, type = "l") ## good
k <- 1 / (1 + vasi_forecast(c(a_em_ls, r_em_ls, sig_em_ls), 120, ER[ ,2 ] [1]))
lines(k)
product1 <- df6m * ftp_2y
product2 <- df1y * ftp_2y
product3 <- df3y * ftp_2y
product4 <- df5y * ftp_2y
product5 <- df10y * ftp_2y
plot(df10y, df1y)
emp_cdf_product1 <- rank(product1) / (length(product1) + 1)
emp_cdf_product2 <- rank(product2) / (length(product2) + 1)
emp_cdf_product3 <- rank(product3) / (length(product3) + 1)
emp_cdf_product4 <- rank(product4) / (length(product4) + 1)
emp_cdf_product5 <- rank(product5) / (length(product5) + 1)
plot(emp_cdf_product1, emp_cdf_product5)
cor(dif_2y, diff(er, lag = 1))
summary(lm(dif_2y~diff(er, lag = 1)))
adf.test(log(ftp_2y + 1))
adf.test(log(er + 1))
summary(lm(log(ftp_2y + 1)~log(er + 1)))
cor(diff(log(ftp_2y +1),lag = 1), diff(log(er + 1), lag = 1)) 
cor(dif_2y, diff(er, lag = 1))
plot(log(er +1 ))
summary(lm(diff(log(ftp_2y +1),lag = 1)~ diff(log(er + 1), lag = 1)))
cor(dif_2y, diff(er, lag = 1))
plot(ftp_2y, er)
plot(dif_2y, diff(er, lag = 1))
cor(ftp_2y, er, method = "spearman")
res_ftp_mod <- lm(ftp_2y[2:120]~ftp_2y[1:119])$residuals
res_er_mod <- lm(er[2:120]~er[1:119])$residuals
emp_res_er <- rank(res_er_mod) / (length(res_er_mod) + 1)
emp_res_ftp <- rank(res_ftp_mod) / (length(res_ftp_mod) + 1)
plot(emp_res_ftp, emp_res_er)
