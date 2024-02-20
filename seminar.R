library("readxl")
library(tseries)
library(Metrics)
library(sandwich)
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
    print(sim[i] == expected_forecast[i])
  }
  return(sim)
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
par(mfrow = c(1, 2))
acf(ftp_2y)
pacf(ftp_2y)
#------ with consideration of interest rate leave it later can be ignored now 
FTP_i <- cbind((FTP[, 1] - 0.001), FTP[ , 2:3], FTP[ , 5])
ER_i <- ER[ , 2:5]
plot(ftp_2y , type = "l" , ylim = c(-0.01 , 0.02))
lines(ER[ , 1], col = "blue")
plot(diff(ftp_2y), type = "l" , ylim = c(-0.01 , 0.01))
lines(diff(ER[ , 1]) , col = "blue")
pp <- lm(diff(ftp_2y)~diff(ER [ , 1]))
pp
plot(diff(ftp_2y), type = "l")
#------ use for autocorrelation it need to be garch later on, it shows high autocorrelation
lines(diff(ER[ , 1]), col = "red")
kk <- lm(ftp_2y~ ER[ , 1])
library("aTSA")
adf.test(diff(ftp_2y))
library(tseries)
print(acf(ftp_2y^2))
#------- trying to do normal draws from dif_ftp to predict
plot(density(diff(ftp_2y)))
library(mixR)
mix_norm1 <- mixfit(dif_2y, ncomp= 1, family = "norm")
plot(mix_norm1)
xxx<- dif_2y[2:length(dif_2y)] # xt
xx <- dif_2y[1:(length(dif_2y)-1)] #xt-11
length(xxx)
y <- lm(xxx~xx)
summary(y)
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
##-----AR(1) model
library(lmtest)
AR1 <- lm(train_ftp[2 : length(train_ftp)] ~ train_ftp[1 : (length(train_ftp) - 1)])
AR_mean_res <- mean(AR1$residuals)
AR_sd_res <- sd(AR1$residuals)
plot(AR1$residuals)
summary(AR1)
bptest(AR1)
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
##-------random walk 
rw <- rep(ftp_2y[train + 1], test)
plot(test_ftp, type = "l")
lines( rw, col = "blue")
rmse_rw <- rmse(test_ftp, rw)
rmse_normal <- rmse(test_dif, norm_predict)
plot(norm_predict, type = "l")
plot(test_dif, type = "l") ## this is very bad
pacf(dif_2y)
acf(dif_2y)
#----ARMA(1, 1)
library(forecast)
train_ftp_ts <- ts(train_ftp, start = c(2012, 1), frequency = 12)
ftp_2y_ts <- ts(ftp_2y,start = c(2012, 1), frequency = 12 )
arma11 <- Arima(train_ftp_ts, order = c(1, 0, 1), method = "CSS", include.constant = T)
sim_ar11 <- forecast(arma11, h = 29)
plot(test_ftp)
plot(sim_ar11)
lines(ftp_2y_ts, col = "black")
lines(sim_ar11$fitted, col = "cornflowerblue", lty = 2, lwd = 1.2)
legend("bottomleft", legend = c("fitted value","actual value", "predicted value"), col = c("cornflowerblue","black", "cyan2"), lty = c(2 , 1, 1), lwd = c(1, 1, 2))
plot(sim_ar11$residuals)
plot(predict(sim_ar11))
rmse_arma11 <- rmse(test_ftp, sim_ar11$mean)
rmse_arma11
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
arima0112 <- Arima(train_ftp_ts, order = c(0, 1, 12), include.mean = T, method = "CSS")
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
#---- auto ARMA model 
auto_arma_model <- auto.arima(train_ftp_ts)
sim_arma_auto <- forecast(auto_arma_model, h = 29)
plot(sim_arma_auto, main = "Forecasts from Random Walk")
lines(ftp_2y_ts, col = "black")
lines(auto_arma_model$fitted, col = "cornflowerblue", lty = 2, lwd = 1.2)
legend("bottomleft", legend = c("fitted value","actual value", "predicted value"), col = c("cornflowerblue","black", "cyan2"), lty = c(2 , 1, 1), lwd = c(1, 1, 2))
rmse_arma_auto <- rmse(test_ftp, sim_arma_auto$mean)
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
#---- normal muxture of ftp directly
mix_norm1_ftp <- mixfit(ftp_2y, ncomp= 1, family = "norm")
mix_norm2_ftp <- mixfit(ftp_2y, ncomp= 2, family = "norm")
mix_norm3_ftp <- mixfit(ftp_2y, ncomp= 3, family = "norm")
mix_norm4_ftp <- mixfit(ftp_2y, ncomp= 4, family = "norm")
mix_norm1_ftp$aic
mix_norm2_ftp$aic
mix_norm3_ftp$aic
mix_norm4_ftp$aic
plot(mix_norm3_ftp)
par(mfrow = c(2 , 1))
acf(dif_2y)
pacf(dif_2y)
#----GARCH directly to FTP
# Assuming 'returns' is your financial return time series
library(rugarch)
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)))
fit <- ugarchfit(spec, data = train_ftp_ts)
sim_GARCH <- ugarchforecast(fit, n.ahead = 29)
#fitted(sim_GARCH)[1 : 29]
rmse_GARCH <- rmse(test_ftp, fitted(sim_GARCH))
plot(sim_GARCH)
lines(ftp_2y_ts)
#-----Kalman Filter
library(FKF)
y <- train_ftp_ts
a0 <- y[1]
P0 <- matrix(1)
dt <- ct <- matrix(0)
Zt <- Tt <- matrix(1)
fit.fkf <- optim(c(HHt = var(y) * .5,
                   GGt = var(y) * .5),
                 fn = function(par, ...)
                   -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                 yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                 Zt = Zt, Tt = Tt)
HHt <- as.numeric(fit.fkf$par[1])
GGt <- as.numeric(fit.fkf$par[2])
y_fkf <- fkf(a0, P0, dt, ct, Tt, Zt,
             HHt = matrix(HHt), GGt = matrix(GGt),
             yt = rbind(y))
y_kalman = as.numeric(y_fkf$att)
par(mfrow = c(1, 1))
plot(train_ftp, type = "l")
lines(y_kalman, col = "blue")
#----- check Hamilton for FTP
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
while(count < 100000){    
    #M-step 
  p11 <- (sum(p_ij_t[1 , ]) - p_ij_t[1 , 1]) / (sum(p_ij_t[1 , ]) - p_ij_t[1, T] + sum(p_ij_t[2, ]) - p_ij_t[2, T])
  p22 <- (sum(p_ij_t[4 , ]) - p_ij_t[4 , 1]) / (sum(p_ij_t[4 , ]) - p_ij_t[4, T] + sum(p_ij_t[3, ]) - p_ij_t[3, T])
  mu_1 <- sum((p_ij_t[1 , ] + p_ij_t[2, ]) * ftp_2y) / (sum(p_ij_t[1 , ]) + sum(p_ij_t[2 , ]))
  mu_2 <- sum((p_ij_t[4 , ] + p_ij_t[3, ]) * ftp_2y) / (sum(p_ij_t[4 , ]) + sum(p_ij_t[3 , ]))
  sig_1 <- sum((p_ij_t[1 , ] + p_ij_t[2, ]) * (ftp_2y - mu_1)^2) /  (sum(p_ij_t[1 , ]) + sum(p_ij_t[2 , ]))
  sig_2 <- sum((p_ij_t[4 , ] + p_ij_t[3, ]) * (ftp_2y - mu_2)^2) /  (sum(p_ij_t[4 , ]) + sum(p_ij_t[3 , ]))
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
sp2 <- (T-1)*(sig_1^2 + sig_2^2) / (2*T - 2)
t_test <- (mu_1 - mu_2)/sqrt(sp2/T + sp2/T)
t_test <- (mu_1 - mu_2) / sqrt(sig_1 / T + sig_2 / T)
t_test
p_value <- 2 * (1 - pt(abs(t_test),  ( 2*T))) # this shows there is no effect 
p_value
#--- function of dFTP 
dFTP_func <- function(tau, data, T){
  #T for maturity , tau is cancellation time
  
}
summary(lm(dif_2y~ diff(ER_i[ , 1])))
#-- vasicek on dftp 
n2 <- length(train_dif)
a_vasi2 <- -1 / 12 * log ((n2 * sum(train_dif[1: (n2 - 1)] * train_dif[2 : n2]) - sum(train_dif[2 : n2]) * sum(train_dif[ 1: (n2 - 1)]))/(n2 * sum(train_dif[1 : (n2 - 1)]^2) - sum(train_dif[1 : (n2 - 1)])^2))
lr_mu_vasi2 <- 1 / (n2 * (1 - exp(-a_vasi2 / 12))) * (sum(train_dif[2 : n2]) - exp(a_vasi2 / 12) * sum(train_dif[1 : (n2 - 1)]))
var_vasi2 <- 2 * a_vasi2 / (n2 * (1 - exp(-a_vasi2 / 6))) * sum((train_dif[2 : n] - train_dif[1 : (n -1)] * exp(-a_vasi2/12) - lr_mu_vasi2 * (1 - exp(-a_vasi2/12)))^2)
expected_forecast2 <- rep(0 , length(test_dif))
expected_forecast2 [1] <- train_dif[n2] * exp(-a_vasi2) + lr_mu_vasi2 * (1 - exp(-a_vasi2 / 12))
for(i in 2 : test){
  expected_forecast[i] <- expected_forecast[(i - 1)] * exp(-a_vasi) + lr_mu_vasi * (1 - exp(-a_vasi / 12))
}
plot(test_dif, type = "l")
plot(expected_forecast, col = "blue")
rmse_vasi_paper2 <- rmse(test_dif, expected_forecast)
#-----
shapiro.test(sarima01012$residuals, type = "Ljung-Box")
Box.test(sarima01012$residuals, type = "Ljung-Box")
mu <- mean(sarima01012$residuals)
#------ simulation for different model: # only generate n value tho
#1. SARIMA(1, 0, 0)(0, 0, 1)[4]; 2. GARCH(1, 1); 3. AR(1); 4. ARIMA(0, 1, 1) 5. RW; 6. Vasicek
ftp_forecast <- function(n, type){
  if(type == 1){
    sample_sarima <- arima(ftp_2y_ts, order = c(1, 0, 0), seasonal = c(0,0,1,4), method = "ML")
    sim_1 <- forecast(sample_sarima, h = n)
    return(sim_1$mean)
  }
  else if(type == 2){
    spec1 <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)))
    fit1 <- ugarchfit(spec1, data = ftp_2y_ts)
    sim_2 <- ugarchforecast(fit1, n.ahead = n)
    return(fitted(sim_2))
  }
  else if(type == 3){
    sample_sarima <- arima(ftp_2y_ts, order = c(1, 0, 0), method = "ML")
    sim_1 <- forecast(sample_sarima, h = n)
    return(sim_1$mean)
  }
  else if(type == 4){
    sample_sarima <- arima(ftp_2y_ts, order = c(0, 1, 1), method = "ML")
    sim_1 <- forecast(sample_sarima, h = n)
    return(sim_1$mean)
  }
  else if(type == 5){
    res <- rnorm(n, mean = mean(dif_2y), sd = sd(dif_2y))
    sim <- rep(0, n)
    sim[1] <- res[1] + ftp_2y_ts[length(ftp_2y_ts)]
    for(i in 2 : n){
      sim[i] <- res[i] + sim[(i - 1)] 
    }
    return(sim)
  }
  else if (type == 6){
    return(0)
  }
}

