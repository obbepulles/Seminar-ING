library("readxl")
library('tseries')
library('Metrics')
library('sandwich')
library("aTSA")
library('forecast')
library('generics')


#--- Function implementing the Historical simulation of the option cost of a client given parameters:
#---    ts:         Time series data of 'used amount' from start until cancellation
#---    start:      Date class in DD-MM-YYYY format, start of contract
#---    mat:        Date class in DD-MM-YYYY format, maturity of contract
#---    cancel:     Date class in DD-MM-YYYY format, cancellation date of contract
#---    euri:       Time series data of available Euribor rate 
#---    euri_sim:   Time series data of simulated Euribor rate (dates are after last observation of 'euri' data)
#---    ftp:        Time series data of available FTP data for 2 year maturtiy   
#---    pool_coef:  vector of coefficients of pooled TOBIT result from main file
hist_sim_option_cost_simple <- function(ts, start, mat, cancel, euri, euri_sim, ftp, pool_coef){
  
  #Work in month-intervals
  T <- length(seq(from = start, to = mat, by = 'month')) 
  tau <- length(ts)
  #length of simulated TS should be at least 2, cancelling at maturity is not allowed
  r <- (T - tau) + 1
  ts_sim <- rep(0,r)
    
  #Classify type of client (I or II), 
  if(var(ts) == 0){
    ts_sim <- rep(ts[tau], r)
  }
  else{
    #implement pooled OLS prediction if #obs < 24
    if(tau < 24){
      #cat('Using Pooled TOBIT')
      ts_sim <- simulate_utilization_ts_pool(ts_sim, ts, r, tau, pool_coef)
    }
    #implement AR(1) tobit expectation of #obs >=24
    else{
      ts_sim <- simulate_utilization_ts(ts_sim, ts, r, tau)
    }
  
  }
  
  start_ftp <- (ftp_data %>% filter(format(start, "%m-%Y") == month_year))[2]
  cancel_ftp <- (ftp_data %>% filter(format(cancel, "%m-%Y") == month_year))[2]
  
  ftp_0 <- start_ftp + 0.0001 * (T / 12 - 2)
  ftp_tau <- cancel_ftp + 0.0001 * (r / 12 - 2)
  
  dFTP_tau <- as.numeric(ftp_tau - ftp_0)
  start_pos <- as.numeric(rate_data %>% summarize(which(format(start, "%m-%Y") == month_year)))
  cancel_pos <- as.numeric(rate_data %>% summarize(which(format(cancel, "%m-%Y") == month_year)))
  stop_pos <- cancel_pos + r - 1
  
  rates <-  c(euri,euri_sim)[cancel_pos:stop_pos]
  for(i in 1:length(rates)){
    rates[i] <- exp(-rates[i] * (i + tau) / 12)
  }
  
  option_cost <- dFTP_tau * sum(ts_sim * rates)
  return(option_cost[[1]])
}

#--- Simulate future expected utilization based on the POOLED OLS TOBIT model (if few observations for a client are available)
simulate_utilization_ts_pool <- function(ts_sim, ts, r, tau, pool_coef){
  ts_sim[1:2] <- ts[(tau-1):tau]
  comb_ts <- na.omit(cbind(ts, lag(ts)))
  sigmasq <- exp(pool_coef[3])^2
  if(r > 2){
    for(i in 3:r){
      c <- pool_coef[1] + ts_sim[i-1] * pool_coef[2]
      ts_sim[i] <- min(max((c * (pnorm(1 - c, mean = 0, sd = sqrt(sigmasq)) - pnorm(-c,mean = 0, sd = sqrt(sigmasq))) + 1 - pnorm(1 - c, mean = 0, sd = sqrt(sigmasq))), 0), 1)
      
    }
  }
  return(ts_sim)
}

#--- Simulate future expected utilization based on AR(1) (TOBIT) model for the client based on their past utilization
simulate_utilization_ts <- function(ts_sim, ts, r, tau){
  
  comb_ts <- na.omit(cbind(ts, lag(ts)))
  coef_arima <- tryCatch(
    {
      coef_arima <- coef(censReg(comb_ts[,1] ~ comb_ts[,2], left = 0, right = 1))
    },
    error = function(e) {
      coef_arima <- coef(Arima(ts, order = c(1,0,0), method = "ML"))
      coef_arima
    }
  )
  
  error <- comb_ts[,1] - coef_arima[1] - coef_arima[2]*comb_ts[,2]
  sigmasq <- var(error)
  #Need data from t = tau-1 until T-1, t = tau-1 & t = tau are known
  ts_sim[1:2] <- ts[(tau-1):tau]
  #Simulate the rest
  if(r > 2){
    for(i in 3:r){
      c <- coef_arima[1] + ts_sim[i-1] * coef_arima[2]
      ts_sim[i] <- min(max((c * (pnorm(1 - c, mean = 0, sd = sqrt(sigmasq)) - pnorm(-c,mean = 0, sd = sqrt(sigmasq))) + 1 - pnorm(1 - c, mean = 0, sd = sqrt(sigmasq))), 0), 1)
    }
  }
  return(ts_sim)
}

#--- Compute option cost of ongoing client
#--- ts:        Time series of utilization until latest date
#--- start:     Start date
#--- mat:       Maturity date
#--- euri:      Euribor rates until December 2021
#--- euri_sim:  Simulated Euribor rates (01-2022 --- 01-2032)
#--- ftp:       FTP rates until December 2021
#--- ftp_sim:   Simulated FTP (01-2022 --- 01-2032)
#--- 
ongoing_option_cost <- function(ts, start, mat, euri, euri_sim, ftp, ftp_sim, pool_coef, p_cancel, p_typeone, beta_a, beta_b, sdmodel){
  
  #Work in month-intervals
  T <- length(seq(from = start, to = mat, by = 'month')) - 1
  #Simulate the rest of the time utilizaiton time series for this client
  sim <- simulate_utilization_ongoing(start, mat, ts, pool_coef, p_cancel, p_typeone, beta_a, beta_b, sdmodel)
  ts_sim <- sim[[1]]
  tau <- sim[[2]]
  start_pos <- as.numeric(rate_data %>% summarize(which(format(start, "%m-%Y") == month_year)))
  cancel_pos <- start_pos + tau
  stop_pos <- start_pos + T
  
  start_ftp <- ftp[start_pos]
  cancel_ftp <- (c(ftp, ftp_sim))[cancel_pos]
  
  ftp_0 <- (start_ftp) + 0.0001 * (T / 12 - 2)
  ftp_tau <- (cancel_ftp) + 0.0001 * ((T-tau) / 12 - 2)
  
  if(tau == T){
    ftp_tau = 0
  }
  
  dFTP_tau <- as.numeric(ftp_tau - ftp_0)
  
  rates <-  c(euri,euri_sim)[cancel_pos:stop_pos]
  
  for(i in 1:length(rates)){
    rates[i] <- exp(-rates[i] * (i+tau) / 12)
  }
  u <- c(ts, ts_sim)[(tau-1):(T-1)]
  option_cost <- dFTP_tau * sum(u * rates)

  ratesn <- c(euri,euri_sim)[(start_pos+1):(cancel_pos)]
  for(i in 1:length(ratesn)){
    ratesn[i] <- exp(-ratesn[i] * i / 12)
  }
  eos_denom <- 1 / sum(c(ts,ts_sim)[1:(tau)] * ratesn)
  
  return(list(option_cost, eos_denom))
}

#--- Simulate the utilization for a contract which has not ended yet at the "end" date of the data
simulate_utilization_ongoing <- function(start, mat, ts, pool_coef, p_cancel, p_typeone, beta_a, beta_b, sdmodel){
  
  T <- length(seq(from = start, to = mat, by = 'month')) - 1
  tau <- length(ts)
  r <- T - tau 
  ts_sim <- rep(0,r)
  taustar <- T
  constuse <- TRUE
  
  p_t_less_tau <- pbeta(tau / T, shape1 = beta_a, shape2 = beta_b)
  #Given that the contract survived until time Tau, what is the probability of cancelling in (tau,T)?
  p_cancel_cond <- 1 - (1 - p_cancel) / (1 - p_cancel * p_t_less_tau)
  #Generate the type according to the probability of being type 1 or type 0
  #--- type 0: don't cancel before maturity
  #--- type 1: cancel before maturity
  type <- rbinom(1, 1, prob = p_cancel_cond)
  
  if(length(ts) > 1){
    if(var(ts) > 0){
      constuse <- FALSE
    }
  }
  else if(length(ts) == 1){
    constuse <- ifelse((rbinom(1, 1, prob = (1-p_typeone)) == 1), TRUE, FALSE)
  }
  
  if(tau < 24){
    sigma <- rmixweibull(1, sdmodel$pi, sdmodel$mu, sdmodel$sd)
  }
  else{
    comb_ts <- na.omit(cbind(ts,lag(ts)))
    error <- comb_ts[,1] - pool_coef[1] - pool_coef[2]*comb_ts[,2]
    sigma <- sd(error)
  }
  #If use is varying then make a distinction between clients who cancel and those who dont
  if(type == 1){
    taustar <- simulate_cond_cancel(beta_a, beta_b, tau, T)
    #Simulate random use until taustar, then use expectation for the rest
    
    #If we have constant use then the simulated ts is just the first observation as use does not change over time
    if(constuse){
      ts_sim <- rep(ts[1],r)
      return(list(ts_sim, taustar))
    }
    
    ts_sim[1] <- max(min(pool_coef[1] + pool_coef[2]*ts[tau] + rnorm(1, mean = 0, sd = sigma), 1), 0)
    
    #Until the new cancellation time
    for(i in 2:(taustar - tau)){
      ts_sim[i] <- max(min(pool_coef[1] + pool_coef[2]*ts_sim[i-1] + rnorm(1, mean = 0, sd = sigma), 1), 0)
    }
    #From new cancellation until maturity
    for(i in (taustar - tau + 1):r){
      c <- pool_coef[1] + ts_sim[i-1] * pool_coef[2]
      ts_sim[i] <- max(min((c * (pnorm(1 - c, mean = 0, sd = sigma) - pnorm(-c,mean = 0, sd = sigma)) + 1 - pnorm(1 - c, mean = 0, sd = sigma)), 1), 0)
    }
    
  }
  if(type == 0){
    
    if(constuse){
      ts_sim <- rep(ts[1],r)
      return(list(ts_sim,taustar))
    }
    
    ts_sim[1] <- max(min(pool_coef[1] + pool_coef[2]*ts[tau] + rnorm(1, mean = 0, sd = sigma), 1), 0)
    for(i in 2:r){
      ts_sim[i] <- max(min(pool_coef[1] + pool_coef[2]*ts_sim[i-1] + rnorm(1, mean = 0, sd = sigma), 1), 0)
    }
  }
  
  return(list(ts_sim,taustar))
}

#--- Computes conditional CDF value 
cancel_cond_cdf <- function(taustar, beta_a, beta_b, T, tau){
  num <- pbeta(taustar / T, shape1 = beta_a, shape2 = beta_b) - pbeta(tau / T, shape1 = beta_a, shape2 = beta_b)
  denom <- 1 - pbeta(tau / T, shape1 = beta_a, shape2 = beta_b)
  return(num / denom)
}

#--- Simulates a new cancellation date for a given contract
simulate_cond_cancel <- function(beta_a, beta_b, tau, T){
  u <- runif(1, min = 0, max = 1)
  new_cancel <- (tau + T) / 2
  opt <- optim(new_cancel, fn <- function(x){ ((u-cancel_cond_cdf(x, beta_a, beta_b, T, tau))^2) }, method = "Brent", lower = tau, upper = T)
  return(ceiling(opt$par))
}

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

#--- Simulates a option cost of a completely hypothetical client
#--- GIVEN INFO: FTP_0 & RATE_0 
#--- VARIABLES:
#--- rates:         Simulated rates                         (1Y Maturity)
#--- ftp:           Simulated FTP rates                     (2Y Maturity)
#--- maturity:      Maturity of contract in YEARS           (1, ..., 10)
#--- cancellation:  When the contract is cancelled          (DECIMAL or 1)
#--- Utilization:   Simulated Utilization Time Series       (From initialization until end)
#--- FTP_0:         FTP rate at initiation of the contract  (31-12-2021 KNOWN)
simulate_option_cost <- function(rates, ftp, maturity, cancellation, utilization, FTP_0){
  T <- maturity * 12
  tau <- ceiling(cancellation * T)
  if(tau == 1){
    tau <- 2
  }
  dftp_tau <- (ftp[tau] + 0.0001 * ((T-tau) / 12 - 2)) - FTP_0 
  if(tau == T){
    dftp_tau <- -FTP_0 
  }
  rates <- rates[tau:T]
  for(i in 1:length(rates)){
    rates[i] <- exp(-rates[i] * (i+tau)/12)
  }
  u_after_tau <- utilization[(tau-1):(T-1)]
  if(length(u_after_tau) != length(rates)){
    print(length(u_after_tau))
    print(length(rates))
    print(tau)
    print(T)
  }
  option_cost <- dftp_tau * sum(u_after_tau * rates) 
  return(option_cost)
}

#--- Simulates a hypothetical client
#--- pooled_coefs:  Coefficients from the POOLED TOBIT OLS
#--- cancellation:  When the contract is cancelled              (DECIMAL or 1)
#--- maturity:      Maturity of contract in YEARS               (1, ..., 10)
#--- volatility:    Volatility(SD) of the client                (RANDOM DRAW FROM ONE OF THREE GAMMA DISTRIBUTIONS)
#--- type:          Binary, 0 = constant use, 1 = variable use
simulate_client_utilization <- function(pooled_coefs, cancellation, maturity, volatility, type){
  T <- maturity * 12
  tau <- ceiling(cancellation * T)
  u <- rep(0,T)
  u[1] <- runif(1) 
  if(tau == 1){
    tau <- 2
  }
  
  if(type == 0){
    u <- rep(u[1],T)
    return(u)
  }
  else{
    errors <- rnorm((tau - 1), mean = 0, volatility)
    #Simulate until maturity
    for(i in 2:tau){
      u[i] <- min(max((pooled_coefs[1] + pooled_coefs[2] * u[i - 1] + errors[i-1]), 0), 1)
    }
    if(T > tau){
      for(i in (tau+1):T){
        c <- pooled_coefs[1] + pooled_coefs[2] * u[i - 1]
        u[i] <- min(max((c * (pnorm(1 - c, mean = 0, sd = volatility) - pnorm(-c,mean = 0, sd = volatility)) + 1 - pnorm(1 - c, mean = 0, sd = volatility)), 0), 1)
      }
    }
  }
  return(u)
}

#--- Simulates EOS for a hypothetical client
#--- cancellation:          When the contract is cancelled            (DECIMAL or 1)
#--- maturity:              Maturity of contract in YEARS             (1, ..., 10)
#--- rates:                 Simulated rates                           (From initialization until end)
#--- utilization:           Simulated Utilization Time Series         (From initialization until end)
simulate_eos_denom <- function(cancellation, maturity, rates, utilization){
  tau <- ceiling(cancellation * T)
  if(tau == 1){
    tau <- 2
  }
  rates <- rates[1:tau]
  for(i in 1:length(rates)){
    rates[i] <- exp(-rates[i] * i/12)
  }
  u <- utilization[1:tau]
  eos <- 1 / sum(rates * utilization)
  return(eos)
}



