
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
  T <- length(seq(from = start, to = mat, by = 'month')) - 1
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
    rates[i] <- exp(-rates[i] * i / 12)
  }
  
  option_cost <- dFTP_tau * sum(ts_sim * rates)
  print(option_cost)
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
  if(tau == T){
    return(0)
  }
  
  r <- T - tau + 1
  
  start_ftp <- (ftp_data %>% filter(format(start, "%m-%Y") == month_year))[2]
  cancel_ftp <- ftp_sim[tau]
  
  ftp_0 <- start_ftp + 0.0001 * (T / 12 - 2)
  ftp_tau <- cancel_ftp + 0.0001 * (r / 12 - 2)
  
  dFTP_tau <- as.numeric(ftp_tau - ftp_0)
  start_pos <- as.numeric(rate_data %>% summarize(which(format(start, "%m-%Y") == month_year)))
  cancel_pos <- start_pos + tau
  stop_pos <- start_pos + T
  
  rates <-  c(euri,euri_sim)[cancel_pos:stop_pos]
  for(i in 1:length(rates)){
    rates[i] <- exp(-rates[i] * i / 12)
  }
  
  u <- c(ts,ts_sim)[tau:T]
  option_cost <- dFTP_tau * sum(u * rates)

  return(option_cost[[1]])
}

#--- Simulate the utilization for a contract which has not ended yet at the "end" date of the data
simulate_utilization_ongoing <- function(start, mat, ts, pool_coef, p_cancel, p_typeone, beta_a, beta_b, sdmodel){
  
  T <- length(seq(from = start, to = mat, by = 'month')) - 1
  tau <- length(ts)
  r <- T - tau + 1
  ts_sim <- rep(0,r)
  taustar <- r
  constuse <- TRUE
  if(length(ts) > 1){
    if(var(ts) > 0){
      constuse <- FALSE
    }
  }
  else if(length(ts) == 1){
    constuse <- ifelse((rbinom(1, 1, prob = (1-p_typeone)) == 1), TRUE, FALSE)
  }
  #If we have constant use then the simulated ts is just the first observation as use does not change over time
  if(constuse){
    ts_sim <- rep(ts[1],r)
    return(list(ts_sim,r))
  }
  
  p_t_less_tau <- pbeta(tau / T, shape1 = beta_a, shape2 = beta_b)
  #Given that the contract survived until time Tau, what is the probability of cancelling in (tau,T)?
  p_cancel_cond <- 1 - (1 - p_cancel) / (1 - p_cancel * p_t_less_tau)
  #Generate the type according to the probability of being type 1 or type 0
  #--- type 0: don't cancel before maturity
  #--- type 1: cancel before maturity
  type <- rbinom(1, 1, prob = p_cancel_cond)
  
  if(tau < 24){
    sigma <- rmixgamma(1, sdmodel$pi, sdmodel$mu, sdmodel$sd)
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
    ts_sim[1] <- max(min(pool_coef[1] + pool_coef[2]*ts[tau] + rnorm(1, mean = 0, sd = sigma), 1), 0)
    for(i in 2:r){
      ts_sim[i] <- max(min(pool_coef[1] + pool_coef[2]*ts_sim[i-1] + rnorm(1, mean = 0, sd = sigma), 1), 0)
    }
  }
  print(taustar)
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

