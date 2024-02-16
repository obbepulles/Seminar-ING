
#--- Function implementing the Historical simulation of the option cost of a client given parameters:
#---    ts:         Time series data of 'used amount' from start until cancellation
#---    start:      Date class in DD-MM-YYYY format, start of contract
#---    mat:        Date class in DD-MM-YYYY format, maturity of contract
#---    cancel:     Date class in DD-MM-YYYY format, cancellation date of contract
#---    euri:       Time series data of available Euribor rate 
#---    euri_sim:   Time series data of simulated Euribor rate (dates are after last observation of 'euri' data)
#---    ftp:        Time series data of available FTP data   
#---    p:          Probability of being a type I client (non-const utilization)
#---    pool_coef:  vector of coefficients of pooled TOBIT result from main file


hist_sim_option_cost_simple <- function(ts, start, mat, cancel, euri, euri_sim, ftp, p, pool_coef){
  
  #Work in month-intervals
  T <- length(seq(from = start, to = mat, by = 'month'))-1
  tau <- length(ts)
  #length of simulated TS should be at least 2, cancelling at maturity is not allowed
  r <- (T - tau) + 1
  ts_sim <- rep(0,r)
    
  #Classify type of client (I or II), if only 1 observation, different procedure
  if(tau == 1){
    #(1) simulate a time series tau:T where U_t varying with weight p (type I) 
    #(2) simulate a time series tau:T where U_t constant with weight 1-p (type II) 
    #(3) combine into U_t* = p*TS_1 + (1-p)*TS_2 (unbiased)
    #(4) return U_t*
  }
  else if(var(ts) == 0){
    ts_sim <- rep(ts[tau], r)
  }
  else{
    #implement pooled OLS prediction if #obs < 24
    if(tau < 24){
      cat('Using Pooled TOBIT')
      ts_sim <- simulate_utilization_ts_pool(ts_sim, ts, r, tau, pool_coef)
    }
    #implement AR(1) tobit expectation of #obs >=24
    else{
      ts_sim <- simulate_utilization_ts(ts_sim, ts, r, tau)
    }
  
  }

  return(ts_sim)
}


simulate_utilization_ts_pool <- function(ts_sim, ts, r, tau, pool_coef){
  ts_sim[1:2] <- ts[(tau-1):tau]
  if(r > 2){
    for(i in 3:r){
      ts_sim[i] <- min(max(pool_coef[1] + ts_sim[i-1]*pool_coef[2],0),1)
    }
  }
  return(ts_sim)
}

simulate_utilization_ts <- function(ts_sim, ts, r, tau){
  #prepare AR(1)
  comb_ts <- na.omit(cbind(ts, lag(ts)))
  tryCatch(
    {
      coef <- coefficients(censReg(comb_ts[,1] ~ comb_ts[,2], left = 0, right = 1))
      
    },
    error = function(e) {
      coef <- coefficients(arima(ts, order = c(1,0,0), method = "ML"))
      cat("Using ARIMA")
      NA
    }
  )
  #need from t = tau-1 until T-1, tau-1 & tau are known
  ts_sim[1:2] <- ts[(tau-1):tau]
  #Simulate if necessary
  if(r > 2){
    for(i in 3:r){
      ts_sim[i] <- min(max(coef[1] + ts_sim[i-1]*coef[2],0),1)
    }
  }
  return(ts_sim)
}

