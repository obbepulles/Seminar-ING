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

