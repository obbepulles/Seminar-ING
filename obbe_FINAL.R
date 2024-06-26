library('dplyr')
library('censReg')
library('moments')
library('ggplot2')
library('ggpubr')
library("mixR")
library('forecast')
library("sn")
library('vars')
library('mixtools')
library('fitdistrplus')
library('corrplot')
library('tidyverse')

data <- readxl::read_xlsx("hypothetical_data_set.xlsx",2, skip = 1)
data <- data[,2:9]
data$`Reporting date` <- as.Date(data$`Reporting date`, "%d-%m-%Y")
data$`Maturity date` <- as.Date(data$`Maturity date`, "%d-%m-%Y")
data$`Cancellation date` <- as.Date(data$`Cancellation date`, "%d-%m-%Y")
data$`Start date` <- as.Date(data$`Start date`, "%d-%m-%Y")
rate_data <- readxl::read_xlsx("hypothetical_data_set.xlsx",1, skip = 1, range = "H2:M122")
ftp_data <- readxl::read_xlsx("hypothetical_data_set.xlsx",1, skip = 1, range = "A2:F122")
ftp_2y <- as.vector(unname(ftp_data[, 2]))[[1]]


data <- data %>% mutate(month_year = format(`Reporting date`, "%m-%Y")) 
rate_data <- rate_data %>% mutate(month_year = format(`Date`, "%m-%Y")) 
ftp_data <- ftp_data %>% mutate(month_year = format(`Date`, "%m-%Y"))

merged_df <- inner_join(data, rate_data, by = "month_year")
merged_df <- inner_join(merged_df, ftp_data, by = "month_year")
data <- merged_df[,c(-9,-10,-16)]
#
plot(x = rate_data$Date, y = rate_data$`6M`, type = 'l', col = 'blue2', 
     xlab = "Time", ylab = "Euribor Rate", ylim = c(-0.0075,0.03),
     main = "Euribor Rates")
lines(x = rate_data$Date, y = rate_data$`1Y`, type = 'l', col = 'red3')
lines(x = rate_data$Date, y = rate_data$`3Y`, type = 'l', col = 'green2')
lines(x = rate_data$Date, y = rate_data$`5Y`, type = 'l', col = 'orange2')
lines(x = rate_data$Date, y = rate_data$`10Y`, type = 'l', col = 'purple3')
legend("topright", legend = colnames(rate_data[,2:6]), col = c("blue2", "red3", "green2", "orange2", "purple3"), lty = 1, title = "Maturity")

#Initial use seems to follow uniform distr
first_utilization <- data %>%
  group_by(Client) %>% summarize(first(`Used amount`))

hist(first_utilization$`first(\`Used amount\`)`, main = "Histogram of first utilization in a contract", col = 'lightblue', xlab = "Used amount as fraction")
mean(first_utilization$`first(\`Used amount\`)`)
var(first_utilization$`first(\`Used amount\`)`)

######################################################################################################
#--- First Tobit result, only use non-constant data and length > 5
filtered_df <- data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 5)) %>%
  ungroup()

lag_df <- filtered_df %>% group_by(Client) %>% mutate(`Lag used` = lag(`Used amount`))
Y <- lag_df$`Used amount`[1:38567]
X <- lag_df$`Lag used`[1:38567]
modeldata <- na.omit(as.data.frame(cbind(Y,X)))

tobit_model <- censReg(Y ~ X, left = 0, right = 1, data = modeldata)
summary(tobit_model)
error <- modeldata$Y - modeldata$X * tobit_model$estimate[2] - tobit_model$estimate[1]

######################################################################################################
#--- Some things for plotting purposes, commented out right now
#Plot some clients utilization over time
# for(i in 0:120){
 dataSubset <- data %>% filter(Client >= 295 & Client < 300)
 
 plot <- ggplot(dataSubset, aes(x = as.Date(`Reporting date`,"%d-%m-%Y"), y = `Used amount`, group = Client, color = as.factor(Client))) +
   geom_line() +
   geom_point() +
   labs(title = "Utilization Time Series for a subset of clients",
        x = "Time",
        y = "Used Amount", color = "Clients") +
   theme_minimal()
 print(plot)
 #Sys.sleep(4)
# }

#--- Autocorrealtion part of the analysis of utilisation
#Filter out constant U_t's
filtered_df <- data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 23)) %>%
  ungroup()

lagged_use_filtered <- filtered_df %>% 
  group_by(Client) %>%
  mutate(dU = `Used amount` - lag(`Used amount`))

#Autocorrelation U_t 
autoCorr <- filtered_df %>% group_by(Client) %>% 
  summarise(Autocorrelation = acf(`Used amount`, plot = FALSE)$acf[2])
summary(autoCorr[2])
hist(autoCorr$Autocorrelation)

#Autocorrelation dU_t
lagged_use_filtered <- lagged_use_filtered %>% filter(!is.na(dU))
autoCorrdU <- lagged_use_filtered %>% group_by(Client) %>% summarise(Autocorrelation = acf(dU, plot = FALSE)$acf[2])
summary(autoCorrdU[2])
hist(autoCorrdU$Autocorrelation)
qqnorm(autoCorrdU$Autocorrelation)
qqline(autoCorrdU$Autocorrelation)

pcaf_pvals <- rep(0,10)
for(i in 1:10){
  pacf_results <- filtered_df %>%
    group_by(Client) %>%
    summarize(PACF = list(pacf(`Used amount`, plot = FALSE)$acf[i]))
  corls <- as.numeric(pacf_results$PACF)
  pcaf_pvals[i] <- t.test(x = corls, alternative = "two.sided")$p.value
}

######################################################################################################
#---Filter those clients whose contracts have already matured
#Last reporting date is 31-12-2021

#Clients whose contracts are either finished or cancelled are included
data_fin_can <- data %>% 
                group_by(Client) %>% 
                filter((last(`Reporting date`) != as.Date('31-12-2021', "%d-%m-%Y")) | (last(`Reporting date`) == last(`Maturity date`)))

hist(data_fin_can$`Used amount`, freq = FALSE)
#Filter out constant use
data_nonconst <- data_fin_can %>%
  group_by(Client) %>%
  filter(length(unique(`Used amount`)) > 1) %>%
  ungroup()

num_clients_filter <- data_nonconst %>% summarise(NumClients = n_distinct(Client))
num_clients <- data_fin_can %>% summarise(numClients = n_distinct(Client))

#Probability of using line constantly (type 2)
data_nonconst <- data_fin_can %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 11)) %>%
  ungroup()

lag_nonconst <- data_nonconst %>% group_by(Client) %>% mutate(`Lag used` = lag(`Used amount`))
Y <- lag_nonconst$`Used amount`
X <- lag_nonconst$`Lag used`
modeldata <- (as.data.frame(cbind(Y,X)))

tobit_model <- censReg(Y ~ X, left = 0, right = 1, data = modeldata)
summary(tobit_model)
error <- modeldata$X - modeldata$Y * tobit_model$estimate[2] - tobit_model$estimate[1]
lag_nonconst <- cbind(lag_nonconst, "Error" = error)
t.test((na.omit(lag_nonconst$Error)))

qqnorm(y = error)
qqline(y = error, col = 2)

#Error analysis of tobit model
c6 <- lag_nonconst %>% group_by(Client) %>% summarize(vars = sqrt(var(na.omit(Error)))) 
mod2_weibull <- mixfit(sqrt(c6$vars), ncomp = 3, family = 'weibull', max_iter = 10000)
mod2_weibull
plot(mod2_weibull)

#analyse if covid had an effect
covid_date <- as.Date("01-3-2020", "%d-%m-%Y")
lag_nonconst_covid <- lag_nonconst %>% mutate(covid = if_else(`Reporting date` >= covid_date,1,0))

#We observe no interaction between utiization and euribor / covid
tobit_model_covid <- censReg(data = lag_nonconst_covid, `Lag used` ~ `Used amount` + covid + `6M` + `1Y` + `3Y.x` + `5Y.x` + `10Y.x`, left = 0, right = 1)
summary(tobit_model_covid)

######################################################################################################
#--- Calculate probability of cancelling early
non_cancel_matured <- data %>% group_by(Client) %>% 
                        filter((is.na(last(`Cancellation date`))) & (last(`Reporting date`) == first(`Maturity date`))) %>%
                          summarize(ind = n_distinct(Client)) %>%
                            summarize(ncm = sum(ind))
cancel <- data %>% group_by(Client) %>% 
  filter(!is.na(last(`Cancellation date`))) %>%
  summarize(ind = n_distinct(Client)) %>%
  summarize(ncm = sum(ind))                  

non_matured <- data %>% group_by(Client) %>% 
  filter((is.na(last(`Cancellation date`))) & (last(`Reporting date`) != first(`Maturity date`))) %>%
  summarize(ind = n_distinct(Client)) %>%
  summarize(ncm = sum(ind))
#From the result in the paper
p_cancel <- 0.1766 
######################################################################################################
#--- Calculate probability of being type 1 (constant utilization)

non_single_utilization <- data %>% group_by(Client) %>%
                            filter(length(`Used amount`) > 1)
n_non_single_const <- non_single_utilization %>% group_by(Client) %>%
                          filter(n_distinct(`Used amount`) == 1) %>%
                            summarize(ind = n_distinct(Client)) %>% summarize(n = sum(ind))
n_non_single_var <- non_single_utilization %>% group_by(Client) %>%
                      filter(n_distinct(`Used amount`) > 1) %>%
                        summarize(ind = n_distinct(Client)) %>% summarize(n = sum(ind))
n_single <- 1200 - n_non_single_const - n_non_single_var

p_typeone <- (n_non_single_var / (n_non_single_const + n_non_single_var))[[1]]
######################################################################################################
#--- Check relation between risk factors using scatterplots, later extend to OLS/etc
cancelled_summary <- data %>% filter(!is.na(`Cancellation date`)) %>% filter(`Reporting date` == `Cancellation date`)
fraction <- cancelled_summary %>% group_by(Client) %>% summarize(fraction = ((as.numeric(`Cancellation date`)- as.numeric(`Start date`)) / (as.numeric(`Maturity date`) - as.numeric(`Start date`))))
mean_util <- data %>% filter(!is.na(`Cancellation date`)) %>% group_by(Client) %>% summarize(mean = var(`Used amount`))
plot(y = fraction$fraction, x = cancelled_summary$`Used amount`, xlab = "Utilisation at cancellation date", ylab = "Survival rate of contract")
plot(y = fraction$fraction, x = mean_util$mean, xlab = "Mean of utilisation of client", ylab = "Survival rate of contract")
plot(y = fraction$fraction, x = mean_util$mean, xlab = "Variance of utilisation of client", ylab = "Survival rate of contract")
hist(fraction$fraction)
asd <- data[,c(7,9,10,11,12,13,14,15,16,17,18)]

asd <- rename(asd, "Utilisation" = `Used amount`,
       "6M Euribor" = `6M`,
       "1Y Euribor" = `1Y`,
       "3Y Euribor" = `3Y.x`,
       "5Y Euribor" = `5Y.x`,
       "10Y Euribor" = `10Y.x`,
       "2Y FTP" = `2Y`,
       "3Y FTP" = `3Y.y`,
       "5Y FTP" = `5Y.y`,
       "7Y FTP" = `7Y`,
       "10Y FTP" = `10Y.y`)
corrplot(cor(asd), method = "shade")

x <- fraction$fraction
beta_fit <- fitdistrplus::fitdist(x, "beta")
beta_coef <- coef(beta_fit)
######################################################################################################
data2 <- readxl::read_xlsx("hypothetical_data_set.xlsx", skip = 1)
FTP <- as.data.frame(data2[, 2 : 6])
DFTP <- as.data.frame(data2[ , 1 ])
ER <- as.data.frame(data2[ , 9 : 13])
DER <- as.data.frame(data2[ , 8])
ftp_2y <- FTP[ , 1]
er <- ER[, 2]
#--functions writing here-----
jointpar_all <- ls_vasi_joint(er, ftp_2y)
######################################################################################################
#--- Procedure for finding best average ARIMA(p,d,q) model
filtered_df <- data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 23)) 

param_grid <- as.matrix(unname(expand.grid(p = 0:2, d = 0:1, q = 0:2)))
aic_vector <- rep(0,18)
n <- filtered_df %>% summarize(x = n_distinct(Client))
arima_model <- rep(0,18)
error_counter <- rep(0,18)

#--- Minimum ARIMA procedure
# for(i in 1:(length(n$Client))){
#   util <- filtered_df %>% filter(Client == n$Client[i])
#   util <- util$`Used amount`
#   for(j in 1:18){
# 
#     tryCatch(
#       {
#         arima_model[j] <- AIC(arima(util, order = param_grid[j, ], method = "ML"))
# 
#       },
#       error = function(e) {
#         error_counter[j] <<- error_counter[j] + 1
#         NA
#       }
#     )
#   }
# 
#   aic_vector <- aic_vector + arima_model
#
# }
#AR(1) model best on average
# aic_vector <- aic_vector / (length(n$Client) - error_counter)
# pos <- param_grid[which(aic_vector == min(aic_vector)), ]

######################################################################################################
#--- Option cost and EOS calculations start here
#--- Some results require manual adjusting of values for plots, we indicate where this is the case
cancelled_data <- data %>% group_by(Client) %>% filter(!is.na(`Cancellation date`))

filtered_df <- (data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 1)) %>% ungroup())[,1:7]

pooled_df <- filtered_df %>% group_by(Client) %>%
  mutate(lU = lag(`Used amount`)) %>%
  ungroup()
pooled_df <- pooled_df %>% filter(!is.na(lU))
pooled_reg <- censReg(`Used amount` ~ lU, data = pooled_df, left = 0, right = 1)
pool_coef <- coef(pooled_reg)

cancelled_data <- cancelled_data %>% filter(`Maturity date` <= as.Date("31-12-2021", "%d-%m-%Y"))

clients <- unique(as.numeric(cancelled_data$Client))

sim_joint <- vasi_joint_sim(jointpar_all, 120, er[120], ftp_2y[120])
ftp_sim_joint <- sim_joint[, 2]
euri_sim_joint <- sim_joint[, 1]

hs_costs <- rep(0,length(clients)) 
hs_eos <- rep(0,length(clients))
count <- 1
client_var <- rep(0,length(clients))

for(i in clients){
  client <- cancelled_data %>% filter(Client == i) 
  hs_costs[count] <- hist_sim_option_cost_simple(client$`Used amount`, client$`Start date`[1], client$`Maturity date`[1], client$`Cancellation date`[1], rate_data$`1Y`, euri_sim_joint, ftp_data, pool_coef)
  client_var[count] <- var(diff(client$`Used amount`))
  count <- count + 1
}

hist(hs_costs , main = "Histogram of historical simulation option costs", xlab = "NPV Cost", col = 'green4')
summary(hs_costs)

#Check hs_costs vs maturity in years
sum_data_cancel <- cancelled_data %>% group_by(Client) %>% slice(1) %>% ungroup()
sum_data_cancel <- cbind(sum_data_cancel, hs_costs)
names(sum_data_cancel)[19] <- "Option cost"

ggplot(sum_data_cancel, aes(x = Client, y = `Option cost`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Client", y = "Historical Option Cost", color = "Maturity Length (Years)") +
  theme_minimal()

ggplot(sum_data_cancel, aes(x = Maturity, y = `Option cost`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Maturity", y = "Historical Option cost", color = "Maturity Length (Years)", title = "Historical Option Cost v.s. Maturity") +
  theme_minimal()

#Get NPV for specified alpha quantile per maturity
npv_alphas <- rep(0,10)
for(i in 1:10){
  npv_alpha_iy <- sum_data_cancel %>% filter(Maturity == i) %>% reframe(sort(-`Option cost`))
  npv_alphas[i] <- npv_alpha_iy[ceiling(0.90*length(npv_alpha_iy[[1]])),1][1]
}

npv_alphas
eos <- rep(0, length(hs_costs))

#This part calculates historical EOS
count <- 1
for(i in clients){
  client <- cancelled_data %>% filter(Client == i) 
  ts <- client$`Used amount`[1:(length(client$`Used amount`) - 1)]
  mat_years <- client$Maturity[1]
  start_pos <- as.numeric(rate_data %>% summarize(which(format(client$`Start date`[1], "%m-%Y") == month_year)))
  cancel_pos <- as.numeric(rate_data %>% summarize(which(format(client$`Cancellation date`[1], "%m-%Y") == month_year)))
  rates <- rate_data$`1Y`[(start_pos + 1):cancel_pos]
  for(j in 1:length(rates)){
    rates[j] <- exp(-rates[j] * j / 12)
  }
  
  eos[count] <- npv_alphas[mat_years] / sum(rates * ts)
  count <- count + 1
}

plot(sort(eos))
summary(eos)
hist(eos)

sum_data_cancel <- cbind(sum_data_cancel, eos)
names(sum_data_cancel)[20] <- "EOS"

#Plots which serve the analysis
ggplot(sum_data_cancel, aes(x = Client, y = `EOS`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Client", y = "EOS", color = "Maturity Length (Years)") +
  theme_minimal()

ggplot(sum_data_cancel, aes(x = Maturity, y = `EOS`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Maturity", y = "EOS", color = "Maturity Length (Years)") +
  theme_minimal()



#---Check if variance determines option cost
sum_data_cancel <- cbind(sum_data_cancel, client_var)
names(sum_data_cancel)[21] <- "Client var"

ggplot(sum_data_cancel, aes(x = `Client var`, y = `Option cost`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Client Var", y = "Option cost", color = "Maturity Length (Years)") +
  theme_minimal()

ggplot(sum_data_cancel, aes(x = `Client var`, y = `EOS`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Client var", y = "EOS", color = "Maturity Length (Years)") +
  theme_minimal()

######################################################################################################
#--- Option cost & EOS for ongoing clients
par(mfrow=c(2,2))

set.seed(1)
data_ongoing <- data %>% group_by(Client) %>% 
  filter(last(`Reporting date`) != last(`Maturity date`)) %>% 
  #To replicate paper results, change the `Maturity` value in the line below to be between 1 and 10 and rerun everything from line 354 to 383
  filter(is.na(`Cancellation date`) & `Maturity` == 1)
ongoing_clients <- unique(as.numeric(data_ongoing$Client))
ongoing_option_cost_df <- matrix(0,nrow = 100, ncol = length(ongoing_clients))
ongoing_eos_df <- matrix(0,nrow = 100, ncol = length(ongoing_clients))
count <- 1
for(i in ongoing_clients){
  client <- data_ongoing %>% filter(Client == i)
  for(j in 1:100){
    #This uses the JOINT simulation for Euribor and FTP
    x <- ongoing_option_cost(client$`Used amount`, client$`Start date`[1], client$`Maturity date`[1], rate_data$`1Y`, euri_sim_joint, ftp_2y, ftp_sim_joint, pool_coef, p_cancel, p_typeone, beta_coef[1], beta_coef[2], mod2_weibull)
    ongoing_option_cost_df[j, count] <- x[[1]] * 10^6
    ongoing_eos_df[j, count] <- x[[2]] / 10^6
  }
  count <- count + 1
}

NPV_alpha_ongoing <- sort(-as.vector(ongoing_option_cost_df))[0.99 * 100 * length(ongoing_clients)]
summary(as.vector(ongoing_option_cost_df))
summary(as.vector(ongoing_eos_df*NPV_alpha_ongoing))
  NPV_alpha_ongoing
hist(ongoing_option_cost_df, col = "darkgreen",
     main = "Ongoing 1Y maturity NPV option cost",
     xlab = "NPV option cost")
hist(ongoing_eos_df*NPV_alpha_ongoing, col = "blue2",
     main = "Ongoing 1Y maturity EOS, 99% alpha",
     xlab = "EOS")
######################################################################################################
#--- Simulation of completely hypothetical clients with a completely random path for each risk driver for each client
#--- Manually change values to get medians etc for sensitivity analysis and the plots per maturity.
par(mfrow=c(5,2))
set.seed(1)
N <- 20000
NPV_matrix <- matrix(data = NA, nrow = N, ncol = 10)
EOS_denom_matrix <- matrix(data = NA, nrow = N, ncol = 10)
EOS_matrix <- matrix(data = NA, nrow = N, ncol = 10)
npv_alphas_sim <- rep(NA, 10)

#For the sensitivity analysis, we manually changed the values for each parameter and stored the median EOS 
#Into new objects called `medians1`
for(j in 1:10){  
  FTP_0 <- ftp_2y[120] + 0.0001 * (j - 2)
  for(i in 1:N){
  type_var <- rbinom(1, 1, p_typeone)
  cancel_client <- rbinom(1,1, p_cancel) 
  tau <- 1
   
  sim_joint <- vasi_joint_sim(jointpar_all, j*12, er[120], ftp_2y[120])
  ftp_sim_joint <- sim_joint[, 2]
  euri_sim_joint <- sim_joint[, 1]
  
  if(cancel_client == 1){
    tau <- rbeta(1, beta_coef[1], beta_coef[2])
  }
  
  volatility <- rmixgamma(1, mod2_weibull$pi, mod2_weibull$mu, mod2_weibull$sd)
  utilization <- simulate_client_utilization(pool_coef, tau, j, volatility, type_var)
  NPV_matrix[i, j] <- simulate_option_cost(euri_sim_joint, ftp_sim_joint, j, tau, utilization, FTP_0) * 10^6
  EOS_denom_matrix[i, j] <- simulate_eos_denom(tau, j, euri_sim_joint, utilization) / 10^6
  }
}

alpha <- 0.95
medians <- rep(NA,10)
for(j in 1:10){
  npv_alphas_sim[j] <- (sort(-NPV_matrix[, j]))[ceiling(N * alpha)] 
  EOS_matrix[, j] <- npv_alphas_sim[j] * EOS_denom_matrix[, j]  
  medians[j] <- median(EOS_matrix[, j])
}

#This calculates the difference in medians for the sensitivity analysis, manually adjusted
#(medians1 - medians) * 10^5

#Plots
hist(NPV_matrix[,1], col = "darkgreen",
     main = "Simulated 1Y maturity NPV option cost",
     xlab = "NPV option cost")
hist(EOS_matrix[,1], col = "blue2",
     main = "Simulated 1Y maturity EOS, 99% alpha",
     xlab = "EOS")
summary(NPV_matrix)
summary(EOS_matrix)
npv_alphas_sim
######################################################################################################