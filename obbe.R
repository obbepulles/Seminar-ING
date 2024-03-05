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

#Initial use seems to follow uniform distr
first_utilization <- data %>%
  group_by(Client) %>% summarize(first(`Used amount`))

hist(first_utilization$`first(\`Used amount\`)`, main = "Histogram of first utilization in a contract", col = 'lightblue', xlab = "Used amount as fraction")


#(1) We see two types of clients:
#       type 1: Use credit line once, don't touch again (~20% of clients)
#       type 2: Use credit line more often, high autocorrelation (Tobit model?, ~80% of clients)
#(2) Initial amount used seems to be Unif(0,1)
#(3) No ACF in dU
#(4) For type 2 clients, dU seems to follow a normal distr around (-0.2;0.2) but has thick tails
#    Ideas: Mixture of normals / Truncated normal + other distr for tails
#    Problem: U in [0,1], thick tails make sense, could be a result of TOBIT (LIKELY!!!)

#DGP idea: 
#(1) U_0 ~ Unif[0,1]
#(2) eps_t ~ Norm(mu,sigma^2)
#(3) split up U_t into cases: U_t = 
#     (i)    1              , if phi * U_{t-1} + eps_t >= 1
#     (ii)   0              , if phi * U_{t-1} + eps_t <= 0
#     (iii)  phi * U_{t-1} + eps_t, otherwise

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

#plot(density(error^2))
#qqnorm(y = error)
#qqline(y = error, col = 2)
######################################################################################################

#--- Some random things for plotting purposes
modes <- function(d){
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

#Plot some clients utilization over time
# for(i in 0:120){
# dataSubset <- data %>% filter(Client >= 10*i & Client < (10*(i+1)))
# 
# plot <- ggplot(dataSubset, aes(x = as.Date(`Reporting date`,"%d-%m-%Y"), y = `Used amount`, group = Client, color = as.factor(Client))) +
#   geom_line() +
#   geom_point() +
#   labs(title = "Utilization Time Series for a subset of clients",
#        x = "Client",
#        y = "Used Amount") +
#   theme_minimal()
# print(plot)
# Sys.sleep(4)
# }


######################################################################################################
#--- Most up-to-date data analysis starts
#Filter those clients whose contracts have already matured
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
p <- num_clients_filter$NumClients / sum(num_clients$numClients)

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
c6 <- lag_nonconst %>% group_by(Client) %>% summarize(vars = var(na.omit(Error))) 
plot(hist(c6$vars))

x <- sqrt(c6$vars)
mod2_weibull <- mixfit(x, ncomp = 3, family = 'gamma', max_iter = 10000)
mod2_weibull
plot(mod2_weibull)

mod2_weibull
#analyse if covid had an effect
covid_date <- as.Date("01-3-2020", "%d-%m-%Y")
lag_nonconst_covid <- lag_nonconst %>% mutate(covid = if_else(`Reporting date` >= covid_date,1,0))

#We observe no interaction between utiization and euribor / covid
tobit_model_covid <- censReg(data = lag_nonconst_covid, `Lag used` ~ `Used amount` + covid + `6M` + `1Y` + `3Y.x` + `5Y.x` + `10Y.x`)
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
p_cancel <- (cancel / (non_cancel_matured + cancel))[[1]]
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
mean_util <- data %>% filter(!is.na(`Cancellation date`)) %>% group_by(Client) %>% summarize(mean = mean(`Used amount`))
plot(y = fraction$fraction, x = cancelled_summary$`Used amount`)
plot(y = fraction$fraction, x = mean_util$mean)
hist(fraction$fraction)

x <- fraction$fraction
beta_fit <- fitdistrplus::fitdist(x, "beta")
beta_coef <- coef(beta_fit)
######################################################################################################
filtered_df <- data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 23)) 

param_grid <- as.matrix(unname(expand.grid(p = 0:2, d = 0:1, q = 0:2)))
aic_vector <- rep(0,18)
n <- filtered_df %>% summarize(x = n_distinct(Client))
arima_model <- rep(0,18)
error_counter <- rep(0,18)

for(i in 1:(length(n$Client))){
 util <- filtered_df %>% filter(Client == n$Client[i])
 util <- util$`Used amount`
 for(j in 1:18){
   
   tryCatch(
     {
       arima_model[j] <- AIC(arima(util, order = param_grid[j, ], method = "ML"))
       
     },
     error = function(e) {
       error_counter[j] <<- error_counter[j] + 1
       NA
     }
   )
 }
 
 aic_vector <- aic_vector + arima_model
 
}

#AR(1) model best on average
aic_vector <- aic_vector / (length(n$Client) - error_counter)
pos <- param_grid[which(aic_vector == min(aic_vector)), ]


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

varmod <- VAR(ts(rate_data[,2:6]),p = 1)
varcoef <- coef(varmod)
forecast <- predict(varmod, n.ahead = 120)
euri_sim <- forecast[["fcst"]][["X1Y"]]
euri_sim <- euri_sim[,1]


hs_costs <- rep(0,length(clients)) 
hs_eos <- rep(0,length(clients))
count <- 1
client_var <- rep(0,length(clients))

for(i in clients){
  client <- cancelled_data %>% filter(Client == i) 
  hs_costs[count] <- hist_sim_option_cost_simple(client$`Used amount`, client$`Start date`[1], client$`Maturity date`[1], client$`Cancellation date`[1], rate_data$`1Y`, euri_sim, ftp_data, pool_coef)
  client_var[count] <- var(diff(client$`Used amount`))
  
  count <- count + 1
  
}

hist(hs_costs , main = "Histogram of historical simulation option costs", xlab = "NPV Cost", col = 'green4')
summary(hs_costs)
#check hs_costs vs maturity in years

sum_data_cancel <- cancelled_data %>% group_by(Client) %>% slice(1)  
sum_data_cancel <- cbind(sum_data_cancel, hs_costs)
names(sum_data_cancel)[19] <- "Option cost"

ggplot(sum_data_cancel, aes(x = Client, y = `Option cost`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Client", y = "Historical Option Cost", color = "Maturity Length (Years)") +
  theme_minimal()

ggplot(sum_data_cancel, aes(x = Maturity, y = `Option cost`, color = factor(Maturity))) +
  geom_point() +
  labs(x = "Maturity", y = "Historical Option cost", color = "Maturity Length (Years)") +
  theme_minimal()


npv_alpha <- sort(-hs_costs)[ceiling(0.99 * length(hs_costs))]
eos <- rep(0, length(hs_costs))

count <- 1
for(i in clients){
  client <- cancelled_data %>% filter(Client == i) 
  ts <- client$`Used amount`[1:(length(client$`Used amount`) - 1)]
  
  start_pos <- as.numeric(rate_data %>% summarize(which(format(client$`Start date`[1], "%m-%Y") == month_year)))
  cancel_pos <- as.numeric(rate_data %>% summarize(which(format(client$`Cancellation date`[1], "%m-%Y") == month_year)))
  rates <- rate_data$`1Y`[(start_pos + 1):cancel_pos]
  for(j in 1:length(rates)){
    rates[j] <- (1 + rates[j]/12) ^ (-j)
  }
  
  eos[count] <- npv_alpha / sum(rates * ts)
  count <- count + 1
}

plot(sort(eos))
summary(eos)
hist(eos)

sum_data_cancel <- cbind(sum_data_cancel, eos)
names(sum_data_cancel)[20] <- "EOS"

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
data_ongoing <- data %>% group_by(Client) %>% filter(last(`Reporting date`) != last(`Maturity date`)) %>% filter(is.na(`Cancellation date`))
ongoing_clients <- unique(as.numeric(data_ongoing$Client))
ongoing_option_cost_df <- matrix(0,nrow = 100, ncol = length(ongoing_clients))
ongoing_eos_df <- matrix(0,nrow = 100, ncol = length(ongoing_clients))
count <- 1
for(i in ongoing_clients){
  client <- data_ongoing %>% filter(Client == i)
  for(j in 1:100){
    x <- ongoing_option_cost(client$`Used amount`, client$`Start date`[1], client$`Maturity date`[1], rate_data$`1Y`, euri_sim, ftp_2y, ftp_sim, pool_coef, p_cancel, p_typeone, beta_coef[1], beta_coef[2], mod2_weibull)
    ongoing_option_cost_df[j, count] <- x[[1]]
    ongoing_eos_df[j, count] <- x[[2]]
  }
  #ongoing_option_cost(client$`Used amount`, client$`Start date`[1], client$`Maturity date`[1], rate_data$`1Y`, euri_sim, ftp_data, ftp_sim, pool_coef, p_cancel, p_typeone, beta_coef[1], beta_coef[2], mod2_weibull)
  count <- count + 1
}

NPV_alpha_ongoing <- sort(-as.vector(ongoing_option_cost_df))[0.999 * 100 * length(ongoing_clients)]
summary(as.vector(ongoing_eos_df))
hist(ongoing_option_cost_df)
hist(ongoing_eos_df*NPV_alpha_ongoing)

######################################################################################################