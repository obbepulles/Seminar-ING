library('dplyr')
library('censReg')
library('moments')
library('ggplot2')
library('ggpubr')
library("mixR")
library('forecast')

data <- readxl::read_xlsx("hypothetical_data_set.xlsx",2, skip = 1)
data <- data[,2:9]
data$`Reporting date` <- as.Date(data$`Reporting date`, "%d-%m-%Y")
data$`Maturity date` <- as.Date(data$`Maturity date`, "%d-%m-%Y")
data$`Cancellation date` <- as.Date(data$`Cancellation date`, "%d-%m-%Y")
rate_data <- readxl::read_xlsx("hypothetical_data_set.xlsx",1, skip = 1, range = "H2:M122")
ftp_data <- readxl::read_xlsx("hypothetical_data_set.xlsx",1, skip = 1, range = "A2:F122")


data <- data %>% mutate(month_year = format(`Reporting date`, "%m-%Y")) 
rate_data <- rate_data %>% mutate(month_year = format(`Date`, "%m-%Y")) 
ftp_data <- ftp_data %>% mutate(month_year = format(`Date`, "%m-%Y"))

merged_df <- inner_join(data, rate_data, by = "month_year")
data <- merged_df[,-9]

#calculated dU_t
lagged_use <- data %>% 
    group_by(Client) %>%
    mutate(dU = `Used amount` - lag(`Used amount`))

months_used <- tally(data%>%group_by(Client))$n
maturities <- data %>% group_by(Client) %>% slice(1)
plot(hist(maturities$Maturity))

hist(months_used)
hist(lagged_use$dU)
hist(data$`Used amount`, freq = FALSE)

#Initial use seems to follow uniform distr
first_utilization <- data %>%
  group_by(Client) %>% summarize(first(`Used amount`))


  
hist(first_utilization$`first(\`Used amount\`)`, main = "Histogram of first utilization in a contract", col = 'lightblue', xlab = "Used amount as fraction")

#Filter out constant U_t's
filtered_df <- data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 1)) 

lagged_use_filtered <- filtered_df %>% 
  group_by(Client) %>%
  mutate(dU = `Used amount` - lag(`Used amount`))
hist((lagged_use_filtered$`dU`), freq = FALSE)
plot(density(na.omit(lagged_use_filtered$dU)))
qqnorm(lagged_use_filtered$dU)
qqline(lagged_use_filtered$dU)
kurtosis(na.omit(lagged_use_filtered$dU))

#How many clients use the credit line more than once?
num_clients_filter <- filtered_df %>% summarise(NumClients = n_distinct(Client))
num_clients <- data %>% summarise(numClients = n_distinct(Client))
#num_clients - num_clients_filter
#Single observation clients
num_single <- data %>% group_by(Client) %>% filter(n() == 1) %>% summarise(numClients = n_distinct(Client)) %>% count(numClients)
#Estimated percentage of credit line users
#p <- as.numeric((num_clients - num_clients_filter + num_single$n) / (num_clients - num_single$n))

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

pcaf_pvals <- rep(0,11)
for(i in 1:11){
  pacf_results <- filtered_df %>%
    group_by(Client) %>%
    summarize(PACF = list(pacf(`Used amount`, plot = FALSE)$acf[i]))
  corls <- as.numeric(pacf_results$PACF)
  pcaf_pvals[i] <- t.test(x = corls, alternative = "two.sided")$p.value
  #hist(corls)
  #Sys.sleep(4)
}

# arma_results <- filtered_df %>% group_by(Client) %>%
#                   summarize(ARMA = list(arimaorder(auto.arima(`Used amount`, ic = 'aic'))))
# 
# arimasum <- c(0,0,0)
# for(i in 1:length(arma_results$Client)){
#   arimasum <- arimasum + (as.numeric(arma_results$ARMA[i][[1]]))
# }
# arimasum <- arimasum / length(arma_results$Client)



#(1) We see two types of clients:
#       type 1: Use credit line once, don't touch again (~20% of clients)
#       type 2: Use credit line more often, high autocorrelation (Tobit model?, ~80% of clients)
#(2) Initial amount used seems to be Unif(0,1)
#(3) No ACF in dU
#(4) For type 2 clients, dU seems to follow a normal distr around (-0.2;0.2) but has thick tails
#    Ideas: Mixture of normals / Truncated normal + other distr for tails
#    Problem: U in [0,1], thin tails make sense, could be a result of TOBIT

#DGP idea: 
#(1) U_0 ~ Unif[0,1]
#(2) eps_t ~ Norm(mu,sigma^2)
#(3) split up U_t into cases: U_t = 
#     (i)    1              , if phi * U_{t-1} + eps_t >= 1
#     (ii)   0              , if phi * U_{t-1} + eps_t <= 0
#     (iii)  phi * U_{t-1} + eps_t, otherwise

filtered_df <- data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 1)) %>%
  ungroup()
lagged_use_filtered <- filtered_df %>% 
  group_by(Client) %>%
  mutate(dU = `Used amount` - lag(`Used amount`))

# U_t = phi*U_{t-1} + epsilon_t

lag_df <- lagged_use_filtered %>% group_by(Client) %>% mutate(`Lag used` = lag(`Used amount`))
Y <- lag_df$`Used amount`[1:38567]
X <- lag_df$`Lag used`[1:38567]
modeldata <- na.omit(as.data.frame(cbind(Y,X)))


#ARMA?
#a <- auto.arima(modeldata$X, ic = "aic")
#arima(Y, order = c(1,0,0))

tobit_model <- censReg(Y ~ X, left = 0, right = 1, data = modeldata)
summary(tobit_model)
error <- modeldata$Y - modeldata$X * tobit_model$estimate[2] - tobit_model$estimate[1]

plot(density(error^2))
qqnorm(y = error)
qqline(y = error, col = 2)
###########################################

#Some of the type 2 clients seem more active in their use
#Plot of means does not give anything special, SD plot does (mixture of weibull/other)
filtered_df <- filtered_df %>% 
  group_by(Client) %>%
  mutate(dU = `Used amount` - lag(`Used amount`))

filtered_df <- filtered_df %>% 
  filter(!is.na(`dU`))

c3 <- filtered_df %>% group_by(Client) %>% mutate(sds = sd(`Used amount`)) %>% slice(1)
c4 <- filtered_df %>% group_by(Client) %>% 
        filter(!is.na(dU)) %>% filter(n() > 1) %>% 
          mutate(sds = sd(dU)) %>% slice(1)
means <- c4
hist(means$sds^2)
plot(density(means$sds^2))

modes <- function(d){
  i <- which(diff(sign(diff(d$y))) < 0) + 1
  data.frame(x = d$x[i], y = d$y[i])
}

abline(v =modes(density(means$sds^2))$x)
abline(v = 0)

#Plot some clients utilization over time
# for(i in 0:(2000/5)){
# dataSubset <- lagged_use %>% filter(Client >= 345 & Client < 350)
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

#Model the variance with mixR mixed weibulls, performs best in all 3 cirteria
x <- means$sds^2
mod1_weibull <- mixfit(x, ncomp = 2, family = 'weibull')
mod1_weibull
plot(mod1_weibull)
#Within type 2 clients, we make further distinctions:
# Type 1: Low average variance 
# Type 2: High average variance
# IDEA: model this as state models (state 1 = low var, state 2 = high var)


#########################################################################
#Filter those clients whose contracts have already matured
#Last reporting date is 31-12-2021

#Clients whose contracts are either finished or cancelled are included
data_fin_can <- data %>% 
                group_by(Client) %>% 
                filter((last(`Reporting date`) != as.Date('31-12-2021', "%d-%m-%Y")) | (last(`Reporting date`) == last(`Maturity date`)))
#add change in used amount 
data_fin_can <- data_fin_can %>% 
  group_by(Client) %>%
  mutate(dU = `Used amount` - lag(`Used amount`))

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

c5 <- data_nonconst %>% group_by(Client) %>% mutate(vars = var(`Used amount`)) %>% slice(1)
plot(density(c5$vars))

#Gamma performs better on 3 comp, weibull on 2 comp, analysis done on VARIANCE
x <- (c5$vars)
mod1_weibull <- mixfit(x, ncomp = 2, family = 'weibull')
mod1_weibull
plot(mod1_weibull)


lag_nonconst <- data_nonconst %>% group_by(Client) %>% mutate(`Lag used` = lag(`Used amount`))
Y <- lag_nonconst$`Used amount`
X <- lag_nonconst$`Lag used`
modeldata <- (as.data.frame(cbind(Y,X)))

tobit_model <- censReg(X ~ Y, left = 0, right = 1, data = modeldata)
summary(tobit_model)
error <- modeldata$X - modeldata$Y * tobit_model$estimate[2] - tobit_model$estimate[1]
lag_nonconst <- cbind(lag_nonconst, error)

qqnorm(y = error)
qqline(y = error, col = 2)

#Error analysis of tobit model
c6 <- (lag_nonconst %>% group_by(Client) %>% mutate(vars = var(`Used amount`)) %>% slice(1))
plot(hist(c6$vars))

x <- (c6$vars)
mod2_weibull <- mixfit(x, ncomp = 3, family = 'gamma', max_iter = 10000)
mod2_weibull
plot(mod2_weibull)


#analyse if covid had an effect
covid_date <- as.Date("01-3-2020", "%d-%m-%Y")
lag_nonconst_covid <- lag_nonconst %>% mutate(covid = if_else(`Reporting date` >= covid_date,1,0))

#We observe no interaction between utiization and euribor / covid
tobit_model_covid <- censReg(data = lag_nonconst_covid, `Lag used` ~ `Used amount` + covid + `6M` + `1Y` + `3Y` + `5Y` + `10Y`)
summary(tobit_model_covid)


########################################################################
#calculate probability of cancelling early
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
########################################################################
#Calculate probability of being type 1 (constant utilization)

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
########################################################################

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
 util <- util$`Used amount` + 1
 for(j in 1:18){
   
   tryCatch(
     {
       arima_model[j] <- AIC(arima(util, order = param_grid[j, ]))
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


