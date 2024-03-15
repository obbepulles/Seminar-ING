#Filter out constant U_t's
filtered_df <- data %>%
  group_by(Client) %>%
  filter((length(unique(`Used amount`)) > 1) & (length(`Used amount`) > 23)) %>%
  ungroup()

lagged_use_filtered <- filtered_df %>% 
  group_by(Client) %>%
  mutate(dU = `Used amount` - lag(`Used amount`))

#--- NOT RELEVANT ANYMORE AS WE KNOW AR(1) IS BETTER THAN ARIMA(0,1,0)
# hist((lagged_use_filtered$`dU`), freq = FALSE)
# summary(lagged_use_filtered$dU)
# plot(density(na.omit(lagged_use_filtered$dU)))
# qqnorm(lagged_use_filtered$dU)
# qqline(lagged_use_filtered$dU)
# kurtosis(na.omit(lagged_use_filtered$dU))

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
random_t <- mean(autoCorrdU$Autocorrelation) + sd(autoCorrdU$Autocorrelation)*rt(800,df = length(autoCorrdU$Autocorrelation))
qqnorm(random_t)

pcaf_pvals <- rep(0,11)
for(i in 1:11){
  pacf_results <- filtered_df %>%
    group_by(Client) %>%
    summarize(PACF = list(pacf(`Used amount`, plot = FALSE)$acf[i]))
  corls <- as.numeric(pacf_results$PACF)
  pcaf_pvals[i] <- t.test(x = corls, alternative = "two.sided")$p.value
  hist(corls)
Sys.sleep(4)
}
