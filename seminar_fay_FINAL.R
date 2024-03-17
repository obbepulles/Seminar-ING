library(readxl)
library(lubridate)
library(dplyr)
library(MASS)
library(fitdistrplus)
library(RColorBrewer)
library(ggplot2)
library(tidyr)
library(timeDate)
set.seed(20189)

excel_file <- "hypothetical_data_set.xlsx"
data <- read_excel(excel_file,sheet = 2,skip = 1)

data$`Reporting date` <- as.Date(data$`Reporting date`, format = "%d-%m-%Y")
data$`Cancellation date` <- as.Date(data$`Cancellation date`, format = "%d-%m-%Y")
data$`Start date` <- as.Date(data$`Start date`, format = "%d-%m-%Y")
data$`Maturity date` <- as.Date(data$`Maturity date`, format = "%d-%m-%Y")

time_difference_in_months <- interval(data$`Start date`[5], 
                                      data$`Cancellation date`[5]) %/% months(1)


###### How many clients cancel their Subscription?
###### What's the utilized time span for each client?

#### Now let's make a table
data_easy <- data %>%
  group_by(Client) %>%
  slice(n()) # Select the last row for each group



data_easy$Cancel <- ifelse(!is.na(data_easy$`Cancellation date`), 1, 0)
# If not canceled, how many months does the subscription last
data_easy$Expected_time <- interval(data_easy$`Start date`, 
                                    data_easy$`Maturity date`) %/% months(1)
data_easy$Matured <- ifelse(data_easy$Cancel == 0 & data_easy$`Maturity date` == data_easy$`Reporting date`, 1, 0)

# If canceled, how many months does the subscription last
data_easy$Time_lasting <- data_easy$Expected_time
data_easy$Time_lasting[data_easy$Matured==0] <- interval(data_easy$`Start date`[data_easy$Matured==0], data_easy$`Reporting date`[data_easy$Matured==0]) %/% months(1)
data_easy$Time_lasting <- as.numeric(data_easy$Time_lasting)
# Calculate the difference between Reporting date and Maturity date
data_easy$Difference <- interval(data_easy$`Reporting date` , data_easy$`Maturity date`) %/% months(1)
# For rows where Matured == 1, set the Difference to NA
data_easy$Difference[data_easy$Matured == 1 | data_easy$Cancel == 1] <- NA

data_easy$Matured[data_easy$Difference == 1] <- 1
data_easy$Ended <- ifelse(data_easy$Cancel== 1 | data_easy$Matured == 1, 1, 0)

# Average utility and variance
u_bar <- c()
v_bar <- c()
i=0
while (i <= 1199){
  subdata <- filter(data, Client==i)
  u <- mean(subdata$`Used amount`)
  v <- sqrt(var(subdata$`Used amount`))
  u_bar <- append(u_bar, u)
  v_bar <- append(v_bar, v)
  i = i+1
}

table_client_based <- data.frame(
  Client = data_easy$Client,
  Canceled = data_easy$Cancel,
  Matured = data_easy$Matured,
  Ended = data_easy$Ended,
  Time_Lasting = data_easy$Time_lasting,
  Expected_maturity = data_easy$Expected_time,
  Used_time = data_easy$Time_lasting/data_easy$Expected_time,
  Average_used_amount = u_bar,
  std_used_amount = v_bar
)


###### How many matured | Should I conclude the almost matured ones???
sum(table_client_based$Matured)

# Data overview
hist(table_client_based$Used_time, probability = TRUE)
k = hist(table_client_based$Used_time[table_client_based$Ended== 1],prob = TRUE)

list_of_tables <- list()
for (q in 1:10) {
  table <- table_client_based[table_client_based$Expected_maturity==12*q,]
  list_of_tables[[q]] <- table
}

x_axis <- seq(1,10,1) ## lose maturity=9&10
y <- c()
for (i in x_axis) {
  yy <- sum(table_client_based$Canceled[table_client_based$Expected_maturity==i*12])/sum(list_of_tables[[i]]$Ended)  
  y <- append(y,yy)
}

y_axis <- y
plot(x_axis, y_axis) ##Looks like an exponential function

#### what distribution of the used_time? Beta seems good
x = table_client_based$Used_time[table_client_based$Canceled==1]
hist(x, freq = FALSE, main = "Histogram with Fitted Distributions")


exp_fit <- fitdist(x, "exp")
beta_fit <- fitdist(x, "beta")
gamma_fit <- fitdist(x, "gamma")
weibull_fit <- fitdist(x, "weibull")
plot(beta_fit, col = "red", add = TRUE)

# Create histogram with true data
hist_data <- ggplot(data = data.frame(x), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "lightblue", color = "black") +
  labs(title = "Histogram with Fitted Distributions", x = "Used Time", y = "Density") +
  scale_color_manual(values = c("red", "purple", "magenta"),
                     labels = c("Beta", "Gamma", "Weibull"))

# Overlay density curves for fitted distributions
hist_with_fitted <- hist_data +
  stat_function(fun = dbeta, args = list(shape1 = beta_fit$estimate[1], shape2 = beta_fit$estimate[2]), aes(color = "Beta")) +
  stat_function(fun = dgamma, args = list(shape = gamma_fit$estimate[1], rate = gamma_fit$estimate[2]), aes(color = "Gamma")) +
  stat_function(fun = dweibull, args = list(shape = weibull_fit$estimate[1], scale = weibull_fit$estimate[2]), aes(color = "Weibull")) +
  theme_minimal()

# Plot the histogram with fitted distributions
print(hist_with_fitted)


# Calculate AIC values
aic_values <- c(exponential = exp_fit$aic,
                beta = beta_fit$aic, 
                gamma = gamma_fit$aic, 
                weibull = weibull_fit$aic)

# Create a data frame to store the results
aic_table <- data.frame(Distribution = names(aic_values),
                        AIC = aic_values)

# Print the table
print(aic_table)

#######################################
# Function to calculate maturity rate and cancellation rate
calculate_rate <- function(maturity, field) {
  sum(table_client_based[[field]][table_client_based$Expected_maturity == maturity * 12]) / 120
}

# Calculate maturity and cancellation rates
x_axis <- seq(1, 10)
maturity_rate <- sapply(x_axis, function(i) calculate_rate(i, "Matured"))
cancellation_rate <- sapply(x_axis, function(i) calculate_rate(i, "Canceled"))
ongoing_rate <- 1 - maturity_rate - cancellation_rate

# Create data frames for ggplot
maturity_df <- data.frame(Maturity = x_axis, Rate = maturity_rate, Type = "Maturity Rate")
cancellation_df <- data.frame(Maturity = x_axis, Rate = cancellation_rate, Type = "Cancellation Rate")
ongoing_df <- data.frame(Maturity = x_axis, Rate = ongoing_rate, Type = "Ongoing Rate")

# Combine data frames
combined_df <- rbind(maturity_df, cancellation_df, ongoing_df)

print(combined_df)
# Create plots with adjusted axis ranges
plots <- ggplot(combined_df, aes(x = Maturity, y = Rate, color = Type)) +
  geom_line() +
  geom_point() +
  labs(x = "Expected Maturity", y = "Rate") +
  theme_minimal()

# Adjust plot dimensions to maintain the desired aspect ratio
plots + theme(
  aspect.ratio = 7/8,
  plot.margin = unit(c(1, 1, 1, 1), "cm")
)



############################
weighted_mature <- c()
for (i in 1:10) {
  weighted_mature[i] <- sum(table_client_based$Matured[table_client_based$Expected_maturity == i * 12]) / 423
}
average_cancellation <- weighted.mean(cancellation_rate,weighted_mature)

theta_values <- data.frame(Index = 1:10, Theta = weighted_mature)

print(theta_values)
##############################
#### To see if the maturity affects the cancellation prob.
model <- glm(Canceled ~ log(Time_Lasting) + Average_used_amount 
             + std_used_amount, data = table_client_based[table_client_based$Expected_maturity < 8*12,], family = binomial)
summary(model)

### Expected_maturity significant
hist(table_client_based$Expected_maturity, breaks = 10, probability = TRUE)
x_axis <- seq(1,10,1) ## lose maturity=9&10
y <- c()
for (i in x_axis) {
  yy <- sum(table_client_based$Canceled[table_client_based$Expected_maturity==i*12])/sum(list_of_tables[[i]]$Ended)  
  y <- append(y,yy)
}

y_axis <- y
plot(x_axis, y_axis) ##Looks like an exponential function


#### Model the cancellation rate
exp_function <- function(x, a, b) {
  a * exp(b * x)
}

fit <- nls(y_axis ~ exp_function(x_axis, a, b), 
           start = list(a = 0.1, b = 0.01))
summary(fit)
curve(exp_function(x, coef(fit)["a"], coef(fit)["b"]), 
      from = 1, 
      to = 10, 
      add = TRUE,
      col = "red")

prob_of_cancellation <- function(maturity_in_month) {
  p <- 0.1 * exp(0.02 * maturity_in_month)
  return(p)
}



### Time_lasting significant
hist(table_client_based$Time_Lasting[table_client_based$Ended==1&table_client_based$Expected_maturity==12], breaks = 10, probability = TRUE)
x_axis <- seq(1,120)
y <- c()
yy <- c()
for (i in x_axis) {
  yy <- sum(table_client_based$Canceled[table_client_based$Time_Lasting == i])/191
  y <- append(y,yy)
}

y_axis <- y
plot(x_axis, y_axis) ##Looks like 

##### Or... for marturity = 24
x_axis <- seq(1,120)
y <- c()
for (i in x_axis) {
  yy <- sum(table_client_based$Canceled[table_client_based$Expected_maturity==120&table_client_based$Used_time==i/120])/sum(table_client_based$Canceled[table_client_based$Expected_maturity==120])
  y <- append(y,yy)
}

y_axis <- y
plot(x_axis, y_axis) ##### Not significant

####### Plot used_time and maturity

x_axis <- seq(1,10)
y <- c()
for (i in x_axis) {
  yy <- mean(table_client_based$Used_time[table_client_based$Expected_maturity==i*12 & table_client_based$Canceled==1])
  y <- append(y,yy)
}

y_axis <- y
plot(x_axis, y_axis)

##### Simulate the \tau for ongoing client
calculate_tau <- function(T_in_month){
  random_c <- runif(1)
  if (random_c >= average_cancellation){
    tau <- T_in_month
  }
  else tau <- round(rbeta(1, shape1 = beta_fit$estimate[1], shape2 = beta_fit$estimate[2])*T_in_month)
  
  if (tau == 0){
    tau = 1
  }
  return(tau)
}

table_ongoing <- table_client_based[table_client_based$Ended==0,1:6]
table_ongoing$time2marturity <- table_ongoing$Expected_maturity - table_ongoing$Time_Lasting

next_values <- numeric(nrow(table_ongoing))

# Iterate over each row of table_ongoing and calculate next_ value
for (i in 1:nrow(table_ongoing)) {
  time2marturity <- table_ongoing$Expected_maturity[i] - table_ongoing$Time_Lasting[i]
  next_values[i] <- calculate_tau(time2marturity)
}

# Add the calculated values to table_ongoing
table_ongoing$next_ <- next_values


##### Transform it in original format
data_1200_completed <- data[, c(2:4, 6:8)]

# Iterate through each row of table_ongoing
for (i in 1:nrow(table_ongoing)) {
  # Get the client number
  client_number <- table_ongoing$Client[i]
  print(client_number)
  # Find the corresponding rows in data_1200_completed
  client_rows <- which(data_1200_completed$Client == client_number)
  
  # Determine the number of rows to add based on table_ongoing$next_
  num_new_rows <- table_ongoing$next_[i]
  
  # Get the last reporting date for the client
  start_date = tail(data_1200_completed$`Start date`[client_rows], n = 1)
  maturity_date = tail(data_1200_completed$`Maturity date`[client_rows], n = 1)
  last_reporting_date <- tail(data_1200_completed$`Reporting date`[client_rows], n = 1)
  
  # Generate new reporting dates based on the last reporting date
  new_reporting_dates <- seq(as.Date(last_reporting_date)+1, length=num_new_rows+1, by="months") - 1
  #new_reporting_dates <- as.Date(format(new_reporting_dates, "%Y-%m-%d")) - 1
  print(new_reporting_dates)
  # Create new rows and append to data_1200_completed
  new_rows <- data.frame(`Reporting date` = tail(new_reporting_dates,num_new_rows),
                         `Cancellation date` = rep(NA, num_new_rows),
                         Client = rep(client_number, num_new_rows),
                         `Start date` = rep(start_date, num_new_rows),
                         `Maturity date` = maturity_date,
                         `Used amount` = rep(NA, num_new_rows),
                         stringsAsFactors = FALSE)
  # Rename the columns of new_rows
  names(new_rows) <- c("Reporting date", "Cancellation date", "Client", "Start date", "Maturity date", "Used amount")
  
  # Now, new_rows should have the same column names as data_1200_completed
  
  data_1200_completed <- rbind(data_1200_completed, new_rows)
}

# Make sure to reorder rows by client and reporting date
data_1200_completed <- data_1200_completed[order(data_1200_completed$Client, 
                                                 data_1200_completed$`Reporting date`), ]


######## Fix the cancellation date
for (i in 1:nrow(table_ongoing)) {
  client_number <- table_ongoing$Client[i]
  # Find the corresponding rows in data_1200_completed
  client_rows <- which(data_1200_completed$Client == client_number)
  maturity_date = tail(data_1200_completed$`Maturity date`[client_rows], n = 1)
  last_reporting_date <- tail(data_1200_completed$`Reporting date`[client_rows], n = 1)
  if (maturity_date != last_reporting_date){
    data_1200_completed$`Cancellation date`[data_1200_completed$Client==client_number] <- last_reporting_date
  }
}
