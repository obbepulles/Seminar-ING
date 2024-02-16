library(readxl)
library(lubridate)
library(dplyr)
library(MASS)
library(fitdistrplus)
library(RColorBrewer)


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
  Time_Lasting = data_easy$Time_lasting,
  Expected_maturity = data_easy$Expected_time,
  Used_time = data_easy$Time_lasting/data_easy$Expected_time,
  Average_used_amount = u_bar,
  std_used_amount = v_bar
)


###### How many matured | Should I conclude the almost matured ones???
sum(table_client_based$Matured)


# Calculate the difference between Reporting date and Maturity date
data_easy$Difference <- interval(data_easy$`Reporting date` , data_easy$`Maturity date`) %/% months(1)
# For rows where Matured == 1, set the Difference to NA
data_easy$Difference[data_easy$Matured == 1 | data_easy$Cancel == 1] <- NA

# Data overview
hist(table_client_based$Used_time, probability = TRUE)
k = hist(table_client_based$Used_time[table_client_based$Canceled== 1],prob = TRUE)
k



#### what distribution of the used_time? Beta seems good
x = table_client_based$Used_time[table_client_based$Canceled==1]
hist(x, freq = FALSE, main = "Histogram with Fitted Distributions")
plot(beta_fit, col = "red", add = TRUE)

# Create histogram with true data
hist_data <- ggplot(data = data.frame(x), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 10, fill = "lightblue", color = "black") +
  labs(title = "Histogram with Fitted Distributions", x = "Used Time", y = "Density")

# Overlay density curves for fitted distributions
hist_with_fitted <- hist_data +
  stat_function(fun = dbeta, args = list(shape1 = beta_fit$estimate[1], shape2 = beta_fit$estimate[2]), color = "red") +
  theme_minimal()

# Plot the histogram with fitted distributions
print(hist_with_fitted)

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








##############################
#### To see if the maturity affects the cancellation prob.
model <- glm(Canceled ~ Expected_maturity + Time_Lasting + Average_used_amount 
             + std_used_amount, data = table_client_based[table_client_based$Expected_maturity < 8*12,], family = binomial)
summary(model)



### Expected_maturity significant
hist(table_client_based$Expected_maturity, breaks = 10, probability = TRUE)
x_axis <- seq(12,120,12) ## lose maturity=9&10
y <- c()
for (i in x_axis) {
  yy <- sum(table_client_based$Canceled[table_client_based$Expected_maturity==i])/120  
  y <- append(y,yy)
}

y_axis <- y
plot(x_axis, y_axis) ##Looks like an exponential function
hist(y,breaks=20)

#### Model the cancellation rate
exp_function <- function(x, a, b) {
  a * exp(b * x)
}

fit <- nls(y_axis ~ exp_function(x_axis, a, b), 
           start = list(a = 0.1, b = 0.01))
summary(fit)
curve(exp_function(x, coef(fit)["a"], coef(fit)["b"]), 
      from = 12, 
      to = 120, 
      add = TRUE,
      col = "red")

prob_of_cancellation <- function(maturity_in_month) {
  p <- 0.1 * exp(0.02 * maturity_in_month)
  return(p)
}


### Time_lasting significant
hist(table_client_based$Time_Lasting, breaks = 10, probability = TRUE)

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

















