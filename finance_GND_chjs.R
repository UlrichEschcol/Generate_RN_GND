# Load required libraries
library(quantmod)
library(ggplot2)

# ---------- Step 1: Download and Prepare Data ----------

# Helper function to calculate log returns
calculate_log_returns <- function(prices) {
  return(c(NA, diff(log(prices))))
}

# Function to download price data and compute returns
get_stock_returns <- function(tickers, start_date, end_date) {
  stock_returns <- list()
  for (ticker in tickers) {
    getSymbols(ticker, src = "yahoo", from = start_date, to = end_date, auto.assign = TRUE)
    price_data <- Cl(get(ticker))
    log_ret <- na.omit(calculate_log_returns(price_data))
    stock_returns[[ticker]] <- as.numeric(log_ret)
  }
  return(stock_returns)
}

# ---------- Step 2: Density Functions ----------

# GND density function
# GND density function
dgamma_gnd <- function(x, mu, sigma, gamma1) {
  C <- gamma(1/2 + 1) / gamma((gamma1 - 1) / gamma1 + 1) * ((gamma1 - 1) / gamma1)^((gamma1 - 1) / gamma1)
  return((C / (sqrt(pi) * sigma)) * exp(-((gamma1 - 1) / gamma1) * (abs(x - mu) / sigma)^(gamma1 / (gamma1 - 1))))
}

# Estimated GND standard deviation from data
estimate_sigma_gnd <- function(data, mu_est, gamma1) {
  n <- length(data)
  exponent <- gamma1 / (gamma1 - 1)
  power <- 2 * (gamma1 - 1) / gamma1
  sigma2_est <- ((1 / n) * sum(abs(data - mu_est)^exponent))^power
  return(sqrt(sigma2_est))  # Return standard deviation
}

# ---------- Step 3: Plotting Function ----------

plot_return_density <- function(returns, ticker, gamma1) {
  returns <- na.omit(returns)
  mu <- mean(returns)
  sigma_gnd <- estimate_sigma_gnd(data = returns, mu_est = mu, gamma1 = gamma1)
  sigma_normal <- sd(returns)
  x_vals <- seq(min(returns), max(returns), length.out = 1000)
  
  gnd_density <- dgamma_gnd(x_vals, mu, sigma_gnd, gamma1)
  normal_density <- dnorm(x_vals, mean = mu, sd = sigma_normal)
  
  # Create dynamic GN label with ?? value
  gn_label <- paste0("GN (\u03B3 = ", gamma1, ")")
  color_map <- setNames(c("blue", "red", "black"), c("Real Data", gn_label, "Normal"))
  ggplot(data.frame(Return = returns), aes(x = Return)) +
    geom_line(stat = "density", aes(y = ..density.., color = "Real Data"), size = 1) +
    geom_line(data = data.frame(x = x_vals, y = gnd_density),
              aes(x = x, y = y, color = gn_label), size = 1) +
    geom_line(data = data.frame(x = x_vals, y = normal_density),
              aes(x = x, y = y, color = "Normal"), size = 1, linetype = "dashed") +
    scale_color_manual(values = color_map) +
    labs(#title = paste0(ticker, " Returns Density"),
         x = "Daily Returns", y = "Density") +
    theme_minimal() +
    theme(legend.title = element_blank())
}

# ---------- Step 4: Example Usage ----------

# Define tickers and date range
tickers <- c("7203.T","SPM.MI","TRN.MI", "ELET3.SA")  # or c("7203.T","SPM.MI","TRN.MI", "ELET3.SA")
returns_list <- get_stock_returns(tickers, start_date = "2022-01-01", end_date = "2025-05-31")

#returns = returns_list[[tick]]
#mean(returns)
#sd(returns)
#estimate_sigma_gnd(returns,mu_est = mean(returns), gamma1 = 8 )


# Plot ELET3.SA with gamma = 2.6
plot_return_density(returns_list[["ELET3.SA"]], "ELET3.SA", gamma1 = 2.6)

# Plot TRN.MI with gamma = 3.8
plot_return_density(returns_list[["TRN.MI"]], "TRN.MI", gamma1 = 3.8)

# Plot 7203.T with gamma = 4.7
plot_return_density(returns_list[["7203.T"]], "7203.T", gamma1 = 4.7)

# Plot SPM.MI with gamma = 8.0
plot_return_density(returns_list[["SPM.MI"]], "SPM.MI", gamma1 = 8.0)
