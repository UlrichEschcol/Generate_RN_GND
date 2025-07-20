library(stats)
library(dplyr)
library(cubature)
library(mvtnorm)

setwd('C:/Users/... /RScript') # Set the path to the file location
source('functions_gnd_ziggurat_chjs.R')


# Compute the volume under the gamma-order generalized normal PDF in dimension p > 1

# Define the p-dimensional gamma-order generalized normal distribution PDF
dgamma_gnd_p <- function(x, mu, Sigma, gamma) {
  p <- length(mu)  # Dimension of the data
  
  # Compute determinant and inverse of Sigma
  det_Sigma <- det(Sigma)
  inv_Sigma <- solve(Sigma)
  
  # Compute the p-quadratic form Q_??(x)
  Q_theta_x <- t(x - mu) %*% inv_Sigma %*% (x - mu)
  
  # Compute the normalizing constant C^p_??
  C_p_gamma <- pi^(-p/2) * gamma(p/2 + 1) / gamma(p*(gamma-1)/gamma + 1) * ((gamma-1)/gamma)^(p*((gamma-1)/gamma))
  
  # Compute the density
  pdf_value <- C_p_gamma * det_Sigma^(-1/2) * 
    exp(-((gamma - 1) / gamma) * (Q_theta_x)^(gamma / (2*(gamma - 1))))
  
  return(as.numeric(pdf_value))
}

## The multivariate CDF

cdf_gnd_p <- function(X, mu, Sigma, gamma1, tol = 1e-11) {
  p <- length(mu)
  
  # Wrapper for multivariate PDF
  pdf_wrapper <- function(x) {
    if (is.null(dim(x))) x <- matrix(x, nrow = 1)
    apply(x, 1, function(row) dgamma_gnd_p(row, mu, Sigma, gamma1))
  }
  
  # Perform integration
  result <- adaptIntegrate(pdf_wrapper,
                           lowerLimit = rep(-Inf, p),
                           upperLimit = X,
                           tol = tol)
  
  # Return as named list
  return(list(
    Volume = as.numeric(result$integral),
    Error = as.numeric(result$error)
  ))
}


###### Test
# Set parameters
p <- 2
mu <- rep(0, p)
Sigma <- diag(1, p)  # Identity covariance
gamma1 <- 2

Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Covariance matrix

cdf_gnd_p(X = c(0,0), mu = rep(0, p), Sigma = Sigma, gamma1)


n=2
mean <- rep(0, n)
lower <- rep(-Inf, n)
upper <- rep(0, n)
corr <- Sigma
pmvnorm(mean, corr, lower, upper)

############################
## Compute table of CDF
############################

# Parameters
gamma_values <- c(-50, -4, -2.2, -2,-1.8, -1.01, 1.01, 1.8,2,2.2, 4, 50)

gamma_values <- c(-2.2,-1.8, 1.8,2.2)

grid_vals <- c(-3, -2, -1, 0, 1, 2, 3)

#gamma_values <- c(-50, -4)
#grid_vals <- c(-3, -2)

# Mean and covariance
p <- 2
mu <- rep(0, p)
Sigma <- diag(1, p)

Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Covariance matrix


# Points of the form (x, x, x)
diag_points <- lapply(grid_vals, function(v) c(v, v))

#diag_points <- lapply(grid_vals, function(v) c(v, v, v))

# Initialize table_result from scratch to ensure no NA overwriting
table_result <- data.frame(gamma = gamma_values)

# Evaluate CDF for each (x, x)
for (pt in diag_points) {
  pt_label <- paste0("F(x,y <= ", pt[1], ")")
  
  cdf_vals <- numeric(length(gamma_values))  # Pre-allocate numeric vector
  
  for (i in seq_along(gamma_values)) {
    gamma1 <- gamma_values[i]
    
    cat("Evaluating CDF at (x, y) =", paste(pt, collapse = ", "), "with gamma1 =", gamma1, "\n")
    cdf_vals[i] <- cdf_gnd_p(pt, mu, Sigma, gamma1)$Volume
    print(cdf_vals[i])
  }
  
  table_result[[pt_label]] <- cdf_vals
}

print(table_result)

######  UNIVARIATE CASE p = 1

#generalized_normal_cdf (x=0, mu=0, sigma2=1, gamma1=2) 
  
# Gamma values to evaluate
gamma_values <- c(-50, -4, -2.2, -2, -1.8, -1.01, 1.01, 1.8, 2, 2.2, 4, 50)

# Univariate evaluation points
grid_vals <- c(-3, -2, -1, 0, 1, 2, 3)

# Parameters
mu <- 0
sigma2 <- 1

# Initialize results table
table_result <- data.frame(gamma = gamma_values)

# Loop over evaluation points
for (x_val in grid_vals) {
  pt_label <- paste0("F(x ??? ", x_val, ")")
  cdf_vals <- numeric(length(gamma_values))  # Pre-allocate
  
  for (i in seq_along(gamma_values)) {
    gamma1 <- gamma_values[i]
    
    cat("Evaluating CDF at x =", x_val, "with gamma1 =", gamma1, "\n")
    
    cdf_val <- tryCatch({
      generalized_normal_cdf(x_val, mu = mu, sigma2 = sigma2, gamma1 = gamma1)
    }, error = function(e) NA)
    
    cdf_vals[i] <- round(cdf_val,4)
  }
  
  table_result[[pt_label]] <- cdf_vals
}

# View final table
print(table_result)



##################
## Check the values between cdf_gnd_p and pmvnorm
##################

# Parameters
n <- 2
mu <- rep(0, n)
Sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)  # Example covariance matrix
gamma1 <- 2  # Normal case for GND

# Define test points (e.g., (x, x))
grid_vals <- c(-3,-2, -1, 0, 1, 2,3)
X_points <- lapply(grid_vals, function(x) rep(x, n))  # X = (x, x)

# Initialize result table
result_table <- data.frame(
  X = sapply(X_points, function(x) paste0("(", paste(x, collapse = ","), ")")),
  cdf_gnd_p = NA,
  pmvnorm = NA,
  difference = NA
)

# Evaluate at each point
for (i in seq_along(X_points)) {
  x <- X_points[[i]]
  
  # Compute GND CDF with gamma1 = 2
  cdf_gnd_val <- cdf_gnd_p(X = x, mu = mu, Sigma = Sigma, gamma1 = gamma1)$Volume
  
  # Compute classical multivariate normal CDF
  pmvnorm_val <- pmvnorm(lower = rep(-Inf, n), upper = x, mean = mu, sigma = Sigma)
  
  # Fill table
  result_table$cdf_gnd_p[i] <- cdf_gnd_val
  result_table$pmvnorm[i] <- pmvnorm_val
  result_table$difference[i] <- abs(cdf_gnd_val - pmvnorm_val)
}

# Print result
print(result_table)
