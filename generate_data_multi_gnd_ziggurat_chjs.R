library(MASS)
library(matrixcalc)  # For matrix determinant and inverse
library(ggplot2)
library(plotly)
library(viridis)  # Explicitly load viridis for color scales

setwd('C:/Users/... /RScript') # Set the path to the file location
source('functions_gnd_ziggurat_chjs.R')

mvgnd <- function(n = 1, mu, Sigma, gamma1, tol = 1e-6) {
  p <- length(mu)
  if (!all(dim(Sigma) == c(p, p))) stop("incompatible arguments")
  
  # Eigen decomposition to get square root of Sigma
  eS <- eigen(Sigma, symmetric = TRUE)
  ev <- eS$values
  if (!all(ev >= -tol * abs(ev[1L]))) stop("'Sigma' is not positive semi-definite")
  Sigma_sqrt <- eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(eS$vectors)
  
  # Generate Z ~ N^p_gamma(0, I_p)
  Z <- matrix(NA, nrow = n, ncol = p)
  for (j in 1:p) {
    Z[, j] <- myziggurat_gnd(total = n, N = 256, mu = 0, sigma2 = 1, gamma1 = gamma1)
  }
  
  # Apply transformation: X = mu + Sigma^{1/2} * Z
  X <- t(apply(Z, 1, function(z) drop(mu + Sigma_sqrt %*% z)))
  
  # Set names
  nm <- names(mu)
  if (is.null(nm) && !is.null(dn <- dimnames(Sigma))) nm <- dn[[1L]]
  dimnames(X) <- list(NULL, nm)
  
  return(X)
}

##############
##############

# Define parameters
mu <- c(0,0)  # Mean vector
Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Covariance matrix
Sigma <- matrix(c(1, 0.1, 0.1, 1), 2, 2)  # Covariance matrix
Sigma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)  # Covariance matrix
Sigma <- matrix(c(1, 0.5, 0.5, 1), 2, 2)  # Covariance matrix

gamma1 <- 4  # Shape parameter
# Generate multivariate GND samples
samples <- mvgnd(n = 1000, mu = mu, Sigma = Sigma, gamma1=gamma1)

#head(samples)
#var(samples)
#cor(samples)
colnames(samples) = c("X","Y")
samples = as.data.frame(samples)
samples$Z = apply(samples, 1, function(x) mvgnd_pdf(x, mu, Sigma, gamma1))

##############
##############
# Convert samples to dataframe
df <- as.data.frame(samples)
colnames(df) <- c("X", "Y", "Density")

# Create interactive 3D surface plot
plot_ly(df, x = ~X, y = ~Y, z = ~Density, type = "mesh3d",
        intensity = ~Density, colorscale = "Viridis") %>%
  layout(#title = "3D Surface Plot of Generalized Normal Distribution",
         scene = list(xaxis = list(title = "X"),
                      yaxis = list(title = "Y"),
                      zaxis = list(title = "Density")))

##############
##############

df <- as.data.frame(samples)
colnames(df) <- c("X", "Y", "Z")

# Create interactive 3D scatter plot
plot_ly(df, x = ~X, y = ~Y, z = ~Z, type = "scatter3d", mode = "markers",
        marker = list(size = 3, color = df$Z, colorscale = "Viridis")) %>%
  layout(title = "3D Scatter Plot of Bivariate GND Samples",
         scene = list(xaxis = list(title = "X"),
                      yaxis = list(title = "Y"),
                      zaxis = list(title = "Z")))

##############
##############
##  p = 3

# Define parameters
mu <- c(0,0,0)  # Mean vector
Sigma <- matrix(c(1, 0.5, 0.25, 0.5, 1, 0.35, 0.25, 0.35, 1), 3, 3)  # Covariance matrix


gamma1 <- 1.01  # Shape parameter
# Generate multivariate GND samples
samples <- mvgnd(n = 1000, mu = mu, Sigma = Sigma, gamma1=gamma1)

##
colnames(samples) = c("X","Y","Z")
samples = as.data.frame(samples)
samples$Density = apply(samples, 1, function(x) mvgnd_pdf(x, mu, Sigma, gamma1))
##
