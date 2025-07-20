############################################################################
##### Functions to Generate Data Following GND with Ziggurat Algorithm #####
############################################################################

# Define the PDF for the GND
dgamma_gnd <- function(x, mu, sigma2, gamma1) {
  lambda_gamma = gamma(1/2 + 1) / gamma((gamma1 - 1) / gamma1 + 1) * ((gamma1 - 1) / gamma1)^((gamma1 - 1) / gamma1)
  pdf_values <- (lambda_gamma / sqrt(pi * sigma2)) * exp(-((gamma1 - 1) / gamma1) * (abs(x - mu) / sqrt(sigma2))^(gamma1 / (gamma1 - 1)))
  return(pdf_values)
}

# Compute the inverse GND PDF
gnd_inverse <- function(y, mu, sigma2, gamma1) {
  # Compute lambda_gamma inside
  lambda_gamma <- {
    num <- gamma(1/2 + 1)
    den <- gamma((gamma1 - 1) / gamma1 + 1)
    power_term <- ((gamma1 - 1) / gamma1)^((gamma1 - 1) / gamma1)
    (num / den) * power_term
  }
  
  # Compute x
  sigma <- sqrt(sigma2)
  inside_log <- (y * sqrt(pi * sigma2)) / lambda_gamma
  term <- (-gamma1 / (gamma1 - 1)) * log(inside_log)
  adjustment <- term^((gamma1 - 1) / gamma1)
  
  x_pos <- mu + sigma * adjustment
  x_neg <- mu - sigma * adjustment
  
  # Determine which value to return based on the sign of x - mu
  if (x_pos - mu >= 0) {
    return(x_pos)
  } else {
    return(x_neg)
  }
}

# Define the CDF for Generalized Normal Distribution (GND)
generalized_normal_cdf <- function(x, mu, sigma2, gamma1) {
  
  sign_x = sign(x)
  
  if (x < 0 ) { x=-x } 
  
  gamma_val <- (gamma1 - 1) / gamma1
  
  term1 <- 1/(2*gamma(gamma_val)) 
  
  term2 <- pgamma(gamma_val*((x-mu)/sqrt(sigma2))^(1/gamma_val), shape = gamma_val, lower.tail = FALSE) * gamma(gamma_val)
  
  cdf_value <- 1 - term1 * term2 
  
  if (sign_x == - 1 ) { cdf_value = 1 - cdf_value } 
  
  return(cdf_value)
}

# Generalized error function component (Erf_{gamma1/(gamma1−1)})
gamma_inc_erf <- function(a, x) {
  integrand <- function(t) exp(-t^a)
  gamma_inc <- integrate(integrand, lower = 0, upper = x)$value
  result <- gamma_inc * gamma(a + 1) / sqrt(pi)
  return(result)
}
###
generalized_normal_cdf1 <- function(z, gamma1) {
  if (gamma1 == 1) {
    return(pnorm(z))  # Standard normal case
  }
  
  # Sign function
  sgn <- function(x) ifelse(x > 0, 1, ifelse(x < 0, -1, 0))
  
  a <- gamma1 / (gamma1 - 1)
  scale <- ((gamma1 - 1) / gamma1)^((gamma1 - 1) / gamma1) * abs(z)
  erf_term <- gamma_inc_erf(a, scale)
  
  cdf <- 0.5 + (sqrt(pi) * sgn(z)) / (2 * gamma((gamma1-1) / gamma1) * gamma(gamma1 / (gamma1 - 1))) * erf_term
  return(cdf)
}

#####


# Modify the Ziggurat algorithm to use the GND CDF
# Adapted from: https://github.com/BoxiLin/ZigguratAlgorithm
findAxy_gnd <- function(N, mu, sigma2, gamma1) {
  x1 <- numeric(N)
  y1 <- numeric(N)
  flag1 = 0
  A1 = 0
  xcap = 10
  xflo = 0
  xcur = (xcap + xflo) / 2
  flag1 = findyN_gnd(N, xcur, mu, sigma2, gamma1)[[1]]
  
  while (flag1 != 1) {
    if (flag1 == 100) {
      xflo = xcur
    } else if (flag1 == -100) {
      xcap = xcur
    }
    xcur = (xcap + xflo) / 2
    flag1 = findyN_gnd(N, xcur, mu, sigma2, gamma1)[[1]]

  }
  A1 = round(findyN_gnd(N, xcur, mu, sigma2, gamma1)[[2]], digits = 6)
  x1 = round(findyN_gnd(N, xcur, mu, sigma2, gamma1)[[3]], digits = 6)
  y1 = round(findyN_gnd(N, xcur, mu, sigma2, gamma1)[[4]], digits = 6)
  list(A = A1, x = x1, y = y1)
}


findyN_gnd <- function(N, xinit, mu, sigma2, gamma1) {
  x <- numeric(N)
  y <- numeric(N)
  
  flag = 0
  x[1] <- xinit
  tail <- 1 - generalized_normal_cdf(x[1], mu, sqrt(sigma2), gamma1)
  y[1] <- dgamma_gnd(x[1], mu, sigma2, gamma1)
  A <- x[1] * y[1] + tail
  
  x[2] = x[1]
  y[2] = y[1] + A / x[2]
  
  for (i in 2:N - 1) {
    if (y[i] > dgamma_gnd(0, mu, sigma2, gamma1)) {
      flag = -1
      break
    }
    
    x[i + 1] = gnd_inverse(y[i],mu, sigma2, gamma1)  
    y[i + 1] = y[i] + A / x[i]
  }
  
  if (flag == -1) {
    flag = 100
  } else if (abs(y[N] - dgamma_gnd(0, mu, sigma2, gamma1)) <= 1e-6) {
    flag = 1  # FOUND
  } else if (y[N] > dgamma_gnd(0, mu, sigma2, gamma1)) {
    flag = 100
  } else if (y[N] < dgamma_gnd(0, mu, sigma2, gamma1)) {
    flag = -100
  }
  return(list(flag, A, x, y))
}

GenFromTail_gnd <- function(a, mu, sigma2, gamma1) {
  u1 <- runif(1)
  u2 <- runif(1)
  while (u2 > a *  (a^(gamma1 / (gamma1 - 1)) - (gamma1 / (gamma1 - 1)) * log(u1))^( (1 - gamma1) / gamma1 )) {
    u1 <- runif(1)
    u2 <- runif(1)
  }
  #x = sqrt((a^2 - 2 * log(u1)))
  x =   (a^(gamma1 / (gamma1 - 1)) - (gamma1 / (gamma1 - 1)) * log(u1))^( (gamma1 - 1) / gamma1 )
  
  return(x)
}

# Modified Ziggurat algorithm for GND
myziggurat_gnd <- function(total, N, mu, sigma2, gamma1) {
  table <- findAxy_gnd(N, mu, sigma2, gamma1)
  #print("findAxy_gnd")
  
  xtable <- table$x
  ytable <- table$y
  Atable <- table$A
  
  x <- xtable
  y <- ytable
  A <- Atable
  immediate = 0
  rejection = 0
  R0 <- x[1] * y[1]
  data <- numeric(total)
  signpn <- (2 * rbinom(total, 1, 0.5) - 1)
  
  i = 1
  while (i <= total) {
    k <- sample(0:(N - 1), 1)
    u1 <- runif(1)
    
    if (k > 0) {
      xi <- u1 * x[k]
      if (xi < x[k + 1]) {
        data[i] = signpn[i] * xi
        immediate = immediate + 1
        i = i + 1
        
      } else {
        u2 <- runif(1)
        yi <- y[k] + u2 * (y[k + 1] - y[k])
        if (yi <= dgamma_gnd(xi, mu, sigma2, gamma1)) {
          data[i] = signpn[i] * xi
          i = i + 1
          
        } else {
          rejection = rejection + 1
        }
      }
    } else {
      w <- runif(1, max = A)
      if (w <= R0) {
        data[i] = signpn[i] * w / y[1]
        immediate = immediate + 1
        i = i + 1
        
      } else {
        data[i] = signpn[i] * GenFromTail_gnd(x[1], mu, sigma2, gamma1)
        i = i + 1
      }
    }
  }
  
  print(noquote(paste("rejection rate:", rejection / total)))
  print(noquote(paste("immediate accept rate:", immediate / total)))
  return(data)
}



# Define the Multivariate Generalized Normal Density Function
mvgnd_pdf <- function(x, mu, Sigma, gamma1) {
  p <- length(mu)  # Dimension of the distribution
  
  
  # Regularize Sigma if it is singular
  epsilon <- 1e-6  # Small value for numerical stability
  if (!is.positive.definite(Sigma)) {
    Sigma <- Sigma + diag(epsilon, p)
  }
  
  # Compute determinant and inverse of Sigma
  det_Sigma <- det(Sigma)
  inv_Sigma <- solve(Sigma)
  
  # Compute quadratic form Q(x)
  Q_x <- t(x - mu) %*% inv_Sigma %*% (x - mu)
  
  # Compute the normalizing constant C_gamma^p
  C_gamma_p <- (pi^(-p / 2)) * 
    (gamma(p / 2 + 1) / gamma((p * (gamma1 - 1) / gamma1) + 1)) * 
    ((gamma1 - 1) / gamma1)^(p * (gamma1 - 1) / gamma1)
  
  # Compute density function
  density <- C_gamma_p * (det_Sigma^(-1/2)) * 
    exp(-((gamma1 - 1) / gamma1) * Q_x^(gamma1 / (2 * (gamma1 - 1))))
  
  return(density)
}
