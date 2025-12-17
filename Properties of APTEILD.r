library(pracma)
eps<- 1e-12
# CDF of APTEILD
pAPTEILD <- function(x, alpha, theta, beta, delta){
  x <- pmax(x, eps)
  
  expo <- (1 + (theta * beta) / (theta + beta) * x^(-delta)) *
          exp(-theta / x^delta)
  
  if (abs(alpha - 1) < 1e-8) {
    F <- expo
  } else {
    F <- (exp(log(alpha) * expo) - 1) / (alpha - 1)
  }
  
  pmin(pmax(F, 0), 1)
}

# PDF of APTEILD
dAPTEILD <- function(x, alpha, theta, beta, delta){
  x <- pmax(x, eps)
  
  A <- log(alpha) / (alpha - 1)
  B <- delta * theta^2 / (theta + beta)
  C <- (beta + x^delta) / x^(2 * delta + 1)
  D <- exp(-theta / x^delta)
  
  expo <- (1 + (theta * beta) / (theta + beta) * x^(-delta)) * D
  
  dens <- A * B * C * D * alpha^expo
  dens[dens <= 0 | !is.finite(dens)] <- 0
  dens
}


# Normalization check
norm_cond <- function(alpha, theta, beta, delta){
  integrate(function(x)
    dAPTEILD(x, alpha, theta, beta, delta),
    lower = 0, upper = Inf, rel.tol = 1e-10)$value
}

# Survival function
sAPTEILD <- function(x, alpha, theta, beta, delta){
  1 - pAPTEILD(x, alpha, theta, beta, delta)
}

# Hazard rate function
hAPTEILD <- function(x, alpha, theta, beta, delta){
  dAPTEILD(x, alpha, theta, beta, delta) /
    sAPTEILD(x, alpha, theta, beta, delta)
}

# Quantile function (numerical)
qAPTEILD <- function(p, alpha, theta, beta, delta){
  sapply(p, function(pp){
    uniroot(function(x)
      pAPTEILD(x, alpha, theta, beta, delta) - pp,
      lower = eps, upper = 1e6, tol = 1e-10)$root
  })
}

# Ordinary moments
mAPTEILD <- function(k, alpha, theta, beta, delta){
  integral(function(x)
    x^k * dAPTEILD(x, alpha, theta, beta, delta),
    0, Inf, reltol = 1e-10, method = "Simpson")
}


# Incomplete moments
imAPTEILD <- function(k, t, alpha, theta, beta, delta){
  integral(function(x)
    x^k * dAPTEILD(x, alpha, theta, beta, delta),
    0, t, reltol = 1e-10, method = "Simpson")
}
# Mean, variance, CV, skewness, kurtosis
chAPTEILD <- function(alpha, theta, beta, delta){
  mu  <- mAPTEILD(1, alpha, theta, beta, delta)
  mu2 <- mAPTEILD(2, alpha, theta, beta, delta)
  mu3 <- mAPTEILD(3, alpha, theta, beta, delta)
  mu4 <- mAPTEILD(4, alpha, theta, beta, delta)
  
  var <- mu2 - mu^2
  cv  <- sqrt(var) / mu
  
  skew <- (mu3 - 3 * mu * mu2 + 2 * mu^3) / var^(3/2)
  kurt <- (mu4 - 4 * mu * mu3 + 6 * mu^2 * mu2 - 3 * mu^4) / var^2
  
  c(mean = mu,
    variance = var,
    CV = cv,
    skewness = skew,
    kurtosis = kurt)
}
# Moment generating function
mgfAPTEILD <- function(t, alpha, theta, beta, delta){
  integral(function(x)
    exp(t * x) * dAPTEILD(x, alpha, theta, beta, delta),
    0, Inf, reltol = 1e-10, method = "Simpson")
}
