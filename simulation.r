library(pracma)
library(doParallel)
library(foreach)
vCDF_APTEILD <- function(x, alpha, theta, beta, delta){
  x <- pmax(x, 1e-12)
  t1 <- 1 + (theta*beta/(theta+beta)) * (1 / (x^delta))
  t2 <- exp(-theta / (x^delta))
  expo <- t1 * t2
  Aexp <- exp(expo * log(alpha))   
  cdf <- (Aexp - 1) / (alpha - 1)
  cdf[cdf < 0] <- 0
  cdf[cdf > 1] <- 1
  return(cdf)
}

vPDF_APTEILD <- function(x, alpha, theta, beta, delta){
  x <- pmax(x, 1e-12)
  
  A <- log(alpha)/(alpha - 1)
  B <- delta * theta^2/(theta + beta)
  C <- (beta + x^delta)/(x^(2*delta+1))
  D <- exp(-theta/(x^delta))
  expo <- (1 + (theta*beta/(theta+beta))*(1/(x^delta))) * D
  
  pdf <- A * B * C * D * exp(expo * log(alpha))
  
  pdf[!is.finite(pdf)] <- 0
  pdf[pdf < 0] <- 0
  return(pdf)
}

build_inv_cdf_APTEILD <- function(alpha, theta, beta, delta, x_max=50, n_grid=30000){
  xg <- seq(1e-8, x_max, length.out = n_grid)
  cdfg <- vCDF_APTEILD(xg, alpha, theta, beta, delta)
  keep <- !duplicated(cdfg)
  list(x=xg[keep], cdf=cdfg[keep])
}

inv_sample_APTEILD <- function(u, invobj){
  approx(invobj$cdf, invobj$x, xout=u, rule=2, ties="ordered")$y
}

rAPTEILD <- function(n, alpha, theta, beta, delta, invobj=NULL){
  if(is.null(invobj))
    invobj <- build_inv_cdf_APTEILD(alpha,theta,beta,delta)
  u <- runif(n)
  inv_sample_APTEILD(u, invobj)
}


# Objective Functions (4 parameters)

MLE_obj <- function(par, data){
  a=par[1]; t=par[2]; b=par[3]; d=par[4]
  if(any(!is.finite(par)) || any(par <= 0)) return(1e12)
  dens <- vPDF_APTEILD(data,a,t,b,d)
  if(any(dens<=0)) return(1e10 + sum(dens<=0))
  return(-sum(log(dens + 1e-12)))
}

OLS_obj <- function(par, data){
  a=par[1]; t=par[2]; b=par[3]; d=par[4]
  if(any(!is.finite(par)) || any(par<=0)) return(1e12)
  xs <- sort(data)
  n <- length(xs)
  ecdf <- (1:n)/(n+1)
  tcdf <- vCDF_APTEILD(xs,a,t,b,d)
  sum((ecdf - tcdf)^2)
}

WLS_obj <- function(par,data){
  a=par[1]; t=par[2]; b=par[3]; d=par[4]
  xs <- sort(data)
  n <- length(xs)
  ecdf <- (1:n)/(n+1)
  tcdf <- vCDF_APTEILD(xs,a,t,b,d)
  w <- 1/(ecdf*(1-ecdf)+0.01)
  sum(w*(ecdf-tcdf)^2)
}

CVM_obj <- function(par,data){
  a=par[1]; t=par[2]; b=par[3]; d=par[4]
  xs <- sort(data)
  n <- length(xs)
  ecdf <- ((1:n)-0.5)/n
  tcdf <- vCDF_APTEILD(xs,a,t,b,d)
  sum((ecdf-tcdf)^2) + 1/(12*n)
}

MPS_obj <- function(par,data){
  a=par[1]; t=par[2]; b=par[3]; d=par[4]
  xs <- sort(data)
  tcdf <- vCDF_APTEILD(xs,a,t,b,d)
  sp <- c(tcdf[1], diff(tcdf), 1 - tcdf[length(tcdf)])
  if(any(sp <= 0)) return(1e10)
  return(-sum(log(sp)))
}


fast_opt <- function(fn, data, lower, upper){
  starts <- list(c(1.5,1,1,1),
                 c(5, 0.5, 1.5, 1),
                 c(1,2,1,1.5))
  best_val <- Inf; best_par <- NULL
  for(st in starts){
    out <- tryCatch(
      optim(st, fn=function(p) fn(p,data), method="L-BFGS-B",
            lower=lower, upper=upper,
            control=list(maxit=200)),
      error=function(e) NULL
    )
    if(!is.null(out) && out$value < best_val){
      best_val <- out$value
      best_par <- out$par
    }
  }
  return(list(par=best_par, value=best_val))
}

g <- 100
sample_sizes <- c(50,100,200,300,500)
true.alpha <- 2.0
true.theta <- 1.5
true.beta  <- 0.5
true.delta <- 1

lower <- c(0.1,0.01,0.01,0.01)
upper <- c(5,5,5,5)

ncores <- max(1, parallel::detectCores()-1)
cl <- makeCluster(ncores)
registerDoParallel(cl)

invobj_true <- build_inv_cdf_APTEILD(true.alpha,true.theta,true.beta,true.delta)

summary_results <- list()

for(n in sample_sizes){
  cat("\nRunning n =",n,"\n")
  
  res <- foreach(j=1:g, .combine='rbind') %dopar% {
    sim <- rAPTEILD(n,true.alpha,true.theta,true.beta,true.delta,invobj_true)
    
    mle <- fast_opt(MLE_obj, sim, lower, upper)$par
    ols <- fast_opt(OLS_obj, sim, lower, upper)$par
    wls <- fast_opt(WLS_obj, sim, lower, upper)$par
    cvm <- fast_opt(CVM_obj, sim, lower, upper)$par
    mps <- fast_opt(MPS_obj, sim, lower, upper)$par
    
    c(mle, ols, wls, cvm, mps)
  }
  
  res_mat <- matrix(res, ncol=20, byrow=TRUE)
  
  get_est <- function(i) matrix(as.numeric(res_mat[,i:(i+3)]), ncol=4)
  
  mleE <- get_est(1)
  olsE <- get_est(5)
  wlsE <- get_est(9)
  cvmE <- get_est(13)
  mpsE <- get_est(17)
  
  calc_metrics <- function(est){
    est <- est[apply(est,1,function(z) all(is.finite(z))),]
    if(nrow(est)==0) return(NULL)
    
    means <- colMeans(est)
    bias  <- means - c(true.alpha,true.theta,true.beta,true.delta)
    vars  <- apply(est,2,var)
    mse   <- vars + bias^2
    
    list(mean=means, bias=bias, mse=mse)
  }
  
  summary_results[[as.character(n)]] <- list(
    MLE = calc_metrics(mleE),
    OLS = calc_metrics(olsE),
    WLS = calc_metrics(wlsE),
    CVM = calc_metrics(cvmE),
    MPS = calc_metrics(mpsE)
  )
  
  print(summary_results[[as.character(n)]])
}

stopCluster(cl)
cat("\nSimulation Finished.\n").