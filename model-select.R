# Bayesian Final Project - Fall 2017
# Model Selection 

# See model.select() function for running instructions
require(MASS)
cat("To run: model <- model.select(X_train, y_train, iters).\n") 
cat("'model' here is a logical vector of size equal to the number of columns of X_train.\n")
cat("Iters controls the number of samples of the Gibb's Sample.\n")

# Define Likelihood function.
lpy.X <- function(y, X) {
  n <- dim(X)[1];
  p <- dim(X)[2];
  g <- length(y);
  nu_0 <- 1;
  theta <- ginv(t(X) %*% X) %*% t(X) %*% y
  err <- y - X %*% theta
  sigma2_0 <- sum(err^2)/length(err)
  if (p == 0) {
    H_g <- 0;
    sigma2_0 <- mean(y^2);
  }
  if (p > 0) {
    H_g <- (g/(g+1))*X%*%ginv(t(X)%*%X)%*%t(X);
  }
  SSR_g <- t(y)%*%(diag(1, nrow = n)-H_g)%*%y;
  -(1/2)*(n*log(pi)+p*log(1+g)+(nu_0+n)*log(nu_0*sigma2_0+SSR_g)-nu_0*log(nu_0*sigma2_0))+
    lgamma((nu_0+n)/2)-lgamma(nu_0/2);
}

# Define Model Selection
model.select <- function(X, y, iters) {
  # From Hoff, page 168
  nu_0 <- 1;
  # Set sigma with ordinary leaste squares error.
  theta <- ginv(t(X) %*% X) %*% t(X) %*% y
  err <- y - X %*% theta
  sigma2_0 <- sum(err^2)/length(err)
  g <- n
  t <- iters;
  
  # Run MCMC 
  p <- dim(X)[2]
  Z <- matrix(NA, t, p);
  B <- matrix(0, t, p);
  z <- rep(1, p);
  lpy.c <- lpy.X(y, X[, z == 1, drop = FALSE])
  
  cat("Running Gibb's Sampler for model selection!")
  for(i in 1:t) {
    cat("Sample", i, "\n")
    for(j in sample(1:p)) {
      zp <- z;
      zp[j] <- 1-zp[j];
      lpy.p <- lpy.X(y, X[, zp == 1, drop = FALSE]);
      r <- (lpy.p-lpy.c)*(-1)^(zp[j] == 0);
      z[j] <- rbinom(1, 1, 1/(1+exp(-r)));
      if (z[j] == zp[j]) {
        lpy.c <- lpy.p;
      }
      Z[i, ] <- z;
      H_g <- (g/(g+1))*X[, Z[i,] == 1, drop = FALSE]%*%ginv(t(X[, Z[i,] == 1, drop = FALSE])%*%
                                                              X[, Z[i,] == 1, drop = FALSE])%*%t(X[, Z[i,] == 1, drop = FALSE]);
      SSR_g <- t(y)%*%(diag(1, nrow = n)-H_g)%*%y;
      sigma2 <- 1/rgamma(t, (nu_0+n)/2, (nu_0*sigma2_0+SSR_g)/2);
      Vb <- g*ginv(t(X[, Z[i,] == 1, drop = FALSE])%*%X[, Z[i,] == 1, drop = FALSE])/(g+1);
      Eb <- Vb%*%t(X[, Z[i,] == 1, drop = FALSE])%*%y;
      E <- matrix(rnorm(sum(Z[i, ]), 0, sqrt(sigma2)), 1, sum(Z[i, ]));
      B[i, Z[i,] == 1] <- t(t(E%*%chol(Vb, pivot = TRUE))+c(Eb));
    }
  }
  
  model <- rep(0, p)
  for (i in 1:p) {
    model[i] <- sum(B[,i] != 0)/t;
  }
  model.subset <- model > 0.95
  cat("Trained model has", sum(model.subset), "regressors.")
  return(model.subset)
}









