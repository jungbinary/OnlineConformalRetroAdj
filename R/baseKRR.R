# ----- Utility Functions -----
compute_dist_sq <- function(x, y) {
  outer(rowSums(x^2), rowSums(y^2), '+') - 2 * tcrossprod(x, y)
}

.compute_ntk <- function(X1, X2, eps = 1e-12) {
  # ReLU 2-layer NTk
  
  # norms
  n1 <- sqrt(rowSums(X1^2))
  n2 <- sqrt(rowSums(X2^2))
  
  n1n2 <- n1 %o% n2
  n1n2_safe <- pmax(n1n2, eps)
  
  # cosine and angle
  cos_th <- tcrossprod(X1, X2) / n1n2_safe
  cos_th <- pmin(pmax(cos_th, -1), 1)   # clamp to [-1,1]
  th <- acos(cos_th)
  
  # components
  Sigma     <- (n1n2 / pi) * (sin(th) + (pi - th) * cos_th)
  Sigma_dot <- (pi - th) / pi
  
  Sigma[n1n2 < eps] <- 0
  
  K <- Sigma + Sigma_dot
  return(K)
}

# ----- Main Kernel Switch -----
compute_kernel_mat <- function(x, y, sigma = 1.0,
                               kernel = c("rbf", "ntk")) {
  kernel <- match.arg(kernel)
  
  if (kernel == "rbf") {
    # Gaussian/RBF kernel
    D2 <- compute_dist_sq(x, y)
    return(exp(- D2 / (2 * sigma^2)))
  } else if (kernel == "ntk") {
    # 2-layer ReLU NTK
    return(.compute_ntk(x, y))
  }
}

# KRR -----
krr_init <- function(X_init, y_init, kernel = "rbf", sigma = NULL, lambda = NULL){
  stopifnot(is.matrix(X_init) || is.data.frame(X_init))
  X_init <- as.matrix(X_init)
  n <- nrow(X_init); d <- ncol(X_init)
  stopifnot(length(y_init) == n)
  
  # --- sigma grid (RBF 전용) ---
  if (tolower(kernel) == "rbf") {
    if (is.null(sigma)) {
      # median heuristic
      dist_vec <- as.vector(dist(X_init))^2
      sigma_med <- sqrt(median(dist_vec))
      sigma_grid <- sigma_med * 2^seq(-3, 3, length.out = 10)
    } else {
      sigma_grid <- as.numeric(sigma)
    }
  } else {
    sigma_fixed <- sqrt(d)
    sigma_grid <- sigma_fixed
  }
  
  # --- lambda grid ---
  lambda_grid <- if (is.null(lambda)) 10^seq(-5, 0, length.out = 5) else as.numeric(lambda)
  
  best_loocv <- Inf
  best_model <- NULL
  
  for (s in sigma_grid){
    K <- compute_kernel_mat(X_init, X_init, sigma = s, kernel = kernel)
    
    for (l in lambda_grid){
      # Q = (K + lambda I)^{-1}
      Q <- solve(K + l * diag(n))
      
      # LOO residuals: r_i = (Q y)_i / Q_ii
      r_loo    <- as.vector((Q %*% y_init) / diag(Q))
      loocv_er <- mean(r_loo^2)
      
      if (loocv_er < best_loocv) {
        best_loocv <- loocv_er
        best_model <- list(
          X      = X_init,
          y      = y_init,
          Q      = Q,
          lambda = l,
          sigma  = s,
          kernel = kernel,
          cv_mse = loocv_er
        )
      }
    }
  }
  
  return(best_model)
}


krr_update <- function(model, x_new, y_new){
  X       <- model$X
  Q       <- model$Q
  lambda  <- model$lambda
  sigma   <- model$sigma        
  kernel  <- model$kernel      
  
  if (is.null(dim(x_new))) x_new <- matrix(x_new, nrow = 1)
  
  k_x  <- compute_kernel_mat(X, x_new, sigma, kernel = kernel)  
  k_xx <- as.numeric(compute_kernel_mat(x_new, x_new, sigma, kernel = kernel)) 
  
  Qk_x <- Q %*% k_x
  delta <- as.numeric(1 / (k_xx + lambda - t(k_x) %*% Qk_x))
  
  B12 <- -delta * Qk_x
  B11 <- Q + delta * (Qk_x %*% t(Qk_x))
  
  Qnew <- rbind(
    cbind(B11, B12),
    cbind(t(B12), delta)
  )
  
  model$X <- rbind(X, x_new)
  model$y <- c(model$y, y_new)
  model$Q <- Qnew
  return(model)
}


krr_downdate <- function(model){
  
  Q <- model$Q
  
  q11  <- Q[1, 1]
  q12  <- Q[-1, 1, drop = FALSE]
  Q22  <- Q[-1, -1, drop = FALSE]
  
  Qnew <- Q22 - (q12 %*% t(q12)) / q11
  
  model$X <- model$X[-1, ,drop = FALSE]
  model$y <- model$y[-1]
  model$Q <- Qnew
  
  return(model)
}

# Jackknife+ Residual -----
compute_resid_jk <- function(model, x_new){
  
  X <- model$X
  y <- model$y
  Q <- model$Q
  sigma <- model$sigma
  lambda <- model$lambda
  kernel <- model$kernel
  
  S <- diag(nrow(Q)) - lambda * Q
  
  y_hat <- as.numeric(S %*% y)
  S_diag <- diag(S)
  R_signed <- (y - y_hat) / pmax(1e-8, 1 - S_diag)
  
  k_x <- compute_kernel_mat(X, x_new, sigma, kernel = kernel)         # (n x 1)
  if(is.vector(k_x)) k_x <- matrix(k_x, ncol = 1)
  
  S_new <- as.numeric(t(k_x) %*% Q)                  # (1 x n)
  f_new <- as.numeric(S_new %*% y)                   # scalar
  
  f_new_loo <- f_new - S_new * R_signed
  L <- f_new_loo - abs(R_signed)
  U <- f_new_loo + abs(R_signed)
  
  return(list(R = abs(R_signed), L = L, U = U))
}

# Compute Beta with Jackknife+ Residual
compute_beta_jk <- function(y_new, L, U) {
  
  F_L <- mean(L <= y_new)
  F_U <- mean(U <= y_new)
  
  beta_t <- pmax(0, 
                 pmin(1, F_L, (1 -F_U))
  )
  
  return(beta_t)
}

krr_predict <- function(model, x_new){
  with(model, {
    k_x <- compute_kernel_mat(X, x_new, sigma = sigma, kernel = kernel)         # (n x 1)
    if(is.vector(k_x)) k_x <- matrix(k_x, ncol = 1)
    
    S_new <- as.numeric(t(k_x) %*% Q)                  # (1 x n)
    as.numeric(S_new %*% y)
  })
}

