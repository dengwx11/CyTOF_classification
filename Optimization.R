library(truncnorm)

update_w <- function(i, j, X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu) {
  numerator <- X %*% t(H_prev) + lambda1 * A %*% S + lambda2 * A0
  denominator <- W_prev %*% H_prev %*% t(H_prev) + 2 * mu * W_prev
  return(W_prev[i, j] * numerator[i, j] / denominator[i, j])
}

update_h <- function(i, j, W, X, H_prev) {
  numerator <- t(W) %*% X
  denominator <- t(W) %*% W %*% H_prev
  return(H_prev[i, j] * numerator[i, j] / denominator[i, j])
}

norm_h <- function(i, j, H) {
  return(H[i, j] / sum(H[, j]))
}

update_w_mat <- function(X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu) {
  w_row <- nrow(W_prev)
  w_col <- ncol(W_prev)
  W <- matrix(, nrow = w_row, ncol = w_col)
  for (i in 1:w_row) {
    for (j in 1:w_col) {
      W[i, j] <- update_w(i, j, X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu)
      W_prev[i, j] <- W[i, j]
    }
  }
  return(W)
}

update_h_mat <- function(W, X, H_prev) {
  h_row <- nrow(H_prev)
  h_col <- ncol(H_prev)
  H <- matrix(, nrow = h_row, ncol = h_col)
  for (j in 1:h_col) {
    for (i in 1:h_row) {
      H[i, j] <- update_h(i, j, W, X, H_prev)
      H_prev[i, j] <- H[i, j]
    }
    for (i in 1:h_row) {
      H[i, j] <- norm_h(i, j, H)
    }
  }
  return(H)
}

loss_func <- function(X, A, S, A0, H, W, lambda1, lambda2, mu) {
  decov <- 1 / 2 * sum(diag(t(X - W %*% H) %*% (X - W %*% H)))
  trace <- lambda1 * sum(diag(t(W) %*% A %*% S)) + lambda2 * sum(diag(t(W) %*% A0))
  regu <- mu * sum(diag(t(W) %*% W))
  return(decov - trace + regu)
}
 
optim_wh <- function(X, A, S, A0, H0, W0, lambda1, lambda2, mu) {
  H_prev <- H0
  W_prev <- W0
  loss_prev <- loss_func(X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu)
  loss_save <- c()
  repeat {
    W <- update_w_mat(X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu)
    H <- update_h_mat(W, X, H_prev)
    loss <- loss_func(X, A, S, A0, H, W, lambda1, lambda2, mu)
    loss_save <- c(loss_save, loss)
    if (abs(loss - loss_prev) < 1e-4) {
      break
    }
    H_prev <- H
    W_prev <- W
    loss_prev <- loss
  }
  return(list(W, H, loss_save))
}

# Wang inhouse
S <- read.csv('/home/bz234/project/Objects/Optimization//Wang/S_wang_inhouse.csv', row.names = 1)
X <- read.csv('/home/bz234/project/Objects/Optimization/Wang/X_wang_inhouse.csv', row.names = 1)
A <- read.csv('/home/bz234/project/Objects/Optimization/Wang/A_wang_inhouse.csv', row.names = 1)
A0 <- read.csv('/home/bz234/project/Objects/Optimization/Wang/A0_wang_inhouse.csv', row.names = 1)

S <- as.matrix(S)
X <- as.matrix(X)
A <- as.matrix(A)
A0 <- as.matrix(A0)

w_row <- nrow(X)
w_col <- ncol(S)
W0 <- matrix(nrow = w_row, ncol = w_col)
for (i in 1:w_row) {
    W0[i, ] <- rtruncnorm(w_col, a=0, b=Inf, mean = 1, sd = .2)
}

h_row <- ncol(S)
h_col <- ncol(X)
H0 <- matrix(nrow = h_row, ncol = h_col)
for (i in 1:h_row) {
    H0[i, ] <- rtruncnorm(h_col, a=0)
}
H0 <- apply(H0, 2, function(x) x / sum(x))

lambda1 <- .5
lambda2 <- .5
mu <- .5

results <- optim_wh(X, A, S, A0, H0, W0, lambda1, lambda2, mu)
write.csv(results[[1]], '/home/bz234/project/Objects/Optimization/Wang/W_inhouse_optim.csv')
write.csv(results[[2]], '/home/bz234/project/Objects/Optimization/Wang/H_inhouse_optim.csv')
write.csv(results[[3]], '/home/bz234/project/Objects/Optimization/Wang/Loss_inhouse_optim.csv')
