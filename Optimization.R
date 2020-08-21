update_w <- function(i, j, X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu) {
  numerator <- X %*% t(H_prev) + lambda1 * A %*% S + lambda2 * A0
  denominator <- W %*% H_prev %*% t(H_prev) + 2 * mu * W_prev
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
    }
    for (i in 1:h_row) {
      H[i, j] <- norm_h(i, j, H)
    }
  }
  return(H)
}

loss_func <- function(X, A, S, A0, H, W, lambda1, lambda2, mu) {
  decov <- sum(diag(t(X - W %*% H) %*% (X - W %*% H)))
  trace <- lambda1 * sum(diag(t(W) %*% A %*% S)) + lambda2 * sum(diag(t(W) %*% A0))
  regu <- mu * sum(diag(t(W) %*% W))
  return(decov - trace + regu)
}
 
optim_wh <- function(X, A, S, A0, H0, W0, lambda1, lambda2, mu) {
  H_prev <- H0
  W_prev <- W0
  loss_prev <- loss_func(X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu)
  repeat {
    W <- update_w_mat(X, A, S, A0, H_prev, W_prev, lambda1, lambda2, mu)
    H <- update_h_mat(W, X, H_prev)
    loss <- loss_func(X, A, S, A0, H, W, lambda1, lambda2, mu)
    if (abs(loss - loss_prev) < 1e-4) {
      break
    }
    H_prev <- H
    W_prev <- W
    loss_prev <- loss
  }
  return(list(W, H))
}

