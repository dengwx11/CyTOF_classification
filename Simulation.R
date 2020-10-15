set.seed(2020)
library(Seurat)

K = 5 # cell types number ## K could be larger
D = 10 # surface markers number
N = 500 # ADT/CyTOF cell number ## N could be larger
G = 100 # RNA gene number
pi_ber1 = 0.55
pi_ber2 = 0.9

#generate gamma
generate_gamma <- function(D, K, pi_ber1, pi_ber2){
    gamma <- matrix(rbinom(n=D*K, size = 1, prob = pi_ber1), nrow = D, ncol = K)
    gamma <- gamma * (matrix(rbinom(n=D*K, size = 1, prob = pi_ber2), nrow = D, ncol = K)-0.5) * 2
    return(gamma)
}

max_cos <- function(w){
    K <- ncol(w)
    w_cos <- matrix(0, ncol= K, nrow = K)
    for (i in 1:K){
        for (j in 1:K){
            if (i==j){
                w_cos[i,j] = -1
            } else {
                w_cos[i,j] = (sum(w[,i]*w[,j]))/(sqrt(sum((w[,i]^2)))*sqrt(sum((w[,j]^2))))
            }
        }
    }
    max_cos = max(w_cos)
    return(max_cos)
}

simulate_gamma <- function(iteration, K, D){  ## iteration could be set as 1,3,10
    gamma.list <- list()
    cos.list <- list()
    output <- list()
    for (i in 1:iteration){
        gamma <- generate_gamma(D, K, pi_ber1, pi_ber2)
        gamma.list[[i]] <- gamma
        cos.list[[i]] <- max_cos(gamma)
    }
    minindex <- which.min(cos.list)
    print(cos.list[[minindex]])
    # output$w <- w.list[[minindex]]
    output$gamma <- gamma.list[[minindex]]
    output$max_cos <- cos.list[[minindex]]

    return(output)
}

simulate_w_pre <- function(D,K, big_w_mean=2, big_tau_w = 10, small_w_mean=0.5, small_tau_w = 10){
    
    big_w <- matrix(rnorm(D*K, mean = big_w_mean, sd = 1/sqrt(big_tau_w)), nrow = D, ncol = K)
    small_w <- matrix(rnorm(D*K, mean = small_w_mean, sd = 1/sqrt(small_tau_w)), nrow = D, ncol = K)
    
    output <- list()
    output$big_w <- big_w
    output$small_w <- small_w

    return(output)
}

simulate_w <- function(D, K, gamma, p.0 = 0.1, q.0=.2, p.neg1= 0.05, q.neg1 = 0.1){
    gamma_post <- matrix(0, nrow = D, ncol = K)
    gamma_post[which(gamma==0)] <- sample(c(0:2), size = length(which(gamma==0)), 
                                        replace = T, prob = c(1-p.0-q.0, q.0, p.0))
    gamma_post[which(gamma==1)] <- sample(c(0:2), size = length(which(gamma==1)), 
                                        replace = T, prob = c(0, 0, 1))
    gamma_post[which(gamma==-1)] <- sample(c(0:2), size = length(which(gamma==-1)), 
                                        replace = T, prob = c(1-p.neg1-q.neg1, q.neg1, p.neg1))                                                                           

    w_pre <- simulate_w_pre( D, K)
    w <- matrix(0, nrow = D, ncol = K)
    w[which(gamma_post == 2)] <- w_pre$big_w[which(gamma_post == 2)]
    w[which(gamma_post == 1)] <- w_pre$small_w[which(gamma_post == 1)]
    w[which(gamma_post == 0)] <- 0.1

    return(w)
}

simulate_label <- function(D,N,K, prob_k = c(1,2,2,3,3)){
    label <- matrix(sample(c(1:K),size = N, replace = T, prob = prob_k), nrow = 1, ncol = N)
    H <- matrix(0, nrow = D, ncol = N)
    H <- sapply(c(1:N), function(i) {
        h = matrix(0,nrow = K, ncol = 1)
        h[label[1,i]] =  1
        return(h)
    })
    output <- list()
    output$label <- label
    output$H <- H
    return(output)
}

simulate_x <- function(D,N,w,label, mean_var_ratio = 10){ ## 'mean_var_ratio' could be adjusted
    X <- matrix(0, nrow = D, ncol = N)
    for(i in 1:N){
        celltype <- label[1,i]
        for(j in 1:D){
            X[j,i] <- rnorm(1, mean=w[j,celltype], sd = sqrt(w[j,celltype]/mean_var_ratio) )  
        }
    }
    return(X)
}

simulate_AS <- function(D,K,w,corr = 0.5){  ## 'corr' could be adjusted
    AS <- matrix(0, nrow = D, ncol = K)
    for(i in 1:D){
        for(j in 1:K){
            AS[i,j] <- rnorm(1, mean = w[i,j], sd = w[i,j]*0.5/corr)
        }
    }
    return(AS)
}


