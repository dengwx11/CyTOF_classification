#source('Optimization.v1.R')
source('./Optimization.v2.R')


run <- function(X,lambda1,lambda2,mu, eta, AS, A0, D, K, N, epsilon = 10^(-4),fixed_loop = 0) {
    ## Initialization
    W <- matrix(1,nrow = D, ncol = K)
    H <- matrix(1,nrow = K, ncol = N)


    L.current <- compute_L(X,W,H,lambda1,lambda2,mu,eta)
    L.min <- L.current
    L.prev <- 100000
    L.save <- c(L.current)
    k=1
     while(abs(L.prev - L.current) > epsilon){
    #while(abs(L.prev - L.current) > epsilon && L.min >= L.current){
    #for(i in 1:1000){
        W <- update_W(X, W, H, mu, AS, A0, D, K, lambda1, lambda2)
        L.prev <- L.current
        L.current <- compute_L(X,W,H,lambda1,lambda2,mu,eta)
        #print(paste0("update W, L change ", L.current-L.prev))
        H <- update_H(X, W, H, K, N, D, eta)
        L.prev <- L.current
        L.current <- compute_L(X,W,H,lambda1,lambda2,mu,eta)
        #print(paste0("update H, L change ", L.current-L.prev))
        L.save <- c(L.save, L.current)
        if(L.min>L.current) L.min=L.current
        #print(L.current)
        k = k+1
        if(fixed_loop >0 && k>fixed_loop) break
    }

    print(paste0("iteration = ",k))

    rst <- list()
    rst$W <- W
    rst$H <- H
    rst$L <- L.current
    rst$L.save <- L.save

    return(rst)

}

