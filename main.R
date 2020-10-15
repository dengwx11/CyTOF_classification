source('Optimization.v1.R')


run <- function(X,lambda1,lambda2,mu, AS, A0, D, K, N, epsilon = 10^(-4)) {
    ## Initialization
    W <- matrix(1,nrow = D, ncol = K)
    H <- matrix(1,nrow = K, ncol = N)


    L.current <- compute_L(X,W,H,lambda1,lambda2,mu)
    L.prev <- 100000
    L.save <- c(L.current)
    while(abs(L.prev - L.current) > epsilon){
    #for(i in 1:1000){
        W <- update_W(X, W, H, mu, AS, A0, D, K, lambda1, lambda2)
        H <- update_H(X, W, H, K, N)
        L.prev <- L.current
        L.current <- compute_L(X,W,H,lambda1,lambda2,mu)
        L.save <- c(L.save, L.current)
        print(L.current)
    }

    rst <- list()
    rst$W <- W
    rst$H <- H
    rst$L <- L.current
    rst$L.save <- L.save

    return(rst)

}

rst<-run(X,1,1,1,AS,A0,D,K,N)
plot(as.vector(W),as.vector(rst$W))
plot(as.vector(true.H),as.vector(rst$H))