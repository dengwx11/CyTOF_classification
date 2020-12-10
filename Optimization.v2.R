update_w <- function(x, w, H, mu, as, a0, K, lambda1, lambda2){
    term1 <- w %*% (H%*%t(H)+mu*diag(K))
    term2 <- x %*% t(H)
    term3 <- lambda1 * as
    term4 <- lambda2 * a0

    term1.pos <- matrix(0, nrow =1, ncol = K)
    term2.pos <- matrix(0, nrow =1, ncol = K)
    term3.pos <- matrix(0, nrow =1, ncol = K)
    term4.pos <- matrix(0, nrow =1, ncol = K)
    term1.neg <- matrix(0, nrow =1, ncol = K)
    term2.neg <- matrix(0, nrow =1, ncol = K)
    term3.neg <- matrix(0, nrow =1, ncol = K)
    term4.neg <- matrix(0, nrow =1, ncol = K)

    term1.pos[1, which(term1>=0)] <- term1[ which(term1>=0) ]
    term2.pos[1, which(term2>=0)] <- term2[ which(term2>=0) ]
    term3.pos[1, which(term3>=0)] <- term3[ which(term3>=0) ]
    term4.pos[1, which(term4>=0)] <- term4[ which(term4>=0) ]
    term1.neg[1, which(term1<0)] <- -term1[ which(term1<0) ]
    term2.neg[1, which(term2<0)] <- -term2[ which(term2<0) ]
    term3.neg[1, which(term3<0)] <- -term3[ which(term3<0) ]
    term4.neg[1, which(term4<0)] <- -term4[ which(term4<0) ]

    new_w <- w * (term1.neg + term2.pos + term3.pos + term4.pos)/(term1.pos + term2.neg + term3.neg+term4.neg )

    return(new_w)
}

update_h <- function(x,h,W, eta){
    one_k <- matrix(1, nrow = nrow(h),ncol=1)
    one_kk <- matrix(1, nrow = nrow(h), ncol = nrow(h))

    ## model 0
    #new_h <- h * ( t(W) %*% x  )/(t(W) %*% W %*% h + 10^(-16))
    #new_h <- new_h/sum(new_h)

    ## model 1
    new_h <- h * ( t(W) %*% x + eta * one_k )/(t(W) %*% W %*% h + eta * one_kk %*% h )

    ## model 2
    #new_h <- h * ( t(W) %*% x )/(t(W) %*% W %*% h + eta * one_k + 10^(-16))
    
    return(new_h)
}

update_W <- function(X, W, H, mu, AS, A0, D, K, lambda1, lambda2){
    for(i in 1:D){
        x = matrix(X[i,],nrow=1,ncol=N)
        w = matrix(W[i,],nrow=1,ncol=K)
        as = matrix(AS[i,],nrow=1,ncol=K)
        a0 = matrix(A0[i,],nrow=1,ncol=K)
        W[i,] <- update_w(x, w, H, mu, as, a0, K, lambda1, lambda2)
    }


    return(W)
}

update_H <- function(X, W, H, K, N, eta){

    for(j in 1:N){
        x = matrix(X[,j],nrow=D,ncol=1)
        h = matrix(H[,j],nrow=K,ncol=1)
        H[,j] <- update_h(x,h,W,eta)
    }

    return(H)
}

compute_L <- function(X,W,H,lambda1,lambda2,mu,eta){
    one_k <- matrix(1, nrow = 1,ncol= nrow(H))
    one_N <- matrix(1, nrow = 1,ncol = ncol(H))
    L1 = 0.5 * norm(X-W%*%H, type='F')^2
    L2 = lambda1 * sum(diag(t(W)%*%AS))
    L3 = lambda2 * sum(diag(t(W)%*%A0))
    L4 = mu/2 * norm(W,type="F")^2
    L5 = eta/2 * norm(one_k %*% H - one_N, type = "2")^2
    L= L1 - L2 - L3 + L4 + L5
    #print(c(L,L1,L2,L3,L4,L5))
    return(L)
}

compute_L1 <- function(X,rst){
    W = rst$W
    H = rst$H
    L1 = 0.5 * norm(X-W%*%H, type='F')^2
    return(L1)
}