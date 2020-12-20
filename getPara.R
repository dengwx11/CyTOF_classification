set.seed(2020)
source('run_opt.R')
source('preparation.R')
library(ica)
library(mclust)
library(cluster, quietly = TRUE)

### eta initialization based on KKT conditions
getICAInit <- function(X,K){
    imod <- icafast(t(X),K)
    W.ica <- imod$M
    H.ica <- t(imod$S)

    rst = list()
    rst$W = W.ica
    rst$H = H.ica
    return(rst)
}    

## eta: kkt
getEta1 <- function(X,W.init,H.init){

    one_k <- matrix(1, nrow = 1,ncol= K)
    one_N <- matrix(1, nrow = 1,ncol = N)
    eta.kkt <- (t(W.init)%*%W.init%*%H.init - t(W.init)%*%X)/(t(one_k)%*%one_k%*%H.init-t(one_k)%*%one_N)
    eta.kkt <- median(abs(eta.kkt))

    return(eta.kkt)
}

getEta2 <- function(X, depth = 2, lambda_1,lambda_2,mu,eta.init,AS,A0,D,K,N,epsilon,fixed_loop){


    eta.prev <- eta.init
    global.max = eta.prev
    global.min = eta.prev
    rst.prev<-run(X,lambda_1,lambda_2,mu,eta.prev,AS,A0,D,K,N, epsilon,fixed_loop)

    celltype_pred <- apply(rst.prev$H, 2, predict)
    ARI.prev <- adjustedRandIndex(celltype_pred, seur$seurat_clusters)
    Wcorr.prev = cor(as.vector(W),as.vector(rst.prev$W))

    mean.tmp <- mean(rst.prev$W)
    median.tmp <- median(rst.prev$W)
    diff.prev <- abs(mean.X-mean.tmp) + abs(median.X-median.tmp)

    ans=list()
    ans$rst = rst.prev

    for(layer in 1:depth){
        if(eta.prev ==  global.min) eta.left = eta.prev/2
        else eta.left = (eta.prev + eta.prev.left)/2
        if(eta.prev == global.max) eta.right = eta.prev*2
        else eta.right = (eta.prev + eta.prev.right)/2

        candidates = c(eta.left,eta.prev,eta.right)


        local.min = min(candidates)
        local.max = max(candidates)
        if(local.min < global.min) global.min = local.min
        if(local.max > global.max) global.max = local.max

        
        rst.right <- run(X,lambda_1,lambda_2,mu,eta.right,AS,A0,D,K,N,epsilon,fixed_loop)
        rst.left <- run(X,lambda_1,lambda_2,mu,eta.left,AS,A0,D,K,N,epsilon,fixed_loop)  

        celltype_pred.left <- apply(rst.left$H, 2, predict)*2+10
        ARI.left <- adjustedRandIndex(celltype_pred.left, seur$seurat_clusters)
        celltype_pred.right <- apply(rst.right$H, 2, predict)*2+10
        ARI.right <- adjustedRandIndex(celltype_pred.right, seur$seurat_clusters)

        Wcorr.right = cor(as.vector(W),as.vector(rst.right$W))
        Wcorr.left = cor(as.vector(W),as.vector(rst.left$W))

        mean.left <- mean(rst.left$W)
        median.left <- median(rst.left$W)
        diff.left <- abs(mean.X-mean.left) + abs(median.X-median.left)
        mean.right <- mean(rst.right$W)
        median.right <- median(rst.right$W)
        diff.right <- abs(mean.X-mean.right) + abs(median.X-median.right)

        print(paste0("layer = ", layer))
        print(candidates)
        print(c(ARI.left,ARI.prev,ARI.right))
        print(c(Wcorr.left,Wcorr.prev,Wcorr.right))
        print(c(diff.left,diff.prev,diff.right))

        # min.idx = which.min(c(diff.left,diff.prev,diff.right))
        # if(min.idx == 2){   return(mu.prev) }
        # if(min.idx==1){
        #     mu.prev.right = mu.prev
        #     mu.prev = mu.left
        #     ARI.prev = ARI.left
        #     Wcorr.prev = Wcorr.left
        #     diff.prev = diff.left
        # }else if(min.idx == 3){
        #     mu.prev.left = mu.prev
        #     mu.prev = mu.right
        #     ARI.prev = ARI.right
        #     Wcorr.prev = Wcorr.right
        #     diff.prev = diff.right
        # }

        max.idx = which.max(c(ARI.left,ARI.prev,ARI.right))
        if(max.idx == 2){   
            eta.prev.right = eta.right
            eta.prev.left = eta.left
            ans$rst = rst.prev
        }
        if(max.idx==1){
            eta.prev.right = eta.prev
            eta.prev = eta.left
            ARI.prev = ARI.left
            Wcorr.prev = Wcorr.left
            diff.prev = diff.left
            ans$rst = rst.left
        }else if(max.idx == 3){
            eta.prev.left = eta.prev
            eta.prev = eta.right
            ARI.prev = ARI.right
            Wcorr.prev = Wcorr.right
            diff.prev = diff.right
            ans$rst = rst.right
        }



    }

    ans$eta = eta.prev
    ans$ARI = ARI.prev
    


    return(ans)


    
}



# ### the scale between AS and A_0: A_0 = g * AS
# g = sum(AS*A0)/(norm(AS,type='2')^2)

# ### \lambda_1 initialization
# mu.init = 100
# lambda_1.kkt <- (W.ica%*%(H.ica%*%t(H.ica)+mu.init*diag(K))-X%*%t(H.ica))/2/AS
# lambda_1.kkt <- min(abs(lambda_1.kkt))
# lambda_2.kkt <- lambda_1.kkt/g
# print(c(lambda_1.kkt,lambda_2.kkt))


## lambda: ARI with louvian cluster
getLambda1 <- function(X,object = 1,lambda_object.sequence, theotherLambda, mu.init,eta,AS,A0,D,K,N,epsilon,fixed_loop){
    rst.seq <- list()
    L1.seq <- c()
    k=1

    L1.init1 <- 10^(8)
    ARI <- -10
    if(object==1){
        for(lambda_1 in lambda_object.sequence){

            rst.seq[[k]]<-run(X,lambda_1,theotherLambda,mu.init,eta,AS,A0,D,K,N, epsilon,fixed_loop)
            L1.tmp <- compute_L1(X,rst.seq[[k]])
            L1.seq <- c(L1.seq, L1.tmp)

            celltype_pred <- apply(rst.seq[[k]]$H, 2, predict)
            ARI.tmp <- adjustedRandIndex(celltype_pred, seur$seurat_clusters)

            print(paste0("lambda = ", lambda_1,", ARI = ",ARI.tmp))

            if(ARI.tmp > ARI){
                ARI = ARI.tmp
                lambda_1.init1 = lambda_1
                lambda_2.init1 = theotherLambda 
            }

            # if(L1.tmp < L1.init1){
            #     L1.init1 = L1.tmp
            #     lambda_1.init1 = lambda_1
            #     lambda_2.init1 = theotherLambda
            # } 
            k=k+1
        }
    }else{
           
        for(lambda_2 in lambda_object.sequence){
            rst.seq[[k]]<-run(X,theotherLambda,lambda_2,mu.init,eta,AS,A0,D,K,N, epsilon,fixed_loop)
            L1.tmp <- compute_L1(X,rst.seq[[k]])
            L1.seq <- c(L1.seq, L1.tmp)

            celltype_pred <- apply(rst.seq[[k]]$H, 2, predict)
            ARI.tmp <- adjustedRandIndex(celltype_pred, seur$seurat_clusters)

            print(paste0("lambda = ", lambda_2,", ARI = ",ARI.tmp))

            if(ARI.tmp > ARI){
                ARI = ARI.tmp
                lambda_1.init1 = theotherLambda
                lambda_2.init1 = lambda_2 
            }

            # if(L1.tmp < L1.init1){
            #     L1.init1 = L1.tmp
            #     lambda_1.init1 = theotherLambda
            #     lambda_2.init1 = lambda_2
            # } 
            k=k+1
        
    }

    }


    ans <- list()
    ans$lambda_1.init1 = lambda_1.init1
    ans$lambda_2.init1 = lambda_2.init1
    ans$ARI = ARI 
    ans$rst.seq = rst.seq
    return(ans)
}

getLambda2 <- function(X,object = 1,depth = 3,lambda_init, theotherLambda, mu,eta,AS,A0,D,K,N,epsilon,fixed_loop){


    lambda.prev = lambda_init
    global.max = lambda.prev
    global.min = lambda.prev
    if(object == 1){
        rst.prev<-run(X,lambda.prev,theotherLambda,mu,eta,AS,A0,D,K,N, epsilon,fixed_loop)
    }else{
        rst.prev<-run(X,theotherLambda,lambda.prev,mu,eta,AS,A0,D,K,N, epsilon ,fixed_loop)
    }
    celltype_pred <- apply(rst.prev$H, 2, predict)
    ARI.prev <- adjustedRandIndex(celltype_pred, seur$seurat_clusters)
    Wcorr.prev = cor(as.vector(W),as.vector(rst.prev$W))

    ans = list()
    

    for(layer in 1:depth){
        if(lambda.prev ==  global.min) lambda.left = lambda.prev/2
        else lambda.left = (lambda.prev + lambda.prev.left)/2
        if(lambda.prev == global.max) lambda.right = lambda.prev*2
        else lambda.right = (lambda.prev + lambda.prev.right)/2

        candidates = c(lambda.left,lambda.prev,lambda.right)


        local.min = min(candidates)
        local.max = max(candidates)
        if(local.min < global.min) global.min = local.min
        if(local.max > global.max) global.max = local.max

        if(object == 1){
            rst.right <- run(X,lambda.right,theotherLambda,mu,eta,AS,A0,D,K,N,epsilon,fixed_loop)
            rst.left <- run(X,lambda.left,theotherLambda,mu,eta,AS,A0,D,K,N,epsilon=epsilon,fixed_loop=fixed_loop)  
        }else{
            rst.right <- run(X,theotherLambda,lambda.right,mu,eta,AS,A0,D,K,N,epsilon,fixed_loop)
            rst.left <- run(X,theotherLambda,lambda.left,mu,eta,AS,A0,D,K,N,epsilon,fixed_loop)
        }
        celltype_pred.left <- apply(rst.left$H, 2, predict)*2+10
        ARI.left <- adjustedRandIndex(celltype_pred.left, seur$seurat_clusters)
        celltype_pred.right <- apply(rst.right$H, 2, predict)*2+10
        ARI.right <- adjustedRandIndex(celltype_pred.right, seur$seurat_clusters)

        Wcorr.right = cor(as.vector(W),as.vector(rst.right$W))
        Wcorr.left = cor(as.vector(W),as.vector(rst.left$W))

        print(paste0("layer = ", layer))
        print(candidates)
        print(c(ARI.left,ARI.prev,ARI.right))
        print(c(Wcorr.left,Wcorr.prev,Wcorr.right))

        max.idx = which.max(c(ARI.left,ARI.prev,ARI.right))
        if(max.idx == 2){ 
            lambda.prev.right = lambda.right
            lambda.prev.left = lambda.left
            ans$rst = rst.prev
            }
        if(max.idx==1){
            lambda.prev.right = lambda.prev
            lambda.prev = lambda.left
            ARI.prev = ARI.left
            Wcorr.prev = Wcorr.left
            ans$rst = rst.left
        }else if(max.idx == 3){
            lambda.prev.left = lambda.prev
            lambda.prev = lambda.right
            ARI.prev = ARI.right
            Wcorr.prev = Wcorr.right
            ans$rst = rst.right
        }



    }

    ans$lambda = lambda.prev
    ans$ARI = ARI.prev

    return(ans)


    
}


## similarities of median and mean between X and W
getMu1 <- function(X, lambda_1,lambda_2,mu.sequence, eta,AS,A0,D,K,N,epsilon,fixed_loop){
    rst.seq <- list()
    L1.seq <- c()
    sparsity.seq <- c()
    diff <- 10^8
    mu.init = mu.sequence[1]
    k=1

    mean.X <- mean(X)
    median.X <- median(X)
    idx = 1
    for(mu in mu.sequence){
        rst.seq[[k]]<-run(X,lambda_1,lambda_2,mu,eta,AS,A0,D,K,N, epsilon ,fixed_loop)

        L1.tmp <- compute_L1(X,rst.seq[[k]])
        L1.seq <- c(L1.seq, L1.tmp)

        mean.tmp <- mean(rst.seq[[k]]$W)
        median.tmp <- median(rst.seq[[k]]$W)
        diff.tmp <- abs(mean.X-mean.tmp) + abs(median.X-median.tmp)
        print(diff.tmp)
        if(diff.tmp < diff){
            diff = diff.tmp
            mu.init = mu
            idx = k

        }


        k = k+1
    }

    

    ans = list()
    ans$rst.seq = rst.seq
    ans$mu.init = mu.init
    ans$idx = idx
    return(ans)
    
}



getMu2 <- function(X, depth = 2, lambda_1,lambda_2,mu.init,eta,AS,A0,D,K,N,epsilon ,fixed_loop){


    mu.prev <- mu.init
    global.max = mu.prev
    global.min = mu.prev
    rst.prev<-run(X,lambda_1,lambda_2,mu.prev,eta,AS,A0,D,K,N, epsilon,fixed_loop)

    celltype_pred <- apply(rst.prev$H, 2, predict)
    ARI.prev <- adjustedRandIndex(celltype_pred, seur$seurat_clusters)
    Wcorr.prev = cor(as.vector(W),as.vector(rst.prev$W))

    mean.tmp <- mean(rst.prev$W)
    median.tmp <- median(rst.prev$W)
    mean.X <- mean(X)
    median.X <- median(X)
    diff.prev <- abs(mean.X-mean.tmp) + abs(median.X-median.tmp)

    ans=list()
    ans$rst = rst.prev

    for(layer in 1:depth){
        if(mu.prev ==  global.min) mu.left = mu.prev/2
        else mu.left = (mu.prev + mu.prev.left)/2
        if(mu.prev == global.max) mu.right = mu.prev*2
        else mu.right = (mu.prev + mu.prev.right)/2

        candidates = c(mu.left,mu.prev,mu.right)


        local.min = min(candidates)
        local.max = max(candidates)
        if(local.min < global.min) global.min = local.min
        if(local.max > global.max) global.max = local.max

        
        rst.right <- run(X,lambda_1,lambda_2,mu.right,eta,AS,A0,D,K,N,epsilon,fixed_loop)
        rst.left <- run(X,lambda_1,lambda_2,mu.left,eta,AS,A0,D,K,N,epsilon,fixed_loop)  

        celltype_pred.left <- apply(rst.left$H, 2, predict)*2+10
        ARI.left <- adjustedRandIndex(celltype_pred.left, seur$seurat_clusters)
        celltype_pred.right <- apply(rst.right$H, 2, predict)*2+10
        ARI.right <- adjustedRandIndex(celltype_pred.right, seur$seurat_clusters)

        Wcorr.right = cor(as.vector(W),as.vector(rst.right$W))
        Wcorr.left = cor(as.vector(W),as.vector(rst.left$W))

        mean.left <- mean(rst.left$W)
        median.left <- median(rst.left$W)
        diff.left <- abs(mean.X-mean.left) + abs(median.X-median.left)
        mean.right <- mean(rst.right$W)
        median.right <- median(rst.right$W)
        diff.right <- abs(mean.X-mean.right) + abs(median.X-median.right)

        print(paste0("layer = ", layer))
        print(candidates)
        print(c(ARI.left,ARI.prev,ARI.right))
        print(c(Wcorr.left,Wcorr.prev,Wcorr.right))
        print(c(diff.left,diff.prev,diff.right))

        # min.idx = which.min(c(diff.left,diff.prev,diff.right))
        # if(min.idx == 2){   return(mu.prev) }
        # if(min.idx==1){
        #     mu.prev.right = mu.prev
        #     mu.prev = mu.left
        #     ARI.prev = ARI.left
        #     Wcorr.prev = Wcorr.left
        #     diff.prev = diff.left
        # }else if(min.idx == 3){
        #     mu.prev.left = mu.prev
        #     mu.prev = mu.right
        #     ARI.prev = ARI.right
        #     Wcorr.prev = Wcorr.right
        #     diff.prev = diff.right
        # }

        max.idx = which.max(c(ARI.left,ARI.prev,ARI.right))
        if(max.idx == 2){   
            mu.prev.right = mu.right
            mu.prev.left = mu.left
            ans$rst = rst.prev
        }
        if(max.idx==1){
            mu.prev.right = mu.prev
            mu.prev = mu.left
            ARI.prev = ARI.left
            Wcorr.prev = Wcorr.left
            diff.prev = diff.left
            ans$rst = rst.left
        }else if(max.idx == 3){
            mu.prev.left = mu.prev
            mu.prev = mu.right
            ARI.prev = ARI.right
            Wcorr.prev = Wcorr.right
            diff.prev = diff.right
            ans$rst = rst.right
        }



    }

    ans$mu = mu.prev
    ans$ARI = ARI.prev
    


    return(ans)


    
}
