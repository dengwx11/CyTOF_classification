set.seed(2020)
source("getPara.R")



runOptimalPara <- function(X,AS,A0,D,K,N, epsilon = 0.01,fixed_loop=1000,depth=3,
                            lambda1.on =TRUE, lambda2.on = TRUE){
    ## step 1 : get an initial eta by kkt condition
    print("step 1 : get an initial eta by kkt condition")
    if(lambda1.on){
        lambda_1.init = 1
    }else{
        lambda_1.init = 0
    }
    
    if(lambda2.on){
        lambda_2.init = 10
    }else{
        lambda_2.init = 0
    }

    rst = run(X,lambda_1.init,lambda_2.init,50,50,AS,A0,D,K,N, epsilon = epsilon,fixed_loop=fixed_loop)
    plot(as.vector(W),as.vector(rst$W))
    W.init = rst$W
    H.init = rst$H
    eta.kkt = getEta1(X,W.init,H.init)
    print(eta.kkt)

    ## step 2 : get the initial lambdas by screening
    print("step 2 : get the initial lambdas by screening")
    lambda_1.sequence = c(10^(-1),10^(0),10^(1),10^(2))
    lambda_2.sequence = c(2*10^(-1),2*10^(0),2*10^(1),2*10^(2))
    if(lambda2.on){
        rst = getLambda1(X,object = 2,lambda_2.sequence,lambda_1.init, 50,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
        lambda_1.init = rst$lambda_1.init1
        lambda_2.init = rst$lambda_2.init1
        print(c(lambda_1.init,lambda_2.init,rst$ARI))
    }
    if(lambda1.on){
        rst = getLambda1(X,object = 1,lambda_1.sequence,lambda_2.init, 50,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
        lambda_1.init = rst$lambda_1.init1
        lambda_2.init = rst$lambda_2.init1
        print(c(lambda_1.init,lambda_2.init,rst$ARI))
    }
    
    

    ## step 3 : get mu by W reliability
    print("step 3 : get mu by W reliability")
    mu.sequence = c(10^(-1),10^(0),10^(1),10^(2))
    rst =  getMu1(X, lambda_1.init,lambda_2.init,mu.sequence,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
    mu.init = rst$mu.init
    idx = rst$idx
    print(paste0("mu = ",mu.init))


    # ## step 4 : get an second initial eta by kkt condition
    # print("step 4 : get an second initial eta by kkt condition")
    # rst.eta = run(X,lambda_1.init2,lambda_2.init2,mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon, fixed_loop=fixed_loop)
    # plot(as.vector(W),as.vector(rst.eta$W))
    # W.ica = rst.eta$W
    # H.ica = rst.eta$H
    # eta.kkt = getEta1(X,W.ica,H.ica)
    # print(paste0("eta = ",eta.kkt))

    # ## step 5 : select the best lambda
    # print("step 5 : select the best lambda")
    # rst.lambda2  = getLambda2(X,object = 2,depth=3,lambda_2.init2,lambda_1.init2, mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
    # lambda_2.init3 = rst.lambda2$lambda
    # print(lambda_2.init3)
    # rst.lambda1 = getLambda2(X,object = 1,depth=3,lambda_1.init2,lambda_2.init3, mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
    # lambda_1.init3 = rst.lambda1$lambda
    # print(lambda_1.init3)


    ## step 6 : select the best mu
    print("step 6 : select the best mu")
    rst =  getMu2(X, depth = depth, lambda_1.init,lambda_2.init,mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
    mu.init = rst$mu
    print(paste0("mu = ",mu.init))

    ## step 7 : select the best lambda
    print("step 7 : select the best lambda")
    if(lambda2.on){
        rst = getLambda2(X,object = 2,depth=depth,lambda_2.init,lambda_1.init, mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
        lambda_2.init = rst$lambda
        print(lambda_2.init)
    }

    if(lambda1.on){
        rst = getLambda2(X,object = 1,depth=depth,lambda_1.init,lambda_2.init, mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
        lambda_1.init = rst$lambda
        print(lambda_1.init)
    }

    ARI = rst$ARI
    print(paste0("ARI after step 7 is ",ARI))
    
    ## step 8 : select the best mu and lambda one more time if API is too low
    print("step 8 : select the best mu and lambda one more time if ARI is too low")
    if(ARI<0.85){
        print("select the best eta one more time")
        rst =  getEta2(X, depth = depth, lambda_1.init,lambda_2.init,mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
        eta.kkt = rst$eta
        print(paste0("eta = ",eta.kkt))

        print("select the best mu one more time")
        rst =  getMu2(X, depth = depth, lambda_1.init,lambda_2.init,mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
        mu.init = rst$mu
        print(paste0("mu = ",mu.init))

        print("select the best lambda one more time")
        if(lambda2.on){
            rst = getLambda2(X,object = 2,depth=depth,lambda_2.init,lambda_1.init, mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
            lambda_2.init = rst$lambda
            print(lambda_2.init)
        }
        if(lambda1.on){
            rst = getLambda2(X,object = 1,depth=depth,lambda_1.init,lambda_2.init, mu.init,eta.kkt,AS,A0,D,K,N,epsilon = epsilon,fixed_loop=fixed_loop)
            lambda_1.init = rst$lambda
            print(lambda_1.init)
        }
        ARI = rst$ARI
        print(paste0("ARI after step 8 is ",ARI))
    }


    rst$para <- list()
    rst$para$lambda1 <- lambda_1.init
    rst$para$lambda2 <- lambda_2.init
    rst$para$mu <- mu.init
    rst$para$eta <- eta.kkt

    return(rst)
}




##################################


# rst<-run(X,lambda_1.init4,lambda_2.init4,mu.init2,eta.kkt,AS,A0,D,K,N, epsilon = .1,fixed_loop=1500)
# celltype_pred <- apply(rst$H, 2, predict)*2+10
# adjustedRandIndex(celltype_pred, seur$seurat_clusters)
# cor(as.vector(W),as.vector(rst$W))


# rst2<-run(X,1,lambda_2.init4,mu.init2,eta.kkt,AS,A0,D,K,N, epsilon = 10^(-3))
# celltype_pred <- apply(rst2$H, 2, predict)
# adjustedRandIndex(celltype_pred, seur$seurat_clusters)

# k=4
# rst = rst.seq[[k]]
# plot(as.vector(W),as.vector(rst$W))
# hist(as.vector(rst$W),30)
# max(rst$W)

######### 

# dist.matrix <- dist(x = seur@reductions$umap@cell.embeddings)
# clusters <- seur$seurat_clusters
# sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
# median(sil[, 3])
# hist(sil[, 3],30)