## Run on server

options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

which_lambda <- as.character(args[1]) ## 'lambda1'/'lambda2'
lambda <- as.character(args[2])
senario <- as.character(args[3]) ## 1-6

rst.para <- readRDS(paste0('/home/bz234/project/Results/CyTOF/simulation/rst_para/rst_para_', senario, '.rds'))

source("./Simulation.R")
source('./run_para.R')
library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)

if(senario == '1') {
    set.seed(2020)
    K = 5 # cell types number ## K could be larger
    D = 10 # surface markers number
    N = 2000 # ADT/CyTOF cell number ## N could be larger
    G = 100 # RNA gene number
    pi_ber1 = 0.55
    pi_ber2 = 0.9

    ## simulation parameters:
    iteration = 100
    big_w_mean = 2
    big_tau_w = 10
    small_w_mean = 0.5
    small_tau_w = 10
    p.0 = 0.1
    q.0 = 0.2
    p.neg1 = 0.05
    q.neg1 = 0.1
    mean_var_ratio = 10
    corr = 0.5
    prob_k = c(1,2,2,3,3)
}
if(senario == '2') {
    set.seed(2020)
    K = 8 # cell types number ## K could be larger
    D = 10 # surface markers number
    N = 2000 # ADT/CyTOF cell number ## N could be larger
    G = 100 # RNA gene number
    pi_ber1 = 0.55
    pi_ber2 = 0.9

    ## simulation parameters:
    iteration = 100
    big_w_mean = 2
    big_tau_w = 8
    small_w_mean = 1
    small_tau_w = 8
    p.0 = 0.15
    q.0 = 0.3
    p.neg1 = 0.1
    q.neg1 = 0.15
    mean_var_ratio = 5
    corr = 0.2
    #prob_k = c(1,2,2,3,3)
    prob_k = c(3,1,2,3,3,1,2,1)
}
if(senario == '2') {
    set.seed(2020)
    K = 8 # cell types number ## K could be larger
    D = 10 # surface markers number
    N = 2000 # ADT/CyTOF cell number ## N could be larger
    G = 100 # RNA gene number
    pi_ber1 = 0.55
    pi_ber2 = 0.9

    ## simulation parameters:
    iteration = 100
    big_w_mean = 2
    big_tau_w = 8
    small_w_mean = 1
    small_tau_w = 8
    p.0 = 0.15
    q.0 = 0.3
    p.neg1 = 0.1
    q.neg1 = 0.15
    mean_var_ratio = 5
    corr = 0.2
    #prob_k = c(1,2,2,3,3)
    prob_k = c(3,1,2,3,3,1,2,1)
}
if(senario == '3') {
    set.seed(2020)
    K = 8 # cell types number ## K could be larger
    D = 10 # surface markers number
    N = 2000 # ADT/CyTOF cell number ## N could be larger
    G = 100 # RNA gene number
    pi_ber1 = 0.55
    pi_ber2 = 0.9

    ## simulation parameters:
    iteration = 100
    big_w_mean = 2
    big_tau_w = 8
    small_w_mean = 1
    small_tau_w = 8
    p.0 = 0.15
    q.0 = 0.3
    p.neg1 = 0.1
    q.neg1 = 0.15
    mean_var_ratio = 3
    corr = 0.3
    #prob_k = c(1,2,2,3,3)
    prob_k = c(3,1,2,3,3,1,2,1)
}
if(senario == '4') {
    set.seed(2020)
    K = 8 # cell types number ## K could be larger
    D = 10 # surface markers number
    N = 2000 # ADT/CyTOF cell number ## N could be larger
    G = 100 # RNA gene number
    pi_ber1 = 0.55
    pi_ber2 = 0.9

    ## simulation parameters:
    iteration = 100
    big_w_mean = 2
    big_tau_w = 7
    small_w_mean = 1
    small_tau_w = 7
    p.0 = 0.15
    q.0 = 0.3
    p.neg1 = 0.1
    q.neg1 = 0.15
    mean_var_ratio = 3
    corr = 0.2
    #prob_k = c(1,2,2,3,3)
    prob_k = c(3,1,2,3,3,1,2,1)
}
if(senario == '5') {
    set.seed(2020)
    K = 8 # cell types number ## K could be larger
    D = 10 # surface markers number
    N = 2000 # ADT/CyTOF cell number ## N could be larger
    G = 100 # RNA gene number
    pi_ber1 = 0.55
    pi_ber2 = 0.9

    ## simulation parameters:
    iteration = 100
    big_w_mean = 2
    big_tau_w = 7
    small_w_mean = 1
    small_tau_w = 7
    p.0 = 0.2
    q.0 = 0.4
    p.neg1 = 0.1
    q.neg1 = 0.2
    mean_var_ratio = 3
    corr = 0.2
    #prob_k = c(1,2,2,3,3)
    prob_k = c(3,1,2,3,3,1,2,1)
}
if(senario == '6') {
    set.seed(2020)
    K = 8 # cell types number ## K could be larger
    D = 10 # surface markers number
    N = 2000 # ADT/CyTOF cell number ## N could be larger
    G = 100 # RNA gene number
    pi_ber1 = 0.55
    pi_ber2 = 0.9

    ## simulation parameters:
    iteration = 100
    big_w_mean = 2
    big_tau_w = 5
    small_w_mean = 1
    small_tau_w = 5
    p.0 = 0.2
    q.0 = 0.4
    p.neg1 = 0.1
    q.neg1 = 0.2
    mean_var_ratio = 2
    corr = 0.2
    #prob_k = c(1,2,2,3,3)
    prob_k = c(3,1,2,3,3,1,2,1)
}

gamma <- simulate_gamma(iteration, K, D)
A0 <- gamma$gamma
W <- simulate_w(D,K,A0, p.0 = p.0, q.0=q.0, p.neg1= p.neg1, q.neg1 = q.neg1, 
                big_w_mean=big_w_mean, big_tau_w = big_tau_w, 
                small_w_mean=small_w_mean, small_tau_w = small_tau_w)
AS <- simulate_AS(D,K,W,corr=corr)
label.output <- simulate_label(D,N,K,prob_k)
X <- simulate_x(D,N,W,label.output$label, mean_var_ratio = mean_var_ratio)
true.H <- label.output$H
rownames(X) <- c(1:D)
colnames(X) <- c(1:N)
seur <- CreateSeuratObject(X)
seur@meta.data$label <- label.output$label[1,]
seur <- NormalizeData(seur) 
seur@assays$RNA@data = seur@assays$RNA@counts
seur <- FindVariableFeatures(seur, slot='counts', selection.method = "vst", nfeatures = 10)
seur <- ScaleData(seur)
seur <- RunPCA(seur,verbose = TRUE,features = rownames(seur))
seur <- FindNeighbors(seur, dims = 1:6)
seur <- FindClusters(seur)
# seur <- RunUMAP(seur, dims = 1:5)

if(which_lambda == 'lambda1') {
    rst <- run(X,as.numeric(lambda),rst.para$para$lambda2,rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
}
if(which_lambda == 'lambda2') {
    rst <- run(X,rst.para$para$lambda1,as.numeric(lambda),rst.para$para$mu,rst.para$para$eta,
            AS,A0,D,K,N, epsilon = 10^(-3),fixed_loop=2000)
}
truth <- label.output$label
celltype_pred <- apply(rst$H, 2, predict)
truth_onehot <- as.data.frame(t(one_hot(as.data.table(as.factor(truth)))))
celltype_pred_fact <- as.factor(celltype_pred)
celltype_pred_fact <- factor(celltype_pred, levels = as.character(1:K))
pred_onehot <- as.data.frame(t(one_hot(as.data.table(celltype_pred_fact))))

cnt_max <- 0
for(i in 1:length(truth)){
    cnt_max = cnt_max + infer_max(truth[i], rst$H[,i])
}  

## Accuracy
accu <- cnt_max / N
## Cosine similarity
cos_sim <- mean(mapply(cosine, truth_onehot, as.data.frame(rst$H)))
## ARI
ari <- adjustedRandIndex(celltype_pred, truth)
## NMI
nmi <- mean(mapply(NMI, truth_onehot, as.data.frame(pred_onehot)))

results <- c(accu, cos_sim, ari, nmi)
names(results) <- c('accu', 'cos_sim', 'ari', 'nmi')
write.table(results, paste0('/home/bz234/project/Results/CyTOF/simulation/benchmark/s_', senario, '_', which_lambda, '_', lambda, '.txt'))