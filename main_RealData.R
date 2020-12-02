set.seed(2020)
source('run_opt.R')

## Input
A <- as.matrix(read.csv('../write/matrix_A/inhouse_w_hvg_covid_s_35_a.csv',sep=' '))
S <- readRDS('./write/matrix_S/hvgs_inhouse_w_covid_s.rds')
X <- as.matrix(read.table('./write/matrix_X/X_inhouse.txt'))
A0 <- read.csv('./write/matrix_A0/A0_inhouse.csv')
truth_ct = read.table('./write/true_label/celltype_inhouse.txt', sep='\t',header=T)
colnames(A0)=c("Gene","CD4+T cell","CD8+T cell", "naïve B cells","NK cell", "classical monocytes", "non-classical monocytes")
rownames(A0) = A0$Gene
A0 <- as.matrix(A0[-7,-1])
S <- as.matrix(S[,colnames(A0)])
AS <- t(A) %*% S
X <- X[-7,truth_ct$cellID]


## Other para
lambda1 = 0.5
lambda2 = 0.5
mu = 1
D = nrow(X)
K = ncol(AS)
N = ncol(X)

## 
rst<-run(X,2,0,2,AS,A0,D,K,N,epsilon=0.01)

##
H_est <- data.frame(rst$H)
H_est <- cbind(matrix(colnames(A0),ncol=1),H_est)


get_truth <- function(true_ct){
    return(which(as.character(H_est[,1]) == true_ct))
}
infer_max <- function(truth,h){
    return(1*(h[truth]==max(h)))
}

truth = sapply(as.character(truth_ct[,2]), get_truth)
cnt_max = 0
for(i in 1:length(truth)){
    cnt_max = cnt_max + infer_max(truth[i], rst$H[,i])
}  
cnt_max