set.seed(2020)
source('Opt/run_opt.R')

## full penalization
lambda1_vec <- lambda2_vec <- seq(0.3, 3, 0.1)
loss_first <- matrix(, nrow = length(lambda1_vec), ncol = length(lambda2_vec))
for (i in 10:length(lambda1_vec)) {
    for (j in 1:length(lambda2_vec)) {
        rst<-run(X,lambda1_vec[i],lambda2_vec[j],1,AS,A0,D,K,N, epsilon = 10^(-6))
        loss_first[i, j] <- 0.5 * norm(X - rst$W %*% rst$H, type='F')^2
        print(paste0('lambda1:', lambda1_vec[i]))
        print(paste0('lambda2:', lambda2_vec[j]))
    }
}

for (j in 9:length(lambda2_vec)) {
        rst<-run(X,lambda1_vec[i],lambda2_vec[j],1,AS,A0,D,K,N, epsilon = 10^(-6))
        loss_first[i, j] <- 0.5 * norm(X - rst$W %*% rst$H, type='F')^2
        print(paste0('lambda1:', lambda1_vec[i]))
        print(paste0('lambda2:', lambda2_vec[j]))
    }

rownames(loss_first) <- paste('lambda1', seq(0.3, 3, 0.1))
colnames(loss_first) <- paste('lambda2', seq(0.3, 3, 0.1))
heatmap(as.matrix(loss_first), Rowv = NA, Colv = NA, main = 'Loss with lambda1 and lambda2', scale = 'none')

write.table(loss_first, '/Users/mac/Desktop/Yale/Hongyu/PD/loss/loss_first_snr1.txt')

## without AS
loss_first_wo_AS <- vector(, length = length(lambda2_vec))
for (i in 1:length(lambda2_vec)) {
    rst<-run(X,0,lambda2_vec[i],1,AS,A0,D,K,N, epsilon = 10^(-6))
    loss_first_wo_AS[i] <- 0.5 * norm(X - rst$W %*% rst$H, type='F')^2
    print(paste0('lambda2:', lambda2_vec[i]))
}
write.table(loss_first_wo_AS, '/Users/mac/Desktop/Yale/Hongyu/PD/loss/loss_first_wo_AS.txt')

## without A0
loss_first_wo_A0 <- vector(, length = length(lambda1_vec))
for (i in 1:length(lambda1_vec)) {
    rst<-run(X,lambda1_vec[i],0,1,AS,A0,D,K,N, epsilon = 10^(-6))
    loss_first_wo_A0[i] <- 0.5 * norm(X - rst$W %*% rst$H, type='F')^2
    print(paste0('lambda1:', lambda1_vec[i]))
}
write.table(loss_first_wo_AS, '/Users/mac/Desktop/Yale/Hongyu/PD/loss/loss_first_wo_A0.txt')


infer <- function(truth, h){
    return(1*(h[truth]>0))
}
infer_max <- function(truth,h){
    return(1*(h[truth]==max(h)))
}

## best full penalization
lambda1_bst <- lambda1_vec[which(loss_first == min(loss_first, na.rm = T), arr.ind = TRUE)[1]]
lambda2_bst <- lambda2_vec[which(loss_first == min(loss_first, na.rm = T), arr.ind = TRUE)[2]]
rst_full <- run(X,lambda1_bst,lambda2_bst,1,AS,A0,D,K,N, epsilon = 10^(-6))

truth = label.output$label
cnt_full = 0
for(i in 1:length(truth)){
    cnt_full = cnt_full + infer(truth[i], rst_full$H[,i])
}

cnt_full_max = 0
for(i in 1:length(truth)){
    cnt_full_max = cnt_full_max + infer_max(truth[i], rst_full$H[,i])
}  
print(c(cnt_full,cnt_full_max))

## best without AS
rst_wo_AS <- run(X,0,lambda2_vec[which(loss_first_wo_AS == min(loss_first_wo_AS))],1,AS,A0,D,K,N, epsilon = 10^(-6))
cnt_wo_AS = 0
for(i in 1:length(truth)){
    cnt_wo_AS = cnt_wo_AS + infer(truth[i], rst_wo_AS$H[,i])
}

cnt_wo_AS_max = 0
for(i in 1:length(truth)){
    cnt_wo_AS_max = cnt_wo_AS_max + infer_max(truth[i], rst_wo_AS$H[,i])
}  
print(c(cnt_wo_AS,cnt_wo_AS_max))

## best without A0
rst_wo_A0 <- run(X,lambda1_vec[which(loss_first_wo_A0 == min(loss_first_wo_A0))],0,1,AS,A0,D,K,N, epsilon = 10^(-6))
cnt_wo_A0 = 0
for(i in 1:length(truth)){
    cnt_wo_A0 = cnt_wo_A0 + infer(truth[i], rst_wo_A0$H[,i])
}

cnt_wo_A0_max = 0
for(i in 1:length(truth)){
    cnt_wo_A0_max = cnt_wo_A0_max + infer_max(truth[i], rst_wo_A0$H[,i])
}  
print(c(cnt_wo_A0,cnt_wo_A0_max))

barplot(c(cnt_full / N, cnt_wo_AS / N, cnt_wo_A0 / N), main="True postive rate", names.arg=c("full", "without AS", "without A0"))
barplot(c(cnt_full_max / N, cnt_wo_AS_max / N, cnt_wo_A0_max / N), main="True postive rate max", names.arg=c("full", "without AS", "without A0"))
