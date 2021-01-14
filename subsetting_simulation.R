options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

rst.directory <- as.character(args[1])
truth.directory <- as.character(args[2])

rst <- readRDS(rst.directory)
truth <- readRDS(truth.directory)

source("preparation.R")

celltype_pred <- apply(rst$H, 2, predict)
H_est_subsetting <- data.frame(rst$H)
H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))

cutoff.list <- c(1:100)/100                          
df <- data.frame(cutoff = cutoff.list)
df$acc <- df$cos <- df$nmi <- df$ari <- df$sil <- 0                         
                          
for(j in c(1:length(cutoff.list))){
    cnt_max = 0
    N_subsetting = 0
    idx.subsetting = 0
        for(i in 1:length(truth_ct)){
            if(max(H_est_subsetting[,i])>cutoff){
                idx.subsetting = c(idx.subsetting,i)
                N_subsetting = N_subsetting+1
                cnt_max = cnt_max + 1*(truth[i]==celltype_pred[i])
            }

        }  
    #cnt_max
    df$acc[j]<-cnt_max/N_subsetting
}