options(stringsAsFactors=F)
args = commandArgs(trailingOnly=TRUE)

rst.directory <- as.character(args[1])
truth.directory <- as.character(args[2])
seur.directory <- as.character(args[3])

rst <- readRDS(rst.directory)
truth <- readRDS(truth.directory)
seur <- readRDS(seur.directory)

source("preparation.R")
library(mltools)
library(data.table)
library(lsa)
library(aricode)
library(kBET)
library(cowplot)

celltype_pred <- apply(rst$H, 2, predict)
H_est_subsetting <- data.frame(rst$H)
H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))

cutoff.list <- c(1:100)/100                          
df <- data.frame(cutoff = cutoff.list)
df$acc <- df$cos <- df$nmi <- df$ari <- df$sil <- 0   
df$cutoff <- cutoff.list                    
                          
for(j in c(1:length(cutoff.list))){
    ## cnt_max
    cnt_max = 0
    N_subsetting = 0
    idx.subsetting = 0
        for(i in 1:length(truth)){
            if(max(H_est_subsetting[,i])>cutoff.list[j]){
                idx.subsetting = c(idx.subsetting,i)
                N_subsetting = N_subsetting+1
                cnt_max = cnt_max + 1*(truth[i]==celltype_pred[i])
            }

        }
    if(N_subsetting == 0) {
        df$acc[j] <- 0
    } else {
       df$acc[j]<-cnt_max/N_subsetting
    } 

    ## Cosine similarity
    truth_sub <- truth[idx.subsetting]
    truth_fact_sub <- as.factor(truth_sub)
    truth_fact_sub <- factor(truth_fact_sub, levels = as.character(1:K))
    truth_onehot_sub <- as.data.frame(t(one_hot(as.data.table(truth_fact_sub))))
    cos_sim <- mean(mapply(cosine, truth_onehot_sub, as.data.frame(H_est_subsetting[,idx.subsetting])))
    df$cos[j] <- cos_sim

    ## NMI
    celltype_pred_sub <- celltype_pred[idx.subsetting]
    celltype_pred_fact_sub <- as.factor(celltype_pred_sub)
    celltype_pred_fact_sub <- factor(celltype_pred_fact_sub, levels = as.character(1:K))
    pred_onehot_sub <- as.data.frame(t(one_hot(as.data.table(celltype_pred_fact_sub))))
    nmi <- mean(mapply(NMI, truth_onehot_sub, as.data.frame(pred_onehot_sub)))
    df$nmi[j] <- nmi

    ## ARI
    ari <- adjustedRandIndex(celltype_pred_sub, truth_sub)
    df$ari[j] <- ari

    ## Silhouette
    pca.data <- list()
    pca.data$x <- seur@reductions$pca@cell.embeddings[idx.subsetting,]
    sil <- batch_sil(pca.data, celltype_pred_sub)
    df$sil[j] <- sil
}

subset_6 <- melt(df, id.vars = 'cutoff', variable.name = 'Metrics')
write.table(subset_6, '/Users/mac/Desktop/Yale/Hongyu/CyTOF/subset_6.txt')
# subset_1 <- read.table('/Users/mac/Desktop/Yale/Hongyu/CyTOF/subset_1.txt')
ggplot(subset_6, aes(cutoff,value)) + 
geom_point(aes(colour = Metrics), shape = 4) +
labs(y = 'Annotation Metric', x = 'Cutoff', title = 'Scenario 1') + 
geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
# ylim(0.8, 1) +
theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))

pars_lst <- list(subset_1, subset_2, subset_3, subset_4, subset_5, subset_6)
plotlist <- list()
for(i in 1:length(pars_lst)){
  plotlist[[i]] <- ggplot(pars_lst[[i]], aes(cutoff,value)) + 
    geom_point(aes(colour = Metrics), shape = 4) +
    labs(y = 'Annotation Metric', x = 'Cutoff', title = paste('Scenario', i)) + 
    geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
    theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
}
pll <- plot_grid(plotlist = plotlist, ncol=3)
ggsave('/Users/mac/Desktop/Yale/Hongyu/CyTOF/Plots/subset_simu.png', pll, width = 15, height = 8)
    
