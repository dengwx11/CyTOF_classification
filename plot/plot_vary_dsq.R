library(reshape2)
library(ggplot2)
library(cowplot)

corr <- seq(.1, 1, 0.005)
mean_var_ratio <- seq(1, 10, 0.05)
K <- seq(3, 15)
# D <- seq(5, 50)

### Write dsq file
writefunc <- c()
for(i in corr) {
    writefunc <- rbind(writefunc, paste('module load miniconda; source activate r_env; cd /home/bz234/github/CyTOF_classification/; Rscript --vanilla vary_corr_dsq.R', i))
} 
for(i in mean_var_ratio) {
    writefunc <- rbind(writefunc, paste('module load miniconda; source activate r_env; cd /home/bz234/github/CyTOF_classification/; Rscript --vanilla vary_mvr_dsq.R', i))
} 
for(i in K) {
    writefunc <- rbind(writefunc, paste('module load miniconda; source activate r_env; cd /home/bz234/github/CyTOF_classification/; Rscript --vanilla vary_k_dsq.R', i))
} 
# for(i in D) {
#     writefunc <- rbind(writefunc, paste('module load miniconda; source activate r_env; cd /home/bz234/github/CyTOF_classification/; Rscript --vanilla vary_d_dsq.R', i))
# }
# write.table(writefunc, file ='/home/bz234/Scripts/bash/CyTOF/dsq/vary_corr_mvr_k_d.txt', row.names = FALSE, col.names = FALSE, quote=FALSE)

### Read results and plot
files_mvr <- list()
for(i in 1:length(mean_var_ratio)) {
    files_mvr[[i]] <- paste0('/home/bz234/project/Results/CyTOF/simulation/vary_mvr/s2_mvr_', mean_var_ratio[i], '.txt')
}
files_corr <- list()
for(i in 1:length(corr)) {
    files_corr[[i]] <- paste0('/home/bz234/project/Results/CyTOF/simulation/vary_corr/s2_corr_', corr[i], '.txt')
}
files_k <- list()
for(i in 1:length(K)) {
    files_k[[i]] <- paste0('/home/bz234/project/Results/CyTOF/simulation/vary_k/s2_k_', K[i], '.txt')
}

## Read in files
read_error <- function (file) {
  return(tryCatch(read.table(file, header = TRUE), error=function(e) NULL))
}

table_lst_mvr <- lapply(files_mvr, read_error)
table_lst_corr <- lapply(files_corr, read_error)
table_lst_k <- lapply(files_k, read_error)

## Combine results
results_lst_mvr <- matrix(, nrow = length(files_mvr), ncol = 5)
colnames(results_lst_mvr) <- c('mean_var_ratio', rownames(table_lst_mvr[[1]]))
results_lst_mvr[,1] <- mean_var_ratio
for(i in 1:length(files_mvr)) {
    if(is.null(table_lst_mvr[[i]]$x)){
        results_lst_mvr[i,2:5] <- NA
    } else{
    results_lst_mvr[i,2:5] <- table_lst_mvr[[i]]$x
    }
}
results_lst_corr <- matrix(, nrow = length(files_corr), ncol = 5)
colnames(results_lst_corr) <- c('corr', rownames(table_lst_corr[[1]]))
results_lst_corr[,1] <- corr
for(i in 1:length(files_corr)) {
    if(is.null(table_lst_corr[[i]]$x)){
        results_lst_corr[i,2:5] <- NA
    } else{
    results_lst_corr[i,2:5] <- table_lst_corr[[i]]$x
    }
}
results_lst_k <- matrix(, nrow = length(files_k), ncol = 5)
colnames(results_lst_k) <- c('K', rownames(table_lst_k[[1]]))
results_lst_k[,1] <- K
for(i in 1:length(files_k)) {
    if(is.null(table_lst_k[[i]]$x)){
        results_lst_k[i,2:5] <- NA
    } else{
    results_lst_k[i,2:5] <- table_lst_k[[i]]$x
    }
}
## Manually add
results_lst_k[4,2:5] <- c(0.9610000, 0.9337315, 0.9369329, 0.9639048)

## Plot
results_mvr <- melt(as.data.frame(results_lst_mvr), id.vars = 'mean_var_ratio', variable.name = 'Metrics')
p_mvr <- ggplot(results_mvr, aes(mean_var_ratio,value)) + 
    geom_point(aes(colour = Metrics), shape = 4) +
    labs(y = 'Annotation Metric', x = 'Mean variance ratio', title = 'Vary mean variance ratio') + 
    geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
    theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
results_corr <- melt(as.data.frame(results_lst_corr), id.vars = 'corr', variable.name = 'Metrics')
p_corr <- ggplot(results_corr, aes(corr,value)) + 
    geom_point(aes(colour = Metrics), shape = 4) +
    labs(y = 'Annotation Metric', x = 'Corr', title = 'Vary correlation') + 
    geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
    theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
results_k <- melt(as.data.frame(results_lst_k), id.vars = 'K', variable.name = 'Metrics')
p_k <- ggplot(results_k, aes(K,value)) + 
    geom_point(aes(colour = Metrics), shape = 4) +
    labs(y = 'Annotation Metric', x = 'K', title = 'Vary K') + 
    geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
    theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
pll <- plot_grid(plotlist = list(p_mvr, p_corr, p_k), ncol=3)
ggsave('/home/bz234/project/Results/CyTOF/simulation/vary_mvr_corr_k_plot.png', pll, width = 15, height = 4)
 