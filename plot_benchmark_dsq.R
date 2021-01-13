library(reshape2)
library(ggplot2)
library(cowplot)

lambda2 <- seq(0, 60, 0.1)
lambda1 <- seq(0, 5, 0.01)

### Write dsq file
senario <- 1:6
writefunc <- c()
for(s in 1:length(senario)) {
    for(l in 1:length(lambda1)) {
        writefunc <- rbind(writefunc, paste('module load miniconda; source activate r_env; cd /home/bz234/github/CyTOF_classification/; Rscript --vanilla vary_lambda_dsq.R', 'lambda1', lambda1[l], senario[s]))
    } 
}
for(s in 1:length(senario)) {
    for(l in 1:length(lambda2)) {
        writefunc <- rbind(writefunc, paste('module load miniconda; source activate r_env; cd /home/bz234/github/CyTOF_classification/; Rscript --vanilla vary_lambda_dsq.R', 'lambda2', lambda2[l], senario[s]))
    } 
}
write.table(writefunc, file ='/home/bz234/Scripts/bash/CyTOF/dsq/benchmarking_all.txt', row.names = FALSE, col.names = FALSE, quote=FALSE)

### Read results and plot
files_1 <- list()
for(s in 1:6) {
    files_1[[s]] <- list()
    for(i in 1:length(lambda1)) {
        files_1[[s]][[i]] <- paste0('/home/bz234/project/Results/CyTOF/simulation/benchmark/s_', s, '_lambda1_', lambda1[i], '.txt')
    }
}
files_2 <- list()
for(s in 1:6) {
    files_2[[s]] <- list()
    for(i in 1:length(lambda2)) {
        files_2[[s]][[i]] <- paste0('/home/bz234/project/Results/CyTOF/simulation/benchmark/s_', s, '_lambda2_', lambda2[i], '.txt')
    }
}

## Read in files
read_error <- function (file) {
  return(tryCatch(read.table(file, header = TRUE), error=function(e) NULL))
}
table_lst_1 <- list()
for(s in 1:6) {
    table_lst_1[[s]] <- lapply(files_1[[s]], read_error)
}
table_lst_2 <- list()
for(s in 1:6) {
    table_lst_2[[s]] <- lapply(files_2[[s]], read_error)
}

## Combine results
results_lst_1 <- list()
for(s in 1:6) {
    results_lst_1[[s]] <- matrix(, nrow = length(files_1[[s]]), ncol = 5)
    colnames(results_lst_1[[s]]) <- c('lambda1', rownames(table_lst_1[[s]][[1]]))
    results_lst_1[[s]][,1] <- lambda2
    for(i in 1:length(files_1[[s]])) {
        if(is.null(table_lst_1[[s]][[i]]$x)){
            results_lst_1[[s]][i,2:5] <- NA
        } else{
        results_lst_1[[s]][i,2:5] <- table_lst_1[[s]][[i]]$x
            }
    }
}
results_lst_2 <- list()
for(s in 1:6) {
    results_lst_2[[s]] <- matrix(, nrow = length(files_2[[s]]), ncol = 5)
    colnames(results_lst_2[[s]]) <- c('lambda2', rownames(table_lst_2[[s]][[100]]))
    results_lst_2[[s]][,1] <- lambda2
    for(i in 1:length(files_2[[s]])) {
        if(is.null(table_lst_2[[s]][[i]]$x)){
            results_lst_2[[s]][i,2:5] <- NA
        } else{
        results_lst_2[[s]][i,2:5] <- table_lst_2[[s]][[i]]$x
            }
    }
}

## Plot
l1_pll <- list()
for(s in 1:6) {
    results <- melt(as.data.frame(results_lst[[s]]), id.vars = 'lambda1', variable.name = 'Metrics')
    l1_pll[[s]] <- ggplot(results, aes(lambda1,value)) + 
    geom_point(aes(colour = Metrics), shape = 4) +
    labs(y = 'Annotation Metric', x = bquote(lambda[1]), title = paste('Senario', s)) + 
    geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
    theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
}
pll_1 <- plot_grid(plotlist = l1_pll, ncol=3, align = 'h')
ggsave('/home/bz234/project/Results/CyTOF/simulation/benchmark/lamba1_plot.png', pll_1, width = 15, height = 8)
l2_pll <- list()
for(s in 1:6) {
    results_2 <- melt(as.data.frame(results_lst_2[[s]]), id.vars = 'lambda2', variable.name = 'Metrics')
    l2_pll[[s]] <- ggplot(results_2, aes(lambda2,value)) + 
    geom_point(aes(colour = Metrics), shape = 4) +
    labs(y = 'Annotation Metric', x = bquote(lambda[2]), title = paste('Senario', s)) + 
    geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
    theme(axis.text=element_text(size=12),
    axis.title=element_text(size=14,face="bold"),
    plot.title = element_text(size=14))
}
pll_2 <- plot_grid(plotlist = l2_pll, ncol=3, align = 'h')
ggsave('/home/bz234/project/Results/CyTOF/simulation/benchmark/lamba2_plot.png', pll_2, width = 15, height = 8)