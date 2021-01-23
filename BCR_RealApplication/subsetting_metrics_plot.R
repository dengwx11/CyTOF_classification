df_files <- list.files(path = '/home/bz234/project/Results/CyTOF/BCR_Cytof/subset_metrics/', full=T)
titles <- list.files(path = '/home/bz234/project/Results/CyTOF/BCR_Cytof/subset_metrics/')
titles <- substr(titles, 1, nchar(titles) - 4)

df_lst <- list()
for(i in 1:length(df_files)) {
    df_lst[[i]] <- read.table(df_files[i])
    df_lst[[i]] <- df_lst[[i]][-100,]
}
df_lst_long <- list()
for(i in 1:length(df_files)) {
    df_lst_long[[i]] <- melt(df_lst[[i]][,-6], id.vars = 'cutoff', variable.name = 'Metrics')
    
}

## Plot five metrics
plt_lst <- list()
for(i in 1:length(df_files)) {
    plt_lst[[i]] <- ggplot(df_lst_long[[i]], aes(cutoff,value)) + 
                    geom_point(aes(colour = Metrics), shape = 4) +
                    labs(y = 'Annotation Metric', x = 'Cutoff', title = titles[[i]]) + 
                    geom_line(aes(color = Metrics, linetype = Metrics), size=2) +
                    # ylim(0.8, 1) +
                    theme(axis.text=element_text(size=12),
                    axis.title=element_text(size=14,face="bold"),
                    plot.title = element_text(size=14))    
}
options(repr.plot.width = 30, repr.plot.height = 25, repr.plot.res = 300)
plot_grid(plotlist = plt_lst, ncol=4)

## Plot N
plt_lst_n <- list()
for(i in 1:length(df_files)) {
    plt_lst_n[[i]] <- ggplot(df_lst[[i]], aes(cutoff,N)) + 
                        geom_point(shape = 4) +
                        labs(y = 'Total cell number', x = 'Cutoff', title = titles[[i]]) + 
                        geom_line(size=2) +
                        # ylim(0.8, 1) +
                        theme(axis.text=element_text(size=12),
                        axis.title=element_text(size=14,face="bold"),
                        plot.title = element_text(size=14))
}
plot_grid(plotlist = plt_lst_n, ncol=4)