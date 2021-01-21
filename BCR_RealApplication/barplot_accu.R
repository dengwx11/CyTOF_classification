library(reshape2)
library(ggplot2)
library(viridis)

acc_BCR5 <- readRDS("/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/rst_summary/acc_BCR.rds")
acc_ref5 <- readRDS("/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/rst_summary/acc_ref.rds")
acc_BCR6 <- readRDS("/home/bz234/project/Results/CyTOF/BCR_Cytof/bar_accu/acc_bcr_6.rds")
acc_ref6 <- readRDS("/home/bz234/project/Results/CyTOF/BCR_Cytof/bar_accu/acc_ref_6.rds")
acc_BCR7 <- readRDS("/home/bz234/project/Results/CyTOF/BCR_Cytof/bar_accu/acc_bcr_7.rds")
acc_ref7 <- readRDS("/home/bz234/project/Results/CyTOF/BCR_Cytof/bar_accu/acc_ref_7.rds")

acc_BCR_long5 <- melt(acc_BCR5)
acc_BCR_long5$patient.ID <- c(c(1:8),c(1:8))
# options(repr.plot.width = 6, repr.plot.height = 3.5, repr.plot.res = 300)
ggplot(acc_BCR_long5, aes(fill=variable, y=value, x=patient.ID)) + ggtitle("BCR-XL Accuracy 5 Cell types")+
    geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete = T) + 
    scale_x_continuous(breaks=1:8, labels=1:8) +
    ylab('Accuracy') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) + 
    coord_cartesian(ylim=c(0, 1))

acc_ref_long5 <- melt(acc_ref5)
acc_ref_long5$patient.ID <- c(c(1:8),c(1:8))
# options(repr.plot.width = 6, repr.plot.height = 3.5, repr.plot.res = 300)
ggplot(acc_ref_long5, aes(fill=variable, y=value, x=patient.ID)) + ggtitle("Reference Accuracy 5 Cell types")+
    geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete = T) + 
    scale_x_continuous(breaks=1:8, labels=1:8) +
    ylab('Accuracy') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) + 
    coord_cartesian(ylim=c(0, 1))

acc_BCR_long6 <- melt(acc_BCR6)
acc_BCR_long6$patient.ID <- c(c(1:8),c(1:8))
# options(repr.plot.width = 6, repr.plot.height = 3.5, repr.plot.res = 300)
ggplot(acc_BCR_long6, aes(fill=variable, y=value, x=patient.ID)) + ggtitle("BCR-XL Accuracy 6 Cell types")+
    geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete = T) + 
    scale_x_continuous(breaks=1:8, labels=1:8) +
    ylab('Accuracy') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) + 
    coord_cartesian(ylim=c(0, 1))

acc_ref_long6 <- melt(acc_ref6_new)
acc_ref_long6$patient.ID <- c(c(1:8),c(1:8))
# options(repr.plot.width = 6, repr.plot.height = 3.5, repr.plot.res = 300)
ggplot(acc_ref_long6, aes(fill=variable, y=value, x=patient.ID)) + ggtitle("Reference Accuracy 6 Cell types")+
    geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete = T) + 
    scale_x_continuous(breaks=1:8, labels=1:8) +
    ylab('Accuracy') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) +
    coord_cartesian(ylim=c(0, 1))

acc_BCR_long7 <- melt(acc_BCR7)
acc_BCR_long7$patient.ID <- c(c(1:8),c(1:8))
# options(repr.plot.width = 6, repr.plot.height = 3.5, repr.plot.res = 300)
ggplot(acc_BCR_long7, aes(fill=variable, y=value, x=patient.ID)) + ggtitle("BCR-XL Accuracy 7 Cell types")+
    geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete = T) + 
    scale_x_continuous(breaks=1:8, labels=1:8) +
    ylab('Accuracy') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) +
    coord_cartesian(ylim=c(0, 1))

acc_ref_long7 <- melt(acc_ref7_new)
acc_ref_long7$patient.ID <- c(c(1:8),c(1:8))
# options(repr.plot.width = 6, repr.plot.height = 3.5, repr.plot.res = 300)
ggplot(acc_ref_long7, aes(fill=variable, y=value, x=patient.ID)) + ggtitle("Reference Accuracy 7 Cell types")+
    geom_bar(position="dodge", stat="identity") + scale_fill_viridis(discrete = T) + 
    scale_x_continuous(breaks=1:8, labels=1:8) +
    ylab('Accuracy') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) +
    coord_cartesian(ylim=c(0, 1))