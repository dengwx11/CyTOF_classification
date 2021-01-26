library(Seurat)
library(reshape2)
library(ggplot2)
library(cowplot)

## 5 cell types
seur_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern='*seur*', full=T, ignore.case = TRUE)
rst_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern='*rst*', full=T, ignore.case = TRUE)
BCR_XL_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern=glob2rx('*BCR-XL.rds'), full=T, ignore.case = TRUE)
Reference_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/5celltype', pattern=glob2rx('*Reference.rds'), full=T, ignore.case = TRUE)
rst_BCR_files <- intersect(rst_files,BCR_XL_files)
rst_Reference_files <- intersect(rst_files,Reference_files)
seur_BCR_files <- intersect(seur_files,BCR_XL_files)
seur_Reference_files <- intersect(seur_files,Reference_files)

# BCR
ct_cnt_pred <- data.frame(matrix(0,5,nrow=1,ncol=5))
ct_cnt_label <- data.frame(matrix(0,5,nrow=1,ncol=5))
colnames(ct_cnt_pred) <- colnames(ct_cnt_label) <- c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'NK cells')

for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_BCR_files[sample.idx])
    seur <- readRDS(seur_BCR_files[sample.idx])
    truth <- seur$label
    
    ## subsetting
    H_est_subsetting <- rst$H[,-1]
    H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))
    
    idx.subsetting = c()
    for(i in 1:length(seur$label)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
        }

    }  
    # cell type proportions
    seur_sub <- subset(seur, cells=Cells(seur)[idx.subsetting])
    ct_cnt_pred <- rbind(ct_cnt_pred,table(seur_sub$pred)[c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'NK cells')])
    ct_cnt_label <- rbind(ct_cnt_label,table(seur_sub$label)[c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'NK cells')])                                                       
}

# Reference
ct_cnt_ref_pred <- data.frame(matrix(0,5,nrow=1,ncol=5))
ct_cnt_ref_label <- data.frame(matrix(0,5,nrow=1,ncol=5))
colnames(ct_cnt_ref_pred) <- colnames(ct_cnt_ref_label) <- c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'NK cells')

for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_Reference_files[sample.idx])
    seur <- readRDS(seur_Reference_files[sample.idx])
    truth <- seur$label
    
    ## subsetting
    H_est_subsetting <- rst$H[,-1]
    H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))
    
    idx.subsetting = c()
    for(i in 1:length(seur$label)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
        }

    }  
    # cell type proportions
    if(is.null(idx.subsetting)) {
        ct_cnt_ref_pred <- rbind(ct_cnt_ref_pred, NA)
    } else{
        seur_sub <- subset(seur, cells=Cells(seur)[idx.subsetting])
    ct_cnt_ref_pred <- rbind(ct_cnt_ref_pred,table(seur_sub$pred)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                                    'naïve B cells', 'NK cells')])
    }
    ct_cnt_ref_label <- rbind(ct_cnt_ref_label,table(seur_sub$label)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                                       'naïve B cells', 'NK cells')])                          
}

ct_cnt_pred <- ct_cnt_pred[-1,]
ct_cnt_label <- ct_cnt_label[-1,]
ct_cnt_ref_pred <- ct_cnt_ref_pred[-1,]
ct_cnt_ref_label <- ct_cnt_ref_label[-1,]

ct_prop_pred <- apply(ct_cnt_pred,1,function(x) x/sum(x))
ct_prop_label <- apply(ct_cnt_label,1,function(x) x/sum(x))
ct_prop_ref_pred <- apply(ct_cnt_ref_pred,1,function(x) x/sum(x))
ct_prop_ref_label <- apply(ct_cnt_ref_label,1,function(x) x/sum(x))
                           
ct_prop_pred <- ct_prop_pred[,-8]
ct_prop_label <- ct_prop_label[,-8]  
                           
ct_prop_pred <- data.frame(ct_prop_pred)
ct_prop_label <- data.frame(ct_prop_label)  
ct_prop_ref_pred <- data.frame(ct_prop_ref_pred)
ct_prop_ref_label <- data.frame(ct_prop_ref_label)  


ct_prop_pred$celltype <- rownames(ct_prop_pred)
ct_prop_label$celltype <- rownames(ct_prop_label)
ct_prop_ref_pred$celltype <- rownames(ct_prop_ref_pred)
ct_prop_ref_label$celltype <- rownames(ct_prop_ref_label)
                           
ct_prop_pred_long <-  melt(ct_prop_pred, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_label_long <-  melt(ct_prop_label, id.vars = 'celltype', variable.name = 'patient.ID')   
ct_prop_pred_ref_long <-  melt(ct_prop_ref_pred, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_label_ref_long <-  melt(ct_prop_ref_label, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_pred_long$group = "prediction-BCR"
ct_prop_label_long$group = "true-BCR"
ct_prop_pred_ref_long$group = "prediction-Ref"
ct_prop_label_ref_long$group = "true-Ref"
ct_prop_long <- rbind(ct_prop_pred_long,ct_prop_label_long)    
ct_prop_ref_long <- rbind(ct_prop_pred_ref_long,ct_prop_label_ref_long)    

ct_prop_long_all <- rbind(ct_prop_long,ct_prop_ref_long)
options(repr.plot.width = 10, repr.plot.height = 7, repr.plot.res = 300)
ggplot(ct_prop_long_all, aes(x=celltype, y=value,fill=group)) + 
    geom_boxplot()+
    xlab('Cell type') +
    ylab('Proportion') +
    ggtitle('Cell Type Proportion (5 Types)') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) 
                           
## 6 cell types
seur_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/6celltype', pattern='*seur*', full=T, ignore.case = TRUE)
rst_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/6celltype', pattern='*rst*', full=T, ignore.case = TRUE)
BCR_XL_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/6celltype', pattern=glob2rx('*BCR-XL.rds'), full=T, ignore.case = TRUE)
Reference_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/6celltype', pattern=glob2rx('*Reference.rds'), full=T, ignore.case = TRUE)
rst_BCR_files <- intersect(rst_files,BCR_XL_files)
rst_Reference_files <- intersect(rst_files,Reference_files)
seur_BCR_files <- intersect(seur_files,BCR_XL_files)
seur_Reference_files <- intersect(seur_files,Reference_files)

# BCR
ct_cnt_pred <- data.frame(matrix(0,5,nrow=1,ncol=6))
ct_cnt_label <- data.frame(matrix(0,5,nrow=1,ncol=6))
colnames(ct_cnt_pred) <- colnames(ct_cnt_label) <- c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'memory B cells', 'NK cells')

for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_BCR_files[sample.idx])
    seur <- readRDS(seur_BCR_files[sample.idx])
    truth <- seur$label
    
    ## subsetting
    H_est_subsetting <- rst$H[,-1]
    H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))
    
    idx.subsetting = c()
    for(i in 1:length(seur$label)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
        }

    }  
    # cell type proportions
    seur_sub <- subset(seur, cells=Cells(seur)[idx.subsetting])
    ct_cnt_pred <- rbind(ct_cnt_pred,table(seur_sub$pred)[c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'memory B cells', 'NK cells')])
    ct_cnt_label <- rbind(ct_cnt_label,table(seur_sub$label)[c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'memory B cells', 'NK cells')])                          
                              
}
                              
# Reference
ct_cnt_ref_pred <- data.frame(matrix(0,5,nrow=1,ncol=6))
ct_cnt_ref_label <- data.frame(matrix(0,5,nrow=1,ncol=6))
colnames(ct_cnt_ref_pred) <- colnames(ct_cnt_ref_label) <- c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'memory B cells', 'NK cells')

for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_Reference_files[sample.idx])
    seur <- readRDS(seur_Reference_files[sample.idx])
    truth <- seur$label
    
    ## subsetting
    H_est_subsetting <- rst$H[,-1]
    H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))
    
    idx.subsetting = c()
    for(i in 1:length(seur$label)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
        }

    }  
    # cell type proportions
    if(is.null(idx.subsetting)) {
        ct_cnt_ref_pred <- rbind(ct_cnt_ref_pred, NA)
    } else{
        seur_sub <- subset(seur, cells=Cells(seur)[idx.subsetting])
    ct_cnt_ref_pred <- rbind(ct_cnt_ref_pred,table(seur_sub$pred)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                                    'naïve B cells', 'memory B cells', 'NK cells')])
    }
    ct_cnt_ref_label <- rbind(ct_cnt_ref_label,table(seur_sub$label)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                                       'naïve B cells', 'memory B cells', 'NK cells')])                          
}

ct_cnt_pred <- ct_cnt_pred[-1,]
ct_cnt_label <- ct_cnt_label[-1,]
ct_cnt_ref_pred <- ct_cnt_ref_pred[-1,]
ct_cnt_ref_label <- ct_cnt_ref_label[-1,]

ct_prop_pred <- apply(ct_cnt_pred,1,function(x) x/sum(x))
ct_prop_label <- apply(ct_cnt_label,1,function(x) x/sum(x))
ct_prop_ref_pred <- apply(ct_cnt_ref_pred,1,function(x) x/sum(x))
ct_prop_ref_label <- apply(ct_cnt_ref_label,1,function(x) x/sum(x))

ct_prop_pred <- ct_prop_pred[,-3]
ct_prop_label <- ct_prop_label[,-3]  
ct_prop_ref_pred <- ct_prop_ref_pred[,-c(4, 6)]
ct_prop_ref_label <- ct_prop_ref_label[,-c(4, 6)] 
                           
ct_prop_pred <- data.frame(ct_prop_pred)
ct_prop_label <- data.frame(ct_prop_label)  
ct_prop_ref_pred <- data.frame(ct_prop_ref_pred)
ct_prop_ref_label <- data.frame(ct_prop_ref_label)  


ct_prop_pred$celltype <- rownames(ct_prop_pred)
ct_prop_label$celltype <- rownames(ct_prop_label)
ct_prop_ref_pred$celltype <- rownames(ct_prop_ref_pred)
ct_prop_ref_label$celltype <- rownames(ct_prop_ref_label)

saveRDS(ct_prop_pred,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_pred_BCR_6type.rds")
saveRDS(ct_prop_label,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_label_BCR_6type.rds")
saveRDS(ct_prop_ref_pred,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_pred_ref_6type.rds")
saveRDS(ct_prop_ref_label,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_label_ref_6type.rds")

ct_prop_pred_long <-  melt(ct_prop_pred, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_label_long <-  melt(ct_prop_label, id.vars = 'celltype', variable.name = 'patient.ID')   
ct_prop_pred_ref_long <-  melt(ct_prop_ref_pred, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_label_ref_long <-  melt(ct_prop_ref_label, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_pred_long$group = "prediction-BCR"
ct_prop_label_long$group = "true-BCR"
ct_prop_pred_ref_long$group = "prediction-Ref"
ct_prop_label_ref_long$group = "true-Ref"
ct_prop_long <- rbind(ct_prop_pred_long,ct_prop_label_long)    
ct_prop_ref_long <- rbind(ct_prop_pred_ref_long,ct_prop_label_ref_long)    

ct_prop_long_all <- rbind(ct_prop_long,ct_prop_ref_long)
options(repr.plot.width = 10, repr.plot.height = 7, repr.plot.res = 300)
ggplot(ct_prop_long_all, aes(x=celltype, y=value,fill=group)) + 
    geom_boxplot()+
    xlab('Cell type') +
    ylab('Proportion') +
    ggtitle('Cell Type Proportion (6 Types)') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) 

## 7 cell types
seur_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/7celltype', pattern='*seur*', full=T, ignore.case = TRUE)
rst_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/7celltype', pattern='*rst*', full=T, ignore.case = TRUE)
BCR_XL_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/7celltype', pattern=glob2rx('*BCR-XL.rds'), full=T, ignore.case = TRUE)
Reference_files <- list.files(path = '/gpfs/ysm/pi/zhao-data/wd262/new_cytof/write/BCR/7celltype', pattern=glob2rx('*Reference.rds'), full=T, ignore.case = TRUE)
rst_BCR_files <- intersect(rst_files,BCR_XL_files)
rst_Reference_files <- intersect(rst_files,Reference_files)
seur_BCR_files <- intersect(seur_files,BCR_XL_files)
seur_Reference_files <- intersect(seur_files,Reference_files)

# BCR
ct_cnt_pred7 <- data.frame(matrix(0,5,nrow=1,ncol=7))
ct_cnt_label7 <- data.frame(matrix(0,5,nrow=1,ncol=7))
colnames(ct_cnt_pred7) <- colnames(ct_cnt_label7) <- c('CD4 T-cells','CD8 T-cells','monocytes', 'naïve B cells', 
                                                     'memory B cells', 'NK cells', 'DC')

for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_BCR_files[sample.idx])
    seur <- readRDS(seur_BCR_files[sample.idx])
    truth <- seur$label
    
    ## subsetting
    H_est_subsetting <- rst$H[,-1]
    H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))
    
    idx.subsetting = c()
    for(i in 1:length(seur$label)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
        }

    }  
    # cell type proportions
    seur_sub <- subset(seur, cells=Cells(seur)[idx.subsetting])
    ct_cnt_pred7 <- rbind(ct_cnt_pred7,table(seur_sub$pred)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                              'naïve B cells', 'memory B cells', 'NK cells', 'DC')])
    ct_cnt_label7 <- rbind(ct_cnt_label7,table(seur_sub$label)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                                 'naïve B cells', 'memory B cells', 'NK cells', 'DC')])                          
                              
}
                              
# Reference
ct_cnt_ref_pred7 <- data.frame(matrix(0,5,nrow=1,ncol=7))
ct_cnt_ref_label7 <- data.frame(matrix(0,5,nrow=1,ncol=7))
colnames(ct_cnt_ref_pred7) <- colnames(ct_cnt_ref_label7) <- c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                               'naïve B cells', 'memory B cells', 'NK cells', 'DC')

for(sample.idx in 1:length(rst_BCR_files)){
    rst <- readRDS(rst_Reference_files[sample.idx])
    seur <- readRDS(seur_Reference_files[sample.idx])
    truth <- seur$label
    
    ## subsetting
    H_est_subsetting <- rst$H[,-1]
    H_est_subsetting <- apply(H_est_subsetting,2,function(x) x/sum(x))
    
    idx.subsetting = c()
    for(i in 1:length(seur$label)){
        if(max(H_est_subsetting[,i])>0.5){
            idx.subsetting = append(idx.subsetting,i)
        }

    }  
    # cell type proportions
    if(is.null(idx.subsetting)) {
        ct_cnt_ref_pred7 <- rbind(ct_cnt_ref_pred, NA)
    } else{
        seur_sub <- subset(seur, cells=Cells(seur)[idx.subsetting])
    ct_cnt_ref_pred7 <- rbind(ct_cnt_ref_pred7,table(seur_sub$pred)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                                    'naïve B cells', 'memory B cells', 'NK cells', 
                                                                      'DC')])
    }
    ct_cnt_ref_label7 <- rbind(ct_cnt_ref_label7,table(seur_sub$label)[c('CD4 T-cells','CD8 T-cells','monocytes', 
                                                                       'naïve B cells', 'memory B cells', 'NK cells', 
                                                                         'DC')])                          
}

ct_cnt_pred7 <- ct_cnt_pred7[-1,]
ct_cnt_label7 <- ct_cnt_label7[-1,]
ct_cnt_ref_pred7 <- ct_cnt_ref_pred7[-1,]
ct_cnt_ref_label7 <- ct_cnt_ref_label7[-1,]

ct_prop_pred7 <- apply(ct_cnt_pred7,1,function(x) x/sum(x))
ct_prop_label7 <- apply(ct_cnt_label7,1,function(x) x/sum(x))
ct_prop_ref_pred7 <- apply(ct_cnt_ref_pred7,1,function(x) x/sum(x))
ct_prop_ref_label7 <- apply(ct_cnt_ref_label7,1,function(x) x/sum(x))

ct_prop_pred7 <- data.frame(ct_prop_pred7)
ct_prop_label7 <- data.frame(ct_prop_label7)  
ct_prop_ref_pred7 <- data.frame(ct_prop_ref_pred7)
ct_prop_ref_label7 <- data.frame(ct_prop_ref_label7)  

ct_prop_pred7$celltype <- rownames(ct_prop_pred7)
ct_prop_label7$celltype <- rownames(ct_prop_label7)
ct_prop_ref_pred7$celltype <- rownames(ct_prop_ref_pred7)
ct_prop_ref_label7$celltype <- rownames(ct_prop_ref_label7)

saveRDS(ct_prop_pred7,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_pred_BCR_7type.rds")
saveRDS(ct_prop_label7,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_label_BCR_7type.rds")
saveRDS(ct_prop_ref_pred7,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_pred_ref_7type.rds")
saveRDS(ct_prop_ref_label7,"/home/bz234/project/Results/CyTOF/BCR_Cytof/type_prop/sub_prop_label_ref_7type.rds")

ct_prop_pred_long7 <-  melt(ct_prop_pred7, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_label_long7 <-  melt(ct_prop_label7, id.vars = 'celltype', variable.name = 'patient.ID')   
ct_prop_pred_ref_long7 <-  melt(ct_prop_ref_pred7, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_label_ref_long7 <-  melt(ct_prop_ref_label7, id.vars = 'celltype', variable.name = 'patient.ID')
ct_prop_pred_long7$group = "prediction-BCR"
ct_prop_label_long7$group = "true-BCR"
ct_prop_pred_ref_long7$group = "prediction-Ref"
ct_prop_label_ref_long7$group = "true-Ref"
ct_prop_long7 <- rbind(ct_prop_pred_long7,ct_prop_label_long7)    
ct_prop_ref_long7 <- rbind(ct_prop_pred_ref_long7,ct_prop_label_ref_long7)  

ct_prop_long_all7 <- rbind(ct_prop_long7,ct_prop_ref_long7)
options(repr.plot.width = 10, repr.plot.height = 7, repr.plot.res = 300)
ggplot(ct_prop_long_all7, aes(x=celltype, y=value,fill=group)) + 
    geom_boxplot()+
    xlab('Cell type') +
    ylab('Proportion') +
    ggtitle('Cell Type Proportion (7 Types)') +
    theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        plot.title = element_text(size=14)) 

