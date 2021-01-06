rm(list=ls())
require(ncdfFlow)
require(cydar)
require(flowCore)

setwd("/Users/JiehuanSun/LargeData/original data/debarcode/")
prefix = "~/Box Sync/study/Yale/research/IPF projects/cytof/results/Gene_IPFvsCTRL/"

files = list.files(pattern = "*\\.fcs", full.names = TRUE)
file1 = files[grep("CTRL",files)]
file2 = files[grep("IPFH",files)]
file3 = files[grep("IPFL",files)]

files = c(file1,file2,file3)
x = read.ncdfFlowSet(files)
x = x[,-c(1:2,4,6,7,8,9,10,58:61)]
pool.ff = poolCells(x)

myTransform <- function(transformId){
    t = new("transform", .Data = function(x) log10(x+1))
    t@transformationId = transformId
    t
}
myT = myTransform(transformId="LogT")

translist = transformList(colnames(pool.ff), myT) 
proc.ff = transform(pool.ff, translist )

gate1 = outlierGate(proc.ff, "Ir191Di", type="lower")
gate2 = outlierGate(proc.ff, "Ir193Di", type="lower")
gate3 = outlierGate(proc.ff, "Pt195Di", type="upper")

processed.exprs = transform(x, translist)
processed.exprs = Subset(processed.exprs, gate1)
processed.exprs = Subset(processed.exprs, gate2)
processed.exprs = Subset(processed.exprs, gate3)

genes = read.delim(file = "/Users/JiehuanSun/LargeData/original data/genes.txt", sep="")
cellmarkers = read.delim(file = "/Users/JiehuanSun/LargeData/original data/celltype.txt", sep="")
markers = c(as.character(genes[,1]), as.character(cellmarkers[,1]))
names = c(as.character(genes[,2]), as.character(cellmarkers[,2]))
names = strsplit(names,split="_")
names = sapply(names, function(x){x[2]})
processed.exprs = processed.exprs[,match(markers, colnames(processed.exprs))]
colnames(processed.exprs)  = names


samples = sampleNames(processed.exprs)
samples = gsub("BaseFileName_", "", samples)
samples = gsub("\\.fcs", "", samples)

data = NULL
rowname = NULL
for(i in seq_along(samples)){
    tmp = exprs(poolCells(processed.exprs[i]))
    data = rbind(data, tmp)
    rowname = c(rowname, rep(samples[i],nrow(tmp)) )
}

rownames(data) = rowname
colnames(data) = colnames(processed.exprs)

data.tmp = data


set.seed(0824)
data = data[,16:30]
indlist = list()
for( i in 1:ncol(data)){
    print(i)
    xx = data[,i]
    fit = kmeans(xx,centers=3)
    ind = match(fit$cluster,sort(fit$centers,index.return=T)$ix)
    indlist[[i]] = ind
}

indmat = do.call(cbind,indlist)
rownames(indmat) = rownames(data)
colnames(indmat) = colnames(data)

indmat[,1:13] = (indmat[,1:13]>1) + 1

options(stringsAsFactors = FALSE)
types = read.delim(file="/Users/JiehuanSun/LargeData/original data/q-cellSubsetsV3.txt")

rownames(types) = types[,1]
types = types[,-1]
typecode = matrix(0,nrow=nrow(types),ncol=ncol(types))
rownames(typecode) = rownames(types)
typecode[types==""] = NA
typecode[types=="-"] = 1
typecode[types=="+"] = 2
typecode[types=="++"] = 3

res = matrix(0,nrow(data),nrow(typecode))
for(i in 1:nrow(typecode)){
    print(i)
    code = typecode[i,]
    tmp = indmat[,!is.na(code)]
    code  = code[!is.na(code)]
    tt = t(tmp) - code
    res[,i] = apply(tt,2,function(x){all(x==0)})
}
rownames(res) = rownames(data)
colnames(res) = rownames(typecode)

unclass = apply(res,1,sum)
res = res[unclass==1, ]
data.tmp = data.tmp[unclass==1,]

df = split(as.data.frame(res), rownames(res))
count.tmp = sapply(df, function(x){
    apply(x,2,sum)
})
count.res = as.data.frame(t(count.tmp))
count.res = data.frame(patient = rownames(count.res), count.res)
write.table(count.res, file="count_cell.txt",row.names = FALSE, sep='\t',quote=FALSE)

res.tmp = res
require(lme4)
res.all = NULL

for(ii in 1:ncol(res.tmp)){
    data = data.tmp[res.tmp[,ii]==1,]
    group = rep("IPF",nrow(data))
    group[grep("CTRL",rownames(data))] = "CTRL"
    group = factor(group,levels = c("CTRL", "IPF"))
    t.vec = NULL
    for(i in 1:15){
        print(i)
        tmp = data.frame(ID=rownames(data),gene=data[,i],group=group)
        fit = lmer(gene~group+(1|ID),data=tmp)
        t.vec = c(t.vec, summary(fit)$coef[2,3])
    }
    
    p = (1-pnorm(abs(t.vec)))*2
    
    res = data.frame(Subset = colnames(res.tmp)[ii], Gene = colnames(data)[1:15], 
                     stat = round(t.vec,4), p = round(p,4), 
                     fdr = round(p.adjust(p,'fdr'),4) )
    res = res[order(res[,4]),]
    
    res.all = rbind(res.all, res)
}

write.table(res.all, file=paste(prefix,"Cell_Gene_IPFvsCTRL.txt"),
            quote=F,sep="\t",row.names=F)


#### plot each individual figure ####

require(limma)
require(ggplot2)
require(ggpubr)
marker.name = colnames(data.tmp)[1:15]
for(ii in 1:ncol(res.tmp)){
    data = data.tmp[res.tmp[,ii]==1,1:15]
    data = avereps(data)
    nnn = colnames(res.tmp)[ii]
    group = rep("IPF",nrow(data))
    group[grep("CTRL",rownames(data))] = "CTRL"
    group = factor(group,levels = c("CTRL", "IPF"))
    rownames(data) = NULL
    data = data.frame(Group= group, data)
    p = ggboxplot(data = data, x = "Group", y = marker.name, 
                  color = "Group", palette = "npg", add = "jitter",
                  merge='flip', ylab=paste('Expression in ', nnn,sep=''))
    
    df = NULL
    ymax3.vec = NULL
    label.vec = NULL
    for(i in 2:16){
        ymax1 = min(max(data[,i])*1.15)
        ymax2 = ymax1 + max(data[,2:16])*0.02
        ymax3 = ymax2 + max(data[,2:16])*0.05
        ppp = res.all[res.all$Subset==nnn & res.all$Gene==marker.name[i-1],4]
        if(ppp<=0.0001){
            label = "****"
        }else if(ppp<=0.001){
            label = "***"
        }else if(ppp<=0.01){
            label = "**"
        }else if(ppp<=0.05){
            label = "*"
        }else{label = 'ns'}
        #label = paste("p = ", ppp, sep='')
        label.vec = c(label.vec,label)
        ymax3.vec = c(ymax3.vec,ymax3)
        df.tmp = data.frame(a = c(0.75+i-2, 0.75+i-2, 1.25+i-2 ,1.25+i-2), 
                         b = c(ymax1, ymax2, ymax2, ymax1),c=i)
        df = rbind(df,df.tmp)
    }
    
    p + geom_line(data = df, aes(x = a, y = b,group=c)) + 
        annotate("text", x = 1:15, y = ymax3.vec, label = label.vec, size = 4)+
        theme(axis.text.x = element_text(angle = 60, hjust = 1))
    
    ggsave(file=paste(prefix,nnn,"_ Cell_Gene_IPFvsCTRL.pdf",sep=''),
           height = 6, width=15)
}




#### plot each individual figure with significant ones ####

require(limma)
require(ggplot2)
require(ggpubr)

col_npg = pal_npg("nrc")(10)

marker.name = colnames(data.tmp)[1:15]
for(ii in 1:ncol(res.tmp)){
    
    nnn = colnames(res.tmp)[ii]
    ggg = res.all[res.all$Subset==nnn & res.all$p<0.05,2]
    
    if(length(ggg)>0){
        data = data.tmp[res.tmp[,ii]==1, ggg, drop=FALSE]
        data = as.data.frame(avereps(data))
        # data = data[grep("IPF",rownames(data)),,drop=FALSE]
        
        group = rep("IPF",nrow(data))
        group[grep("CTRL",rownames(data))] = "CTRL"
        group = factor(group,levels = c("CTRL", "IPF"))
        rownames(data) = NULL
        data = data.frame(Group = group, data)
        
        p = ggboxplot(data = data, x = "Group", y = ggg, 
                      color = "Group", palette = col_npg[c(2,1)], add = "jitter",
                      merge='flip', ylab=paste('Expression in ', nnn,sep=''),
                      lwd=1.2)
        
        df = NULL
        ymax3.vec = NULL
        label.vec = NULL
        for(i in 2:ncol(data)){
            ymax1 = max(data[,i]) + max(data[,2:ncol(data)])*0.05
            ymax2 = ymax1 + max(data[,2:ncol(data)])*0.02
            ymax3 = ymax2 + max(data[,2:ncol(data)])*0.05
            ppp = res.all[res.all$Subset==nnn & res.all$Gene==colnames(data)[i],4]
            if(ppp<=0.0001){
                label = "****"
            }else if(ppp<=0.001){
                label = "***"
            }else if(ppp<=0.01){
                label = "**"
            }else if(ppp<=0.05){
                label = "*"
            }else{label = 'ns'}
            #label = paste("p = ", ppp, sep='')
            label.vec = c(label.vec,label)
            ymax3.vec = c(ymax3.vec,ymax3)
            if(length(ggg)==1){
                df.tmp = data.frame(a = c(0.75+i-2, 0.75+i-2, 1.25+i-1 ,1.25+i-1), 
                                    b = c(ymax1, ymax2, ymax2, ymax1),c=i)
            }else{
                df.tmp = data.frame(a = c(0.75+i-2, 0.75+i-2, 1.25+i-2 ,1.25+i-2), 
                                    b = c(ymax1, ymax2, ymax2, ymax1),c=i)
            }
            
            df = rbind(df,df.tmp)
        }
        if(length(ggg)==1){
            xpoint = 1.5
        }else{xpoint = 1:(ncol(data)-1)}
        p + geom_line(data = df, aes(x = a, y = b,group=c)) + 
            annotate("text", x = xpoint, y = ymax3.vec, label = label.vec, size = 4)+
            theme(axis.text.x = element_text(angle = 60, hjust = 1))
        
        ggsave(file=paste(prefix,nnn,"_ Cell_Gene_IPFvsCTRL.pdf",sep=''),
               height = 4, width = 2 + ncol(data))
    }
}
