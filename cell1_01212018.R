rm(list=ls())
require(ncdfFlow)
require(cydar)
require(flowCore)


setwd("data0/cyTOF/cytof/data")
prefix = "data0/cyTOF/cytof/data/"

files = list.files(pattern = "*\\.fcs", full.names = TRUE)
#file1 = files[grep("CTRL",files)]
#file2 = files[grep("IPFH",files)]
#file3 = files[grep("IPFL",files)]

contrl<-c(8,9,12,13,18)
low<-c(4:7,11,14,15,19,20)
high<-c(1:3,10,16,17)

#files = c(file1,file2,file3)
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

genes = read.delim(file = "/Users/wenxuandeng/GoogleDrive/sucksalt/SC/multisample_SC/data0/cyTOF/cytof/cytof files/genes.txt", sep="")
cellmarkers = read.delim(file = "/Users/wenxuandeng/GoogleDrive/sucksalt/SC/multisample_SC/data0/cyTOF/cytof/cytof files/celltype.txt", sep="")
markers = c(as.character(genes[,1]), as.character(cellmarkers[,1]))
names = c(as.character(genes[,2]), as.character(cellmarkers[,2]))
names = strsplit(names,split="_")
names = sapply(names, function(x){x[2]})
processed.exprs = processed.exprs[,match(markers, colnames(processed.exprs))]
colnames(processed.exprs)  = names


samples = sampleNames(processed.exprs)
samples = gsub("BaseFileName_", "", samples)
samples = gsub("\\.fcs", "", samples)
samples[contrl]<-paste(samples[contrl],"CTRL",sep="_")
samples[high]<-paste(samples[high],"IPFH",sep="_")
samples[low]<-paste(samples[low],"IPFL",sep="_")

data = NULL
rowname = NULL
for(i in seq_along(samples)){
    tmp = exprs(poolCells(processed.exprs[i]))
    data = rbind(data, tmp)
    rowname = c(rowname, rep(samples[i],nrow(tmp)) )
}

rownames(data) = rowname
colnames(data) = colnames(processed.exprs)
data.total<-data

data.total<-cbind(data.total, batch=rownames(data.total))

write.table(data.total,"data_total.txt",sep='\t',row.names=F,col.names=T,quote=F)

####################################################################


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
types = read.delim(file="/Users/wenxuandeng/GoogleDrive/sucksalt/SC/multisample_SC/data0/cyTOF/cytof/cytof files/q-cellSubsetsV3.txt")

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
res = res[unclass==1,]
# data.tmp = data.tmp[unclass==1,]

require(limma)
xx = avereps(res[,1],ID=rownames(res))
prop = apply(res,2,function(x){avereps(x,ID=rownames(res))})
rownames(prop) = rownames(xx)
# name = gsub("\\d+","",rownames(prop))
# name = factor(name, levels=c("CTRL","IPFL","IPFH"))
# 
# pdf(file=paste(prefix,"cell_prop.pdf",sep=''),height = 6, width=12)
# layout(matrix(1:18,nrow=3,byrow=TRUE))
# for(i in 1:18){
#     boxplot(prop[,i]~name,main=colnames(prop)[i],ylab="Proportions")
# }
# dev.off()

group = rep("IPF",nrow(res))
group[grep("CTRL",rownames(res))] = "CTRL"
group = factor(group,levels = c("CTRL", "IPF"))

require(lme4)
t.vec = NULL
for(i in 1:ncol(res)){
    print(i)
    tmp = data.frame(ID=rownames(res),gene=res[,i],group=group)
    
    # result = tryCatch({
    #     fit = glmer(gene~group+(1|ID),data=tmp,family = "binomial")
    # }, warning = function(w) {
    #     print(colnames(res)[i])
    # })
    
    result = tryCatch({
        fit = glmer(gene~group+(1|ID),data=tmp,family = "binomial")
    }, error = function(ee) {
        print(ee)
        print(colnames(res)[i])
    })
    
    t.vec = c(t.vec, summary(fit)$coef[2,3])
}

p = (1-pnorm(abs(t.vec)))*2

res = data.frame(CellTypes = colnames(res), stat= round(t.vec,4), p =  signif(p,2), 
                 fdr = signif(p.adjust(p,'fdr'),2) )
res = res[order(res[,3]),]
write.table(res, file="CellTypes_IPFvsCTRL.txt",quote=F,sep="\t",row.names=F)


### figures ####
require(ggsci)
col_npg = pal_npg("nrc")(10)

group = rep("IPF",nrow(prop))
group[grep("CTRL",rownames(prop))] = "CTRL"
group = factor(group,levels = c("CTRL", "IPF"))

require(ggplot2)
require(ggpubr)
for(i in 1:ncol(prop)){
    nnn = colnames(prop)[i]
    # pdf(file=paste(prefix,nnn,"_CellTypes_IPFvsCTRL.pdf",sep=''),
    #     height = 6, width=12)
    yy = prop[,i]
    df = data.frame(Proportions=yy,Group=group)
    # boxplot(prop[,i] ~ group, main=nnn, ylab="Proportions", 
    #         bty='L',outline=FALSE,frame=FALSE,
    #         ylim=c(max(c(min(yy)-0.1,0)) ,max(yy)*1.2 ) )
    # 
    # stripchart(prop[,i] ~ group, vertical = TRUE, 
    #            method = "jitter", add = TRUE, pch = 20, col = c('blue','red') )
    p = ggboxplot(df, x = "Group", y = "Proportions", 
                  color = "Group", palette = col_npg[c(2,1)], add = "jitter",lwd=1.2)
    ymax1 = min(c(max(yy)*1.15,1))
    ymax2 = min(c(ymax1 + max(yy)*0.02,1))
    ymax3 = min(c(ymax2 + max(yy)*0.05,1))
    label = paste("p = ", res$p[res$CellTypes==nnn], sep='')
    df1 = data.frame(a = c(1, 1, 2 ,2), b = c(ymax1, ymax2, ymax2, ymax1))
    p + geom_line(data = df1, aes(x = a, y = b)) + 
        annotate("text", x = 1.5, y = ymax3, label = label, size = 4) 
    ggsave(file=paste(prefix,nnn,"_CellTypes_IPFvsCTRL.pdf",sep=''),
           height = 6, width=6)
}



#### plot each individual figure with significant ones ####
prop.tmp = prop
require(limma)
require(ggplot2)
require(ggpubr)

ggg = res[res$p<0.05,1]
prop = prop[,ggg]
group = rep("IPF",nrow(prop))
group[grep("CTRL",rownames(prop))] = "CTRL"
group = factor(group,levels = c("CTRL", "IPF"))
rownames(prop) = NULL
prop = data.frame(Group = group, prop)
colnames(prop)[2:ncol(prop)] = ggg
p = ggboxplot(data = prop, x = "Group", y = ggg, 
              color = "Group", palette = col_npg[c(2,1)], add = "jitter",
              merge='flip', ylab='Proportions',lwd=1.2)

df = NULL
ymax3.vec = NULL
label.vec = NULL
for(i in 2:ncol(prop)){
    ymax1 = max(prop[,i]) + max(prop[,2:ncol(prop)])*0.05
    ymax2 = ymax1 + max(prop[,2:ncol(prop)])*0.02
    ymax3 = ymax2 + max(prop[,2:ncol(prop)])*0.05
    ppp = res[res[,1]==colnames(prop)[i],3]
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
}else{xpoint = 1:(ncol(prop)-1)}
p + geom_line(data = df, aes(x = a, y = b,group=c)) + 
    annotate("text", x = xpoint, y = ymax3.vec, label = label.vec, size = 4)+
    theme(axis.text.x = element_text(angle = 60, hjust = 1))

ggsave(file=paste(prefix,"Signif_CellTypes_IPFvsCTRL.pdf",sep=''),
       height = 6, width = 4 + ncol(prop))




