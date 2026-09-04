##ssGESA
rm(list = ls())
#引用包
library(reshape2)
library(tidyverse)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)
expFile="merge矩阵new.csv"           #表达输入文件
gmtFile="immune.gmt"          #免疫数据集文件
clusterFile="clinicaldata.txt"       #分型输入文件

#读取表达输入文件,并对输入文件整理
rt=read.csv(expFile, header=T, sep=",", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#读取基因集文件
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA分析
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#对ssGSEA打分进行矫正
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#输出ssGSEA打分结果
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="merge_ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#读取分型文件
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#数据合并
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)
colnames(scoreCluster)<-str_replace_all(colnames(scoreCluster),"na","")
#把数据转换成ggplot2输入文件
data=melt(scoreCluster, id.vars=c("condition"))
colnames(data)=c("cluster", "Immune", "Fraction")

#绘制箱线图
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"cluster"])))]
library(tidyverse)
data$Immune=str_replace_all(data$Immune,"na","")
p=ggboxplot(data, x="Immune", y="Fraction", color="cluster", 
            ylab="Immune infiltration",
            xlab="",
            legend.title="cluster",
            palette=bioCol)
p=p+rotate_x_text(50)
pdf(file="merge_ssGSEA_boxplot.pdf", width=8, height=6.5)                          #输出图片文件
p+stat_compare_means(aes(group=cluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),label = "p.signif")
dev.off()


colnames(scoreCluster)<-str_replace_all(colnames(scoreCluster),"na","")
data <- scoreCluster
data=dplyr::select(data,-"condition")
library(corrplot)

#相关性矩阵
M=cor(data)
res1=cor.mtest(data, conf.level = 0.95)

#绘制相关性图形
pdf(file="merge_cor.pdf", width=8, height=8)
corrplot(M,
         order="original",
         method = "circle",
         type = "upper",
         tl.cex=0.8, pch=T,
         p.mat = res1$p,
         insig = "label_sig",
         pch.cex = 1.6,
         sig.level=0.05,
         number.cex = 1,
         col=colorRampPalette(c("blue", "white", "red"))(50),
         tl.col="black")
dev.off()
