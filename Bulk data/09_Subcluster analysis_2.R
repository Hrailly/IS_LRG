rm(list = ls())
#引用包
library(limma)
library(reshape2)
library(ggpubr)
expFile="geneexpr.txt"          #表达输入文件
geneCluFile="Cluster.txt"     #基因分型文件

#读取表达输入文件
rt=read.table(expFile, header=T, sep="\t", check.names=F)
# rownames(rt)=rt[,1]
# rt<-rt[,-1]
rt=as.matrix(rt)
exp=rt
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#读取基因分型文件
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#合并数据
sameSample=intersect(row.names(data), row.names(geneClu))
expClu=cbind(data[sameSample,,drop=F], geneClu[sameSample,,drop=F])

#把数据转换成ggplot2输入文件
data=melt(expClu, id.vars=c("cluster"))
colnames(data)=c("Cluster", "Gene", "Expression")

#设置颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"Cluster"])))]

#绘制箱线图
p=ggboxplot(data, x="Gene", y="Expression", color = "Cluster", 
            ylab="Gene expression",
            xlab="",
            legend.title="Cluster",
            palette = bioCol,
            width=1)
p=p+rotate_x_text(60)
p1=p+stat_compare_means(aes(group=Cluster),
                        symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
                        label = "p.signif")

#输出箱线图
pdf(file="boxplot.pdf", width=6, height=5)
print(p1)
dev.off()

