
rm(list = ls())
library(pheatmap)         #引用包
expFile="geneexpr.txt"         #表达输入文件
clusterFile="Cluster.txt"     #聚类结果文件
cliFile="clinicaldata.txt"           #临床数据文件


#读取输入文件
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
exp=t(exp)
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#合并表达和分型数据
sameSample=intersect(row.names(exp), row.names(cluster))
exp=exp[sameSample, , drop=F]
cluster=cluster[sameSample, , drop=F]
expCluster=cbind(exp, cluster)

Project=gsub("(.*?)\\_.*", "\\1", rownames(expCluster))
library(tidyverse)
expCluster=cbind(expCluster, Project)

#合并临床数据
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(expCluster), row.names(cli))
expCluster=expCluster[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expCluster, cli)

#提取热图数据
data=data[order(data$cluster),]
Type=data[,((ncol(exp)+1):ncol(data))]
Type=Type[,-2]
data=t(data[,1:ncol(exp)])

#聚类颜色
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CluCol=bioCol[1:length(levels(factor(Type$cluster)))]
names(CluCol)=levels(factor(Type$cluster))
ann_colors[["cluster"]]=CluCol

#热图可视化
pdf("heatmap.pdf", height=2.2, width=6)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         cluster_cols =F,
         cluster_rows =F,
         scale="row",
         show_colnames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()










