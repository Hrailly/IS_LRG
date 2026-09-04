rm(list = ls())
load("merge_IS.RDATA")
Merge<-outTab1
genename<-read.table("机器学习交集.txt",header = T,sep = "\t")
AAA<-intersect(genename[,1],rownames(Merge))
geneexpr<-Merge[AAA,]

write.table(geneexpr,"geneexpr.txt",quote = F,row.names =T ,sep = "\t")
save(geneexpr,file = "geneexpr.RDATA")
library(ConsensusClusterPlus)


#读取输入文件
data=geneexpr
data=as.matrix(data)

#聚类
maxK=9
results=ConsensusClusterPlus(data,
                             maxK=maxK,
                             reps=50,
                             pItem=0.8,
                             pFeature=1,
                             clusterAlg="pam",
                             distance="euclidean",
                             seed=123456,
                             plot="png")


#输出分型结果
clusterNum=2      #分几类，根据判断标准判断
cluster=results[[clusterNum]][["consensusClass"]]
cluster=as.data.frame(cluster)
colnames(cluster)=c("cluster")
letter=c("A","B","C","D","E","F","G")
uniqClu=levels(factor(cluster$cluster))
cluster$cluster=letter[match(cluster$cluster, uniqClu)]
clusterOut=rbind(ID=colnames(cluster), cluster)
write.table(clusterOut, file="Cluster.txt", sep="\t", quote=F, col.names=F)

