rm(list = ls())
#引用包
library(limma)
library(sva)
library(tidyverse)
library(limma)
library(ggrepel)
library(ggthemes)
library(tidyverse)
library(pheatmap)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


#均一化处理
express <- read.csv("GSE37587矩阵.csv", row.names = 1)
group_list <- read.table("GSE37587样本.txt", sep = "\t", header = T)
express<-dplyr::select(express,group_list$Accession)
# 针对复杂基因名的状况
# express$ID=NA
# for (i in 1:nrow(express)) {
#   express$ID[i]<-unlist(str_split_fixed(rownames(express)[i],"//",3))[2]
# }
# express$ID<-str_replace_all(express$ID," ","")
# express<-distinct(express,ID,.keep_all = T)
# rownames(express)<-express$ID
# express<-dplyr::select(express,-ID)


#去掉一个探针对应多个基因的结果
express<-express[which(unlist(str_split_fixed(rownames(express),"///",2))[,2]==""),]
express=na.omit(express)
range(express)
# if yes no need to log transfer, if above this range, have to do log transfer.
express <- log2(express+1)
range(express)
express=rbind(geneNames=colnames(express), express)
write.table(express, file="GSE37587.txt", sep="\t", quote=F, col.names=F)

rm(list = ls())


files=c("GSE16561.txt","GSE22255.txt" ,"GSE37587.txt")       

#获取交集基因
geneList=list()
for(i in 1:length(files)){
  inputFile=files[i]
  rt=read.table(inputFile, header=T, sep="\t",check.names=F)
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  geneList[[header[1]]]=as.vector(rt[,1])
}
intersectGenes=Reduce(intersect, geneList)

#数据合并
allTab=data.frame()
batchType=c()
for(i in 1:length(files)){
  inputFile=files[i]
  header=unlist(strsplit(inputFile, "\\.|\\-"))
  #读取输入文件，并对输入文件进行整理
  rt=read.table(inputFile, header=T, sep="\t", check.names=F)
  rt=as.matrix(rt)
  rownames(rt)=rt[,1]
  exp=rt[,2:ncol(rt)]
  dimnames=list(rownames(exp),colnames(exp))
  data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
  rt=avereps(data)
  colnames(rt)=paste0(header[1], "_", colnames(rt))

  #对数值大的数据取log2
  qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
  LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
  if(LogC){
    rt[rt<0]=0
    rt=log2(rt+1)}
  if(header[1] != "TCGA"){
    rt=normalizeBetweenArrays(rt)
  }
  #数据合并
  if(i==1){
    allTab=rt[intersectGenes,]
  }else{
    allTab=cbind(allTab, rt[intersectGenes,])
  }
  batchType=c(batchType, rep(i,ncol(rt)))
}

#对数据进行矫正，输出矫正后的结果
outTab=ComBat(allTab, batchType, par.prior=TRUE)

library(FactoMineR)
library(factoextra)
ddb.pca <- PCA(t(allTab), graph = FALSE)

pheno<-data.frame(ID=colnames(allTab))
pheno$cancer<-pheno$ID
pheno[1:63,2]<-"GSE16561"
pheno[64:103,2]<-"GSE22255"
pheno[104:171,2]<-"GSE37587"

pdf(file = "去除批次效应前.pdf",height=8,width=8)
fviz_pca_ind(ddb.pca,
             geom.ind = "point", # 只显示点
             pointsize =2, # 点的大小
             pointshape = 21, # 点的形状
             fill.ind = pheno$cancer, # 分组颜色
             palette = "lacent", # c("#00AFBB", "#E7B800", "#FC4E07")
             addEllipses = TRUE, # 增加置信椭圆
             legend.title = "Groups", # 图例标题
             title="") +
  theme_bw() + # 和ggplot2对接进行美化
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )
dev.off()



library(FactoMineR)
library(factoextra)
ddb.pca <- PCA(t(outTab), graph = FALSE)
pdf(file = "去除批次效应后.pdf",height=8,width=8)
fviz_pca_ind(ddb.pca,
             geom.ind = "point", # 只显示点
             pointsize =2, # 点的大小
             pointshape = 21, # 点的形状
             fill.ind = pheno$cancer, # 分组颜色
             palette = "lacent", # c("#00AFBB", "#E7B800", "#FC4E07")
             addEllipses = TRUE, # 增加置信椭圆
             legend.title = "Groups", # 图例标题
             title="") +
  theme_bw() + # 和ggplot2对接进行美化
  theme(text=element_text(size=14,face="plain",color="black"),
        axis.title=element_text(size=16,face="plain",color="black"),
        axis.text = element_text(size=14,face="plain",color="black"),
        legend.title = element_text(size=16,face="plain",color="black"),
        legend.text = element_text(size=14,face="plain",color="black"),
        legend.background = element_blank(),
        legend.position=c(0.9,0.1)
  )

dev.off()


library(data.table)
clin=fread("clinicaldata.txt",header = T)
clin=dplyr::filter(clin,condition == "IS")
clin=clin$Accession
colnames(outTab)=str_replace_all(colnames(outTab),"GSE16561_","")
colnames(outTab)=str_replace_all(colnames(outTab),"GSE37587_","")
colnames(outTab)=str_replace_all(colnames(outTab),"GSE22255_","")
outTab1=outTab[,clin]
outTab1=as.data.frame(outTab1)
outTab=as.data.frame(outTab)

save(outTab,file = "merge_all.RDATA")
outTab=rbind(geneNames=colnames(outTab), outTab)
write.table(outTab, file="merge_all.txt", sep="\t", quote=F, col.names=F)

save(outTab1,file = "merge_IS.RDATA")
outTab1=rbind(geneNames=colnames(outTab1), outTab1)
write.table(outTab1, file="merge_IS.txt", sep="\t", quote=F, col.names=F)

