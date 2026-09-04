library(glmnet)
library(tidyverse)
rm(list = ls())
set.seed(1234)

data=read.csv("merge矩阵new.csv",header = T,sep = ",",row.names = 1)
gene=read.table("差异基因.txt",header = T,sep = "\t")
gene=gene[,1]
data=data[gene,]
data=as.data.frame(t(data))
data$ID=rownames(data)
clin=read.table("clinicaldata.txt",header = T,sep = "\t")
colnames(clin)[1]="ID"
aaa=intersect(clin$ID,data$ID)
clin=dplyr::filter(clin, ID %in% aaa)
data=dplyr::filter(data, ID %in% aaa)
data1=inner_join(clin,data,by = "ID")


rt=dplyr::select(data1,gene)
rownames(rt)=data1$ID


#构建模型
x=as.matrix(rt)
y=data1$condition

fit=glmnet(x, y, family = "binomial", alpha=1)
pdf(file="fit.pdf",width=6,height=5.5)
plot(fit)
dev.off()

cvfit=cv.glmnet(x, y, family="binomial", alpha=1,type.measure='deviance',nfolds = 10)
pdf(file="cvfit.pdf",width=6,height=5.5)
plot(cvfit)
dev.off()

#输出筛选的特征基因
coef=coef(fit, s = cvfit$lambda.min)
coef


#提取选中的基因名
active.min = which(coef != 0)-1
active.min = active.min[-1]
geneids <- colnames(x)[active.min]
geneids


#提取选中的基因对应的coefficient
index.min = coef[active.min+1]
index.min

combine <- cbind(geneids, index.min)#合并基因名和coef
colnames(combine)=c("ID","risk")
combine=as.data.frame(combine)
write.csv(combine,"risk.csv",row.names = F)



index=which(coef != 0)
lassoGene=row.names(coef)[index]
lassoGene=lassoGene[-1]
write.table(lassoGene, file="LASSO.gene.txt", sep="\t", quote=F, row.names=F, col.names=F)



#引用包
library(e1071)
library(kernlab)
library(caret)

set.seed(1234)


#读取输入文件
data=x
group=y

#SVM-RFE分析
Profile=rfe(x=data,
            y=as.numeric(as.factor(group)),
            sizes = c(2,4,6,8, seq(10,40,by=3)),
            rfeControl = rfeControl(functions = caretFuncs, method = "cv"),
            methods="svmRadial")

#绘制图形
pdf(file="SVM-RFE.pdf", width=6, height=5.5)
par(las=1)
x = Profile$results$Variables
y = Profile$results$RMSE
plot(x, y, xlab="Variables", ylab="RMSE (Cross-Validation)", col="darkgreen")
lines(x, y, col="darkgreen")
#标注交叉验证误差最小的点
wmin=which.min(y)
wmin.x=x[wmin]
wmin.y=y[wmin]
points(wmin.x, wmin.y, col="blue", pch=16)
text(wmin.x, wmin.y, paste0('N=',wmin.x), pos=2, col=2)
dev.off()

#输出选择的基因
featureGenes=Profile$optVariables
write.table(file="SVM-RFE.gene.txt", featureGenes, sep="\t", quote=F, row.names=F, col.names=F)

#随机森林
data1=arrange(data1,desc(condition))

group=data1$condition
data2=data1
rownames(data2)=data2$ID
data2=data2[,-1]
data2=data2[,-1]
colnames(data2)=str_replace_all(colnames(data2),"-",".")
library(randomForest)
library(ggpubr)

rf=randomForest(as.factor(group)~., data=data2, ntree=500)
pdf(file="forest.pdf", width=6, height=6)
plot(rf, main="Random forest", lwd=2)
dev.off()

#找出误差最小的点
optionTrees=which.min(rf$err.rate[,1])
optionTrees
rf2=randomForest(as.factor(group)~., data=data2, ntree=optionTrees)

#查看基因的重要性
#绘制基因的重要性图
importance=importance(x=rf2)
importance=as.data.frame(importance)
importance$size=rownames(importance)
importance=importance[,c(2,1)]
names(importance)=c("Gene","importance")
importance=arrange(importance,desc(importance))
#展示前10个基因的重要性
af=importance[1:10,]
p=ggdotchart(af, x = "Gene", y = "importance",
             color = "importance", # Custom color palette
             sorting = "descending",                       # Sort value in descending order
             add = "segments",                             # Add segments from y = 0 to dots
             add.params = list(color = "lightgray", size = 2), # Change segment color and size
             dot.size = 6,                        # Add mpg values as dot labels
             font.label = list(color = "white", size = 9,
                               vjust = 0.5),               # Adjust label parameters
             ggtheme = theme_bw()         ,               # ggplot2 theme
             rotate=TRUE                                       )#翻转坐标轴 
p1=p+ geom_hline(yintercept = 0, linetype = 2, color = "lightgray")+
  gradient_color(palette =c(ggsci::pal_npg()(2)[2],ggsci::pal_npg()(2)[1])      ) +#颜色
  grids()   
#保存图片
pdf(file="importance.pdf", width=6, height=8)
print(p1)
dev.off()
#挑选疾病特征基因
rfGenes=importance[order(importance[,"importance"], decreasing = TRUE),]
write.table(rfGenes, file="rfGenes.xls", sep="\t", quote=F, col.names=T, row.names=F)




