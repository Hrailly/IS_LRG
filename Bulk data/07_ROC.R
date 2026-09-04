rm(list = ls())
df <- read.table("mergeROC.txt",head=T,sep="\t",check.names = F)
head(df)
df<-df[,-1]


library("pROC")

#定义足够多的颜色，后面画线时从这里选颜色
mycol <- c("slateblue","seagreen3","dodgerblue","firebrick1","lightgoldenrod","magenta","orange2","grey","green","red")



pdf("ROC.pdf",height=6,width=6)
auc.out <- c()


#先画第一条线，此处是miRNA1
x <- plot.roc(df[,1],df[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, #绘制平滑曲线
              ci=TRUE, 
              main="",
              #print.thres="best", #把阈值写在图上，其sensitivity+ specificity之和最大
              col=mycol[2],#线的颜色
              lwd=2, #线的粗细
              legacy.axes=T)#采用大多数paper的画法，横坐标是“1-specificity”，从0到1

ci.lower <- round(as.numeric(x$ci[1]),3) #置信区间下限
ci.upper <- round(as.numeric(x$ci[3]),3) #置信区间上限

auc.ci <- c(colnames(df)[2],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
auc.out <- rbind(auc.out,auc.ci)


#再用循环画第二条和后面更多条曲线
for (i in 3:ncol(df)){
  x <- plot.roc(df[,1],df[,i],
                add=T, #向前面画的图里添加
                smooth=F,
                ci=TRUE,
                col=mycol[i],
                lwd=2,
                legacy.axes=T)
  
  ci.lower <- round(as.numeric(x$ci[1]),3)
  ci.upper <- round(as.numeric(x$ci[3]),3)
  
  auc.ci <- c(colnames(df)[i],round(as.numeric(x$auc),3),paste(ci.lower,ci.upper,sep="-"))
  auc.out <- rbind(auc.out,auc.ci)
}


# 输出AUC、AUC CI到文件
auc.out <- as.data.frame(auc.out)
colnames(auc.out) <- c("Name","AUC","AUC CI")



#绘制图例
legend.name <- paste(colnames(df)[2:length(df)],"AUC",auc.out$AUC,sep=" ")
legend("bottomright", 
       legend=legend.name,
       col = mycol[2:length(df)],
       lwd = 2,
       bty="n")
dev.off()










