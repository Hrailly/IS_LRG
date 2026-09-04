rm(list = ls())
library(circlize)
library(ggsci)
library(parallel)
library(tidyverse)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor

# 计算相关系数的函数
genecor.parallel <- function(data,gene,cl){
  cl <- makeCluster(cl)
  y <- as.numeric(data[gene,])
  rownames <- rownames(data)
  dataframes <- do.call(rbind, parLapply(cl=cl,rownames, function(x){
    dd  <- cor.test(as.numeric(data[x,]), y, type="spearman")
    data.frame(Gene_1=gene, Gene_2=x, cor=dd$estimate, p.value=dd$p.value)
  }))
  stopCluster(cl)
  return(dataframes)
}

# 画图的函数
genecor_circleplot <- function(x){
  Corr <- data.frame(rbind(data.frame(Gene=x[,1], Correlation=x[,3]), 
                           data.frame(Gene=x[,2], Correlation=x[,3])), stringsAsFactors = F)      
  Corr$Index <- seq(1,nrow(Corr),1) #记录基因的原始排序，记录到Index列
  Corr <- Corr[order(Corr[,1]),] #按照基因名排序
  corrsp <- split(Corr,Corr$Gene)
  corrspe <- lapply(corrsp, function(x){x$Gene_Start<-0
  
  #依次计算每个基因的相关系数总和，作为基因终止位点
  if (nrow(x)==1){x$Gene_End<-1}else{
    x$Gene_End<-sum(abs(x$Correlation))} 
  x})
  GeneID <- do.call(rbind,corrspe)
  GeneID <- GeneID[!duplicated(GeneID$Gene),]
  
  #基因配色
  mycol <- pal_d3("category20c")(20)
  n <- nrow(GeneID)
  GeneID$Color <- mycol[1:n]
  
  #连线的宽度是相关系数的绝对值
  Corr[,2] <- abs(Corr[,2]) 
  corrsl <- split(Corr,Corr$Gene)
  aaaaa <- c()
  corrspl <- lapply(corrsl,function(x){nn<-nrow(x)
  for (i in 1:nn){
    aaaaa[1] <- 0
    aaaaa[i+1] <- x$Correlation[i]+aaaaa[i]}
  bbbbb <- data.frame(V4=aaaaa[1:nn],V5=aaaaa[2:(nn+1)])
  bbbbbb <- cbind(x,bbbbb)
  bbbbbb
  })
  Corr <- do.call(rbind,corrspl)
  
  #根据Index列，把基因恢复到原始排序
  Corr <- Corr[order(Corr$Index),]
  
  #V4是起始位置，V5是终止位置
  #把它写入Links里，start_1和end_1对应Gene_1，start_2和end_2对应Gene_2
  x$start_1 <- Corr$V4[1:(nrow(Corr)/2)]
  x$end_1 <- Corr$V5[1:(nrow(Corr)/2)]
  x$start_2 <- Corr$V4[(nrow(Corr)/2 + 1):nrow(Corr)]
  x$end_2 <- Corr$V5[(nrow(Corr)/2 + 1):nrow(Corr)]
  
  #连线（相关系数）的配色
  #相关系数最大为1，最小-1，此处设置201个颜色
  #-1到0就是前100，0到1就是后100
  color <- data.frame(colorRampPalette(c("#67BE54", "#FFFFFF", "#F82C2B"))(201))
  #根据相关系数的数值，给出相应的颜色
  for (i in 1:nrow(x)){
    x[i,8] <- substring(color[x[i,3] * 100 + 101, 1], 1, 7)
  }
  names(x)[8] <- "color"
  
  #绘图区设置
  #par(mar=rep(0,4))
  circos.clear()
  circos.par(start.degree = 90, #从哪里开始画，沿着逆时针顺序
             gap.degree = 5, #基因bar之间的间隔大小
             track.margin = c(0,0.23), #值越大，基因跟连线的间隔越小
             cell.padding = c(0,0,0,0)
  )
  circos.initialize(factors = GeneID$Gene,
                    xlim = cbind(GeneID$Gene_Start, GeneID$Gene_End))
  
  #先画基因
  circos.trackPlotRegion(ylim = c(0, 1), factors = GeneID$Gene, 
                         track.height = 0.05, #基因线条的胖瘦
                         panel.fun = function(x, y) {
                           name = get.cell.meta.data("sector.index") 
                           i = get.cell.meta.data("sector.numeric.index") 
                           xlim = get.cell.meta.data("xlim")
                           ylim = get.cell.meta.data("ylim")
                           circos.text(x = mean(xlim), y = 1,
                                       labels = name,
                                       cex = 1, #基因ID文字大小
                                       niceFacing = TRUE, #保持基因名的头朝上
                                       facing = "bending", #基因名沿着圆弧方向，还可以是reverse.clockwise
                                       adj = c(0.5, -2.8), #基因名所在位置，分别控制左右和上下
                                       font = 2 #加粗
                           )
                           circos.rect(xleft = xlim[1], 
                                       ybottom = ylim[1],
                                       xright = xlim[2], 
                                       ytop = ylim[2],
                                       col = GeneID$Color[i],
                                       border = GeneID$Color[i])
                           
                           circos.axis(labels.cex = 0.7, 
                                       direction = "outside"
                           )})
  
  #画连线
  for(i in 1:nrow(x)){
    circos.link(sector.index1 = x$Gene_1[i], 
                point1 = c(x[i, 4], x[i, 5]),
                sector.index2 = x$Gene_2[i], 
                point2 = c(x[i, 6], x[i, 7]),
                col = paste(x$color[i], "C9", sep = ""), 
                border = FALSE, 
                rou = 0.7
    )}
  
  #画图例
  i <- seq(0,0.995,0.005)
  rect(-1+i/2, #xleft
       -1, #ybottom
       -0.9975+i/2, #xright
       -0.96, #ytop
       col = paste(as.character(color[,1]), "FF", sep = ""),
       border = paste(as.character(color[,1]), "FF", sep = ""))
  text(-0.97, -1.03, "-1")
  text(-0.51, -1.03, "1")
}

inputtemp <- read.csv("merge矩阵new.csv", header = T,sep = "," ,row.names = 1)
aaa<-read.table("机器学习交集.txt",header = T)
aaa<-aaa[,1]
inputtemp<-as.data.frame(t(inputtemp))
inputtemp<-dplyr::select(inputtemp,aaa)
clin<-read.table("clinicaldata.txt",header = T)
inputtemp=inputtemp[clin$Accession,]

ccc<-inputtemp
ddd<-clin
colnames(ddd)[1]="ID"
ccc$ID<-rownames(ccc)
eee<-inner_join(ddd,ccc,by="ID")
write.table(eee,file = "mergeROC.txt",row.names = F,sep = "\t")

genecorl <- lapply(colnames(inputtemp),function(x){
  ddd <- genecor.parallel(data = t(inputtemp), cl=1, gene=x) #一定要注意cl参数根据自己电脑cpu线程调整
  ddd  
})
genecor <- do.call(rbind, genecorl)

# 删掉p value = 0的，也就是自己跟自己配对
genecorr <- genecor[-which(genecor$p.value==0),]

# 删掉A vs B 和 B vs A其中一个
genecorrr<-genecorr[!duplicated(genecorr$cor),]

# 保存到文件
genecorrr$p.value <- NULL
genecor_circleplot(genecorrr)

pdf("mergecorrelation.pdf", width = 5, height = 5)
genecor_circleplot(genecorrr)
dev.off()

drawdata=eee

library(tidyverse)
library(corrplot)
library(circlize)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor


inputtemp1 <- read.table("merge_ssGSEA.result.txt", header=T,sep = "\t" )
rownames(inputtemp1)=inputtemp1$id
rownames(inputtemp1)=str_replace_all(rownames(inputtemp1),"na","")
inputtemp1=inputtemp1[,-1]
inputtemp1=as.data.frame(t(inputtemp1))
inputtemp1$ID=rownames(inputtemp1)

inputtemp$ID<-rownames(inputtemp)
drawdata<-inner_join(inputtemp1,inputtemp,by="ID")




gene <- "HIST2H2BE"
y <- as.numeric(drawdata[,gene])

drawdata1<-select(drawdata,2:23)

### 第1,写出单次处理的function
mycor = function(x){
  dd = cor.test(as.numeric(drawdata1[, x]),y,method ="spearman",exact=FALSE)
  data.frame(cell=x,cor=dd$estimate,p.value=dd$p.value)
}

colnames(drawdata1)<-gsub("\\."," ",colnames(drawdata1))


### 第2步lapply批量作用于函数，返回list
lapplylist = lapply(colnames(drawdata1),mycor)

### 第3步do.call 转换list
cor_data <- do.call(rbind,lapplylist)

cor_data %>% 
  filter(p.value <0.05) %>% 
  ggplot(aes(cor,forcats::fct_reorder(cell,cor)))+
  geom_segment(aes(xend=0,yend=cell))+
  geom_point(aes(col=p.value,size=abs(cor)))+
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3"))+
  #scale_color_viridis_c(begin = 0.5, end = 1)+
  scale_size_continuous(range =c(2,8))+
  theme_bw()+
  ylab(NULL)+
  xlab(gene)

write.csv(cor_data,file = paste0(gene,"相关性分析结果.csv"))
pdf(file=paste0(gene,"免疫细胞相关性分析.pdf"),width=5,height=7)
cor_data %>% 
  filter(p.value <0.05) %>% 
  ggplot(aes(cor,forcats::fct_reorder(cell,cor)))+
  geom_segment(aes(xend=0,yend=cell))+
  geom_point(aes(col=p.value,size=abs(cor)))+
  scale_colour_gradientn(colours=c("#7fc97f","#984ea3"))+
  #scale_color_viridis_c(begin = 0.5, end = 1)+
  scale_size_continuous(range =c(2,8))+
  theme_bw()+
  ylab(NULL)+
  xlab(gene)
dev.off()


