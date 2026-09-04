library(tidyverse)
library(ggpubr)
library(survminer)
library(survival) 
library(export)
library(survivalROC)
library(pheatmap)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(ReactomePA)
library(igraph)
library(ggraph)
library(ggradar)
rm(list = ls())
FF<-read.csv("merge矩阵new.csv",sep = ",",header = T)
rownames(FF)<-FF[,1]
FF<-FF[,-1]
FF<-as.data.frame(t(FF))



gene = "HIST2H2BE"
y <- as.numeric(FF[,gene])#开始相关性分析
colnames <- colnames(FF)
cor_data_df <- data.frame(colnames)
for (i in 1:length(colnames(FF))){
  test <- cor.test(as.numeric(FF[,i]),y,type="pearson")
  cor_data_df[i,2] <- test$estimate
  cor_data_df[i,3] <- test$p.value
}
names(cor_data_df) <- c("symbol","correlation","pvalue")
write.csv (cor_data_df, file =paste(gene,'批量相关性分析.csv',sep=' '), row.names =FALSE)#将文件导出


#相关性分析结果排序选取P<0.05,top50的做相关性热图，R>0.3的做相关性图
cor_data_df<-na.omit(cor_data_df)
pos<-cor_data_df %>%
  filter(pvalue < 0.05 ) %>%
  arrange(desc(correlation)) %>%
  top_n(50,correlation)
neg<-cor_data_df %>%
  filter(pvalue < 0.05 ) %>%
  arrange(correlation) %>%
  top_n(-49,correlation) 
#正相关热图
phemapdata<-FF %>%
  dplyr::select(pos$symbol)
phemapdata<-arrange(phemapdata,desc(phemapdata[,gene]))

phemapdata<-as.data.frame(t(phemapdata))
p<-pheatmap(phemapdata,scale="row",cluster_row=F,cluster_col=F,legend= T,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
pdf(file = paste(gene,'正相关性热图.pdf',sep=' '),width = 8,height = 6)
print(p)
dev.off()
#负相关热图
phemapdata<-FF %>%
  dplyr::select(gene,neg$symbol)
phemapdata<-arrange(phemapdata,desc(phemapdata[,gene]))


phemapdata<-as.data.frame(t(phemapdata))
p<-pheatmap(phemapdata,scale="row",cluster_row=F,cluster_col=F,legend= T,show_colnames = F,color = colorRampPalette(c("blue", "white", "red"))(50))
dev.off()
pdf(file = paste(gene,'负相关性热图.pdf',sep=' '),width = 8,height = 6)
print(p)
dev.off()

FUJI<-cor_data_df %>%
  filter(pvalue < 0.05 ) %>%
  arrange(desc(correlation)) %>%
  top_n(500,correlation)

genename <- as.character(FUJI[,1]) #提取第一列基因名
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
genelist_input<-gene_map[,2]
genelist_input<-na.omit(genelist_input)


#GO分析
#BP
Go_result_BP <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="BP", pvalueCutoff=1) #基因ID类型为ENSEMBL的ID形式，选择BP功能组，以P值0.05为界限
p1<-dotplot(Go_result_BP, showCategory=20)
p2<-barplot(Go_result_BP, showCategory=20)
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_x_discrete(labels=function(x) str_wrap(x, width=50))
y=Go_result_BP
yy<-as.data.frame(y)
write.csv(as.data.frame(y),paste(gene,"GO-BP.csv",sep = " "),row.names =F)#导出结果至默认路径下。
if(length(rownames(yy))!= 0) {
  ggsave(p1,filename = paste(gene,'富集分析GO_BP气泡图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p1,filename = paste(gene,'富集分析GO_BP气泡图.TIFF',sep=' '),width = 7.23,height = 8)
}
#CC
Go_result_CC <- enrichGO(genelist_input, 'org.Hs.eg.db', ont="CC", pvalueCutoff=1) #基因ID类型为ENSEMBL的ID形式，选择CC功能组，以P值0.05为界限
p1<-dotplot(Go_result_CC, showCategory=20) #气泡图，显示前二十个
p2<-barplot(Go_result_CC, showCategory=20) #条形图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_x_discrete(labels=function(x) str_wrap(x, width=50))
y=Go_result_CC
yy<-as.data.frame(y)
write.csv(as.data.frame(y),paste(gene,"GO-CC.csv",sep = " "),row.names =F)#导出结果至默认路径下
if(length(rownames(yy))!= 0) {
  ggsave(p1,filename = paste(gene,'富集分析GO_CC气泡图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p1,filename = paste(gene,'富集分析GO_CC气泡图.TIFF',sep=' '),width = 7.23,height = 8)
}
#MF
Go_result_MF <- enrichGO(genelist_input, 'org.Hs.eg.db',ont="MF", pvalueCutoff=1) #基因ID类型为ENSEMBL的ID形式，选择MF功能组，以P值0.05为界限
p1<-dotplot(Go_result_MF, showCategory=20) #气泡图，显示前二十个
p2<-barplot(Go_result_MF, showCategory=20) #条形图，显示前二十个
p1<-p1 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_x_discrete(labels=function(x) str_wrap(x, width=50))
y=Go_result_MF
write.csv(as.data.frame(y),file = paste(gene,"GO-MF.csv",sep = " "),row.names =F)#导出结果至默认路径下
yy<-as.data.frame(y)
if(length(rownames(yy))!= 0){
  ggsave(p1,filename = paste(gene,'富集分析GO_MF气泡图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p1,filename = paste(gene,'富集分析GO_MF气泡图.TIFF',sep=' '),width = 7.23,height = 8)
}
#KEGG分析
KEGG_result <- enrichKEGG(genelist_input, keyType = "kegg",pvalueCutoff=1,qvalueCutoff=1,pAdjustMethod = "BH", minGSSize = 5, maxGSSize = 500,organism = "hsa", use_internal_data=T)  #KEGG富集分析

p1<-barplot(KEGG_result, showCategory=20)#绘制条形图
p2<-dotplot(KEGG_result, showCategory=20) #气泡图，显示前二十个
p1<-p1 + scale_x_discrete(labels=function(x) str_wrap(x, width=50))
p2<-p2 + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
y=KEGG_result
yy<-as.data.frame(y)
write.csv(as.data.frame(y), file=paste(gene,"KEGG.csv",sep = " "),row.names =F)#导出结果至默认路径下
if(length(rownames(yy))!= 0){
  ggsave(p2,filename = paste(gene,'富集分析KEGG气泡图.pdf',sep=' '),width = 7.23,height = 8)
  ggsave(p2,filename = paste(gene,'富集分析KEGG气泡图.TIFF',sep=' '),width = 7.23,height = 8)
}
#GSEA分析

GSEAdata<-cor_data_df %>%
  filter(pvalue < 0.05 ) %>%
  arrange(desc(correlation)) %>%
  dplyr::select(symbol,correlation)
genename <- as.character(GSEAdata[,1]) #提取第一列基因名
gene_map <- biomaRt::select(org.Hs.eg.db, keys=genename, keytype="SYMBOL", columns=c("ENTREZID"))
colnames(gene_map)[1]<- 'symbol'
genelist_input<-gene_map %>%
  na.omit() %>%
  inner_join(GSEAdata,by= 'symbol') %>%
  dplyr::select(ENTREZID,correlation)
geneList = genelist_input[,2]
names(geneList) = as.character(genelist_input[,1])
geneList = sort(geneList, decreasing = TRUE)

#GSEA分析——GO
Go_gseresult <- gseGO(geneList, org.Hs.eg.db, keyType = "ENTREZID", ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)
#GSEA分析——KEGG
KEGG_gseresult <- gseKEGG(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1,use_internal_data=T)
#GSEA分析——Reactome
Go_Reactomeresult <- gsePathway(geneList, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=1)

write.csv (Go_gseresult,     file = paste(gene,'GSEA_GO.csv',sep=' '), row.names =F)
save(Go_gseresult,file = paste(gene,'GSEA_GO.RDATA',sep=' '))
write.csv (KEGG_gseresult,   file = paste(gene,'GSEA_KEGG.csv',sep=' '), row.names =F)
save(KEGG_gseresult,file = paste(gene,'GSEA_KEGG.RDATA',sep=' '))
write.csv (Go_Reactomeresult,file = paste(gene,'GSEA_Reactome.csv',sep=' '), row.names =F)
save(Go_Reactomeresult,file = paste(gene,'GSEA_Reactome.RDATA',sep=' '))
#波浪图
p<-ridgeplot(Go_gseresult,20) #输出前十个结果
p<-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename = paste(gene,'GSEA_GO波浪图.pdf',sep=' '),width = 8,height = 8)
ggsave(p,filename = paste(gene,'GSEA_GO波浪图.TIFF',sep=' '),width = 8,height = 8)
p<-ridgeplot(KEGG_gseresult, 20) #输出前十个结果
p<-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename = paste(gene,'GSEA_KEGG波浪图.pdf',sep=' '),width = 8,height = 8)
ggsave(p,filename = paste(gene,'GSEA_KEGG波浪图.TIFF',sep=' '),width = 8,height = 8)
p<-ridgeplot(Go_Reactomeresult, 20) #输出前十个结果
p<-p + scale_y_discrete(labels=function(x) str_wrap(x, width=50))
ggsave(p,filename = paste(gene,'GSEA_Reactomeresult波浪图.pdf',sep=' '),width = 8,height = 8)
ggsave(p,filename = paste(gene,'GSEA_Reactomeresult波浪图.TIFF',sep=' '),width = 8,height = 8)
