library(dplyr) 
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(ggsci)
library(clusterProfiler)
library(org.Hs.eg.db)
#BiocManager::install("rlang")
#BiocManager::install("org.Hs.eg.db")


Proj="New studies-dotline graph"
dirname=paste0("Fig-202408-",Proj)
dir.create(dirname)
#导入数据(lgFC版本)
load("D:/Cycy1/R_Proj/Diet/diet_R/lgFC_combat_FDRnotsel_remove36_NA80_sek.Rdata")
load("D:/Cycy1/R_Proj/Diet/diet_R/lgFC_sample_batch_list_remove36.Rdata")

combat=lgFC_combat
clinical=lgFC_clinical

#加SYMBOL列
s2e <- bitr(rownames(combat), 
            fromType = "UNIPROT",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)#人类

# 表达矩阵行名替换
combat2=combat
combat2$UNIPROT=rownames(combat2)
combat3 <- merge(combat2,s2e)#两个数据
combat4=combat3[!duplicated(combat3$SYMBOL),]
combat5=combat4[,-1]
rownames(combat5)=combat5$SYMBOL

#选择差异基因！！！！
target_protein=c("PPBP","SRGN","MIF","SOD1","IDH1","TXN","FTL","PF4V1","RAP1B","LMNA","FLNA","ARPC2","FTH1")
s2e=s2e[s2e$SYMBOL%in%target_protein,]
colnames(s2e)=c("protein","SYMBOL")

###画图
data=combat5
for (i in 1:nrow(s2e)){
  protein<-s2e$protein[i]
  symbol=s2e$SYMBOL[i]
  #修改文件名！！！！！！
  theme7<-paste0(protein,"_",symbol)
  data2=data[symbol,]
  
  data3=as.data.frame(t(data2))
  data5=data3
  data5$ID2=row.names(data5)
  df1<-merge(clinical,data5)
  df2=df1[,c(3,2,7)]
  colnames(df2)=c("Diet","Time","log2FC")
  df2=na.omit(df2)
  df2$Time=as.numeric(df2$Time)
  df2$log2FC=as.numeric(df2$log2FC)
  
  df2$Diet[which(df2$Diet=="A")]<-"OGTT"
  df2$Diet[which(df2$Diet=="B")]<-"Ketogenic Breakfast"
  df2$Diet[which(df2$Diet=="C")]<-"High Fat Breakfast"
  df2$Diet[which(df2$Diet=="D")]<-"Normal Breakfast"
  df2$Diet[which(df2$Diet=="E")]<-"Low Fat Breakfast"
  
  df2$Time[which(df2$Time=="1")]<-"0"
  df2$Time[which(df2$Time=="2")]<-"30"
  df2$Time[which(df2$Time=="3")]<-"60"
  df2$Time[which(df2$Time=="4")]<-"120"
  df2$Time[which(df2$Time=="5")]<-"180"
  
  a<-ggline(df2, x = "Time", y = "log2FC", 
            add = c("mean_se"),
            order = c("0","30","60","120","180"),size=1.5,
            color = "Diet", palette = c("orange","#DC0000B2","#00A087B2","#4DBBD5B2","#3C5488FF" ),title=theme7)
  df2_filtered <- df2[df2$Time != "0", ]
  b<-a+ stat_compare_means(data = df2_filtered, aes(group = Diet), label = "p.signif", label.y = 2, size = 8)
  c<-ggpar(b,ylim=c(-2.5,2.5),
           xlab = "Time (min)",
           legend = "right",
           yticks.by=1,
           font.title = c(30,"plain","black"),
           font.caption = c(22,"plain","black"),
           font.x = c(22,"plain","black"),
           font.y = c(22, "plain", "black"),
           font.tickslab = c(22, "plain", "black"),
           font.legend =  c(22,"plain","black"))
  ggsave(paste0(dirname,"/",theme7,"-2.pdf"),plot=c,width = 10, height =5)}
