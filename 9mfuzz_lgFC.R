rm(list = ls())
library(readr)
library(plyr)
library(readxl)
library(stringr)
library(magrittr)
#BiocManager::install("Mfuzz")
library(Mfuzz)
source("datamining_library_ge20200306.R")

load("C:/Users/Cycy/Desktop/diet_R/lgFC_combat_FDRnotsel_remove36_NA80_sek.Rdata")
load("C:/Users/Cycy/Desktop/diet_R/lgFC_sample_batch_list_remove36.Rdata")

info=lgFC_clinical
df2=lgFC_combat
#mfuzz################
df3<-data.frame(names(df2),t(df2))
names(df3)<-c("ID2",as.matrix(names(df3[,-1])))
df4<-merge(info,df3,by.x = "ID2")
#选择饮食A
df5<-df4[which(df4$diet=="E"),]
#分组求均值
df6<-aggregate(df5[,colnames(df5)[7:ncol(df5)]],by=list(df5$time),mean,na.rm= TRUE)
row.names(df6)<-c(as.matrix((df6[,1])))
df7<-data.frame(t(df6[,-1]))

set.seed(2021)
a<-ge.mfuzz.cselection(df7,range = seq(3,10,1))
b<-ge.mfuzz.getresult(a,4,"lgFC_E_Time_protein_mfuzz_ge2021type_4")

#anova################
rm(list = ls())
library(readr)
library(plyr)
library(readxl)
library(stringr)
library(magrittr)
source("datamining_library_ge20200306.R")

load("C:/Users/Cycy/Desktop/diet_R/lgFC_combat_FDRnotsel_remove36_NA80_sek.Rdata")
load("C:/Users/Cycy/Desktop/diet_R/lgFC_sample_batch_list_remove36.Rdata")
info=lgFC_clinical
df2=lgFC_combat
df3<-data.frame(names(df2),t(df2))
names(df3)<-c("ID2",as.matrix(names(df3[,-1])))
df4<-merge(info,df3,by.x = "ID2")
df5<-df4[which(df4$diet=="E"),]


Matrix<-df5
for (K in 1:4) {
  
  mfuzz_prot<-read.csv(paste0("lgFC_E_Time_protein_mfuzz_ge2021type_4/mfuzz_",K,".csv"))
  prot<-as.matrix(mfuzz_prot$X)
  
  Matrix2<-data.frame(Matrix[,6],Matrix[,c(which(names(Matrix)%in% prot))])
  
  col<-ncol(Matrix2)
  names(Matrix2)<-c("label",as.matrix(names(Matrix2[,-1])))
  anova<-data.frame(names(Matrix2[,-1]))
  for (i in 2:col) {
    aov<-(summary(aov(Matrix2[,i] ~ Matrix2[,1],Matrix2))[[1]])$`Pr(>F)`[1]
    aov<-as.numeric(aov)
    padj=p.adjust(aov,method = 'BH',n=col)
    # aov<-as.numeric(aov,digits = 4, scientific = F)
    anova[c(i-1),2]<-aov
    anova[c(i-1),3]<-padj
  }
  names(anova)<-c('protein',"pvalue","padj")
  #因为padj太大，所以选了pvalue
  anova2<-data.frame(anova[c(which(anova$pvalue<=0.05)),])
  write.csv(anova2,paste0("lgFC_E_Time_protein_mfuzz_ge2021type_4/",K,"_mfuzz_anova_0.05.csv"),row.names = F)
}
