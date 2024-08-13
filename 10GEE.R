rm(list = ls())
library(readr)
library(plyr)
library(readxl)
library(stringr)
library(magrittr)
library(geepack)


load("C:/Users/ztycy/Desktop/diet_R/combat_FDRnotsel_remove36_NA80_sek.Rdata")
load("C:/Users/ztycy/Desktop/diet_R/sample_batch_list_remove36.Rdata")

info=clinical
df2=combat
df3<-data.frame(names(df2),t(df2))
names(df3)<-c("ID",as.matrix(names(df3[,-1])))
df4<-merge(info,df3,by.x = "ID")

#确定主题##############
theme="A"
gee<-data.frame(names(df4[,-c(1:8)]))
i=9
for (i in 9:891){
  df5=df4[,c(3,5,6,i)]
  df5 <- na.omit(df5)
  df6<-df5[order(df5[,'bianhao']),]
  geeglm3<-geeglm(df6[,4] ~ factor(time), 
                id = bianhao, 
                corstr = "independence", 
                family = "gaussian", 
                data = df6, 
                subset = (diet == theme))

  a<-as.data.frame((summary(geeglm3))[[6]])
  #n记得修改！！！
  gee_time2<-a$`Pr(>|W|)`[2]
  gee_time2<-as.numeric(gee_time2)
  padj_time2=p.adjust(gee_time2,method = 'BH',n=883)
  
  gee_time3<-a$`Pr(>|W|)`[3]
  gee_time3<-as.numeric(gee_time3)
  padj_time3=p.adjust(gee_time3,method = 'BH',n=883)
  
  gee_time4<-a$`Pr(>|W|)`[4]
  gee_time4<-as.numeric(gee_time4)
  padj_time4=p.adjust(gee_time4,method = 'BH',n=883)
  
  gee_time5<-a$`Pr(>|W|)`[5]
  gee_time5<-as.numeric(gee_time5)
  padj_time5=p.adjust(gee_time5,method = 'BH',n=883)

  gee[c(i-8),2]<-gee_time2
  gee[c(i-8),3]<-padj_time2
  
  gee[c(i-8),4]<-gee_time3
  gee[c(i-8),5]<-padj_time3
  
  gee[c(i-8),6]<-gee_time4
  gee[c(i-8),7]<-padj_time4
  
  gee[c(i-8),8]<-gee_time5
  gee[c(i-8),9]<-padj_time5
}
names(gee)<-c('protein',"pvalue_time2","padj_time2","pvalue_time3","padj_time3","pvalue_time4","padj_time4","pvalue_time5","padj_time5")
gee2<-data.frame(gee[c(which(gee$pvalue_time2<0.05|gee$pvalue_time3<0.05|gee$pvalue_time4<0.05|gee$pvalue_time5<0.05)),])
write.csv(gee2,paste0("Result_gee_all/",theme,"_total_gee_P_0.05.csv"),row.names = F)




