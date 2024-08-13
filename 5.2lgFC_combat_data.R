load("C:/Users/Cycy/Desktop/diet_R/lcombat_FDRnotsel_remove36_NA80_sek.Rdata")
load("C:/Users/Cycy/Desktop/diet_R/sample_batch_list_remove36.Rdata")
library(dplyr)
library(RColorBrewer)
library(tidyverse)
df2=combat

#####################
#定义Label
label_diet <- clinical$diet[match(colnames(combat),clinical$ID)]
label_time <- clinical$time[match(colnames(combat),clinical$ID)]
label_bianhao <- clinical$bianhao[match(colnames(combat),clinical$ID)]

#OGTT-log2(foldchange)计算
theme="A"
df8 <- df2
df8[is.na(df8)]<-0

for (i in 2:24) {
  i=as.character(i)
Time_base <- which(label_bianhao==i&label_diet=="A"&label_time=="1")
Time_30min <- which(label_bianhao==i&label_diet=="A"&label_time=="2")
Time_1hour <- which(label_bianhao==i&label_diet=="A"&label_time=="3")
Time_2hours <- which(label_bianhao==i&label_diet=="A"&label_time=="4")
Time_3hours <- which(label_bianhao==i&label_diet=="A"&label_time=="5")
if(length(df8[,Time_base])!=0){
  df8$c1 <- apply(df8,1, function(x) log2((mean(x[Time_base],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
colnames(df8)[ncol(df8)]<-paste0("Time1_",theme,i)}

if(length(df8[,Time_30min])!=0&length(df8[,Time_base])!=0){df8$c2 <- apply(df8,1, function(x) log2((mean(x[Time_30min],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
colnames(df8)[ncol(df8)]<-paste0("Time2_",theme,i)}

if(length(df8[,Time_1hour])!=0&length(df8[,Time_base])!=0){df8$c3 <- apply(df8,1, function(x) log2((mean(x[Time_1hour],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
colnames(df8)[ncol(df8)]<-paste0("Time3_",theme,i)}

if(length(df8[,Time_2hours])!=0&length(df8[,Time_base])!=0){df8$c4 <- apply(df8,1, function(x) log2((mean(x[Time_2hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
colnames(df8)[ncol(df8)]<-paste0("Time4_",theme,i)}

if(length(df8[,Time_3hours])!=0&length(df8[,Time_base])!=0){df8$c5 <- apply(df8,1, function(x) log2((mean(x[Time_3hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
colnames(df8)[ncol(df8)]<-paste0("Time5_",theme,i)}
}

df9=df8[,526:ncol(df8)]

#Keto-log2(foldchange)计算
theme="B"
df8 <- df2
df8[is.na(df8)]<-0

for (i in 2:24) {
  i=as.character(i)
  Time_base <- which(label_bianhao==i&label_diet=="B"&label_time=="1")
  Time_30min <- which(label_bianhao==i&label_diet=="B"&label_time=="2")
  Time_1hour <- which(label_bianhao==i&label_diet=="B"&label_time=="3")
  Time_2hours <- which(label_bianhao==i&label_diet=="B"&label_time=="4")
  Time_3hours <- which(label_bianhao==i&label_diet=="B"&label_time=="5")
  if(length(df8[,Time_base])!=0){
    df8$c1 <- apply(df8,1, function(x) log2((mean(x[Time_base],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
    colnames(df8)[ncol(df8)]<-paste0("Time1_",theme,i)}
  
  if(length(df8[,Time_30min])!=0&length(df8[,Time_base])!=0){df8$c2 <- apply(df8,1, function(x) log2((mean(x[Time_30min],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time2_",theme,i)}
  
  if(length(df8[,Time_1hour])!=0&length(df8[,Time_base])!=0){df8$c3 <- apply(df8,1, function(x) log2((mean(x[Time_1hour],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time3_",theme,i)}
  
  if(length(df8[,Time_2hours])!=0&length(df8[,Time_base])!=0){df8$c4 <- apply(df8,1, function(x) log2((mean(x[Time_2hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time4_",theme,i)}
  
  if(length(df8[,Time_3hours])!=0&length(df8[,Time_base])!=0){df8$c5 <- apply(df8,1, function(x) log2((mean(x[Time_3hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time5_",theme,i)}
}

df10=df8[,526:ncol(df8)]


#HFD-log2(foldchange)计算
theme="C"
df8 <- df2
df8[is.na(df8)]<-0

for (i in 2:24) {
  i=as.character(i)
  Time_base <- which(label_bianhao==i&label_diet=="C"&label_time=="1")
  Time_30min <- which(label_bianhao==i&label_diet=="C"&label_time=="2")
  Time_1hour <- which(label_bianhao==i&label_diet=="C"&label_time=="3")
  Time_2hours <- which(label_bianhao==i&label_diet=="C"&label_time=="4")
  Time_3hours <- which(label_bianhao==i&label_diet=="C"&label_time=="5")
  if(length(df8[,Time_base])!=0){
    df8$c1 <- apply(df8,1, function(x) log2((mean(x[Time_base],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
    colnames(df8)[ncol(df8)]<-paste0("Time1_",theme,i)}
  
  if(length(df8[,Time_30min])!=0&length(df8[,Time_base])!=0){df8$c2 <- apply(df8,1, function(x) log2((mean(x[Time_30min],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time2_",theme,i)}
  
  if(length(df8[,Time_1hour])!=0&length(df8[,Time_base])!=0){df8$c3 <- apply(df8,1, function(x) log2((mean(x[Time_1hour],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time3_",theme,i)}
  
  if(length(df8[,Time_2hours])!=0&length(df8[,Time_base])!=0){df8$c4 <- apply(df8,1, function(x) log2((mean(x[Time_2hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time4_",theme,i)}
  
  if(length(df8[,Time_3hours])!=0&length(df8[,Time_base])!=0){df8$c5 <- apply(df8,1, function(x) log2((mean(x[Time_3hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time5_",theme,i)}
}

df11=df8[,526:ncol(df8)]


#Normal-log2(foldchange)计算
theme="D"
df8 <- df2
df8[is.na(df8)]<-0

for (i in 2:24) {
  i=as.character(i)
  Time_base <- which(label_bianhao==i&label_diet=="D"&label_time=="1")
  Time_30min <- which(label_bianhao==i&label_diet=="D"&label_time=="2")
  Time_1hour <- which(label_bianhao==i&label_diet=="D"&label_time=="3")
  Time_2hours <- which(label_bianhao==i&label_diet=="D"&label_time=="4")
  Time_3hours <- which(label_bianhao==i&label_diet=="D"&label_time=="5")
  if(length(df8[,Time_base])!=0){
    df8$c1 <- apply(df8,1, function(x) log2((mean(x[Time_base],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
    colnames(df8)[ncol(df8)]<-paste0("Time1_",theme,i)}
  
  if(length(df8[,Time_30min])!=0&length(df8[,Time_base])!=0){df8$c2 <- apply(df8,1, function(x) log2((mean(x[Time_30min],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time2_",theme,i)}
  
  if(length(df8[,Time_1hour])!=0&length(df8[,Time_base])!=0){df8$c3 <- apply(df8,1, function(x) log2((mean(x[Time_1hour],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time3_",theme,i)}
  
  if(length(df8[,Time_2hours])!=0&length(df8[,Time_base])!=0){df8$c4 <- apply(df8,1, function(x) log2((mean(x[Time_2hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time4_",theme,i)}
  
  if(length(df8[,Time_3hours])!=0&length(df8[,Time_base])!=0){df8$c5 <- apply(df8,1, function(x) log2((mean(x[Time_3hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time5_",theme,i)}
}

df12=df8[,526:ncol(df8)]

#LFD-log2(foldchange)计算
theme="E"
df8 <- df2
df8[is.na(df8)]<-0

for (i in 2:24) {
  i=as.character(i)
  Time_base <- which(label_bianhao==i&label_diet=="E"&label_time=="1")
  Time_30min <- which(label_bianhao==i&label_diet=="E"&label_time=="2")
  Time_1hour <- which(label_bianhao==i&label_diet=="E"&label_time=="3")
  Time_2hours <- which(label_bianhao==i&label_diet=="E"&label_time=="4")
  Time_3hours <- which(label_bianhao==i&label_diet=="E"&label_time=="5")
  if(length(df8[,Time_base])!=0){
    df8$c1 <- apply(df8,1, function(x) log2((mean(x[Time_base],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
    colnames(df8)[ncol(df8)]<-paste0("Time1_",theme,i)}
  
  if(length(df8[,Time_30min])!=0&length(df8[,Time_base])!=0){df8$c2 <- apply(df8,1, function(x) log2((mean(x[Time_30min],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time2_",theme,i)}
  
  if(length(df8[,Time_1hour])!=0&length(df8[,Time_base])!=0){df8$c3 <- apply(df8,1, function(x) log2((mean(x[Time_1hour],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time3_",theme,i)}
  
  if(length(df8[,Time_2hours])!=0&length(df8[,Time_base])!=0){df8$c4 <- apply(df8,1, function(x) log2((mean(x[Time_2hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time4_",theme,i)}
  
  if(length(df8[,Time_3hours])!=0&length(df8[,Time_base])!=0){df8$c5 <- apply(df8,1, function(x) log2((mean(x[Time_3hours],na.rm=T))/((mean(x[Time_base],na.rm=T)))))
  colnames(df8)[ncol(df8)]<-paste0("Time5_",theme,i)}
}

df13=df8[,526:ncol(df8)]

df14=cbind(df9,df10,df11,df12,df13)

lgFC_combat=df14
save(lgFC_combat,file = "lgFC_combat_FDRnotsel_remove36_NA80_sek.Rdata")

#lgFC_clinical表格
lgFC_clinical<-as.data.frame(colnames(lgFC_combat))
colnames(lgFC_clinical)="ID2"
lgFC_clinical$time=substr(lgFC_clinical$ID2,5,5)
lgFC_clinical$diet=substr(lgFC_clinical$ID2,7,7)
lgFC_clinical$bianhao=substr(lgFC_clinical$ID2,8,length(lgFC_clinical$ID2))
lgFC_clinical$weight=ifelse(lgFC_clinical$bianhao %in% c(6:13,16,18,23),"Normal Weight","Overweight or Obese")
lgFC_clinical$gender=ifelse(lgFC_clinical$bianhao %in% c(2,3,5,8,11,12,14,16,18,23,4),"Male","Female")
rownames(lgFC_clinical)=lgFC_clinical[,1]

save(lgFC_clinical,file="lgFC_sample_batch_list_remove36.Rdata")

