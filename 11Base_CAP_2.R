load("C:/Users/ztycy/Desktop/diet_R/combat_FDRnotsel_remove36_NA80_sek.Rdata")
load("C:/Users/ztycy/Desktop/diet_R/sample_batch_list_remove36.Rdata")
lab<-read.csv("lab.csv")


theme1="0.26"
ogtt<-read.csv(paste0("Result_gee_lgFC/A_gee_lgFC",theme1,"_sel.csv"))
Keto<-read.csv(paste0("Result_gee_lgFC/B_gee_lgFC",theme1,"_sel.csv"))
HFD<-read.csv(paste0("Result_gee_lgFC/C_gee_lgFC",theme1,"_sel.csv"))
Normal<-read.csv(paste0("Result_gee_lgFC/D_gee_lgFC",theme1,"_sel.csv"))
LFD<-read.csv(paste0("Result_gee_lgFC/E_gee_lgFC",theme1,"_sel.csv"))




#Ttest####################
#所有批次的空腹

theme2="1"
#定义A-1
label_diet <- clinical$diet[match(colnames(combat),clinical$ID)]
label_time <- clinical$time[match(colnames(combat),clinical$ID)]

data=combat[,which(label_time==theme2)]
#合并clinical和lab表格
clinical2=clinical[clinical$time==theme2,]
lab2=merge(lab,clinical2,by="bianhao")

#实验组和对照组
theme3="cap238"
library(tidyverse)
lab3=lab2 %>%  mutate(cap_group = if_else(lab2$cap>238,true = "high",false = "normal"))


label_capgroup <- lab3$cap_group[match(colnames(data),lab3$ID)]


sy_type <- which(label_capgroup=="high")
dz_type <- which(label_capgroup=="normal")

#foldchange
df8 <- data
df8[is.na(df8)]<-0
df8$fd <- apply(df8,1, function(x) log2((mean(x[sy_type],na.rm = T)/mean(x[dz_type],na.rm = T))))
x<-c(0.0,0.0)


#pvalue and adjust
df9 <- data
df9[is.na(df9)]<-0
for(i in 1:nrow(df9)){
  x[i] <- t.test(df9[i,sy_type],df9[i,dz_type], paired = F, var.equal = F)$p.value}
df8$P_value<-x
df8$P_value_adjust<-p.adjust(df8$P_value, method="BH")
dep<-df8[,(ncol(df8)-2):ncol(df8)]

logFC_t=0.26
P.Value_t = 0.05
k1 = (dep$P_value< P.Value_t)&(dep$fd < -logFC_t)#下调
k2 = (dep$P_value< P.Value_t)&(dep$fd > logFC_t)#上调
dep <- mutate(dep,change = ifelse(k1,"down",ifelse(k2,"up","stable")))

dep$protein=rownames(dep)
#改名字

library(clusterProfiler)
library(org.Hs.eg.db)
s2e <- bitr(dep$protein, 
            fromType = "UNIPROT",
            toType = "SYMBOL",
            OrgDb = org.Hs.eg.db)#人类
colnames(s2e)[1]="protein"
#有重复！！！
s2e=s2e[!duplicated(s2e$protein),]
dep_2=merge(dep,s2e)
dep_2$newname=paste0(dep_2$protein,"_",dep_2$SYMBOL)

dir.create("Result_cap_2")
dep_3=dep_2[dep_2$change!="stable",]
write.csv(dep_3,file = paste0("Result_cap_2/",theme2,"/",theme2,"_",theme3,"_ttest.csv"))

ogtt_overlap<-intersect(dep_3$protein,ogtt$protein)
length(ogtt_overlap)

Keto_overlap<-intersect(dep_3$protein,Keto$protein)
length(Keto_overlap)

HFD_overlap<-intersect(dep_3$protein,HFD$protein)
length(HFD_overlap)

Normal_overlap<-intersect(dep_3$protein,Normal$protein)
length(Normal_overlap)

LFD_overlap<-intersect(dep_3$protein,LFD$protein)
length(LFD_overlap)

library(plyr)
overlap<-rbind.fill(as.data.frame(ogtt_overlap),as.data.frame(Keto_overlap),as.data.frame(HFD_overlap),as.data.frame(Normal_overlap),as.data.frame(LFD_overlap))

dir.create(paste0("Result_cap_2/",theme2))
write.csv(overlap,file = paste0("Result_cap_2/",theme2,"/",theme2,"_",theme3,"_ttest_overlap.csv"))
