rm(list = ls())

library("readxl")
library("tidyverse")
library("stringr")
library("vioplot")
library("openxlsx") 
library("dplyr")


#需筛选Protein FDR Confidence: Combined


#1先处理单个文件的列名和Low筛选！！后来决定不删Low列，直接跳到2
load("batch_list.Rdata")

result=list()
for (i in 1:length(list)) {
  list[[i]]%>%
    #筛选列
    select(contains("Accession")|contains("Protein FDR Confidence: Combined")|contains("Abundance Ratio: "))%>%
    as.data.frame()->batch
  rownames(batch)=batch[,1]
  batch<-batch[,-1]
  #不删了：筛选掉Low列：可能会导致蛋白有些batch有，有些没有，很多NA
  FDRsel<-which(batch$`Protein FDR Confidence: Combined`!="Low")
  batch<-batch[FDRsel,]%>%
    select(contains("Accession")|contains("Abundance Ratio: "))->batch
  #提取文件名中的batch名
    type_i<-str_replace(list.files("diet_protein_alldata/",pattern = ".xlsx")[i],".xlsx","")
  #修改列名
    colnames(batch)<-c(paste0(type_i,str_replace(colnames(batch),"Abundance Ratio: ","AR"),sep=""))
    #行名作为第一列，有助于下一步合并
    Accession=rownames(batch)
    batch=cbind(Accession,batch)
    #输出结果
    result_name=paste0(type_i,"FDRsel.xlsx",sep="")
    result[[i]]<-write.xlsx(batch,result_name,rowNames=F,colNames=T) }

#重新合并
list.files("diet_protein_alldata_FDRsel/",pattern = ".xlsx",full.names = TRUE)%>%
  lapply(readxl::read_excel)%>%
  reduce(full_join,by="Accession") -> diet_data

options(max.print=1000000)
diet_data<-as.data.frame(diet_data)
rownames(diet_data)=diet_data[,1]
diet_data<-diet_data[,-1]

save(diet_data,file="diet_data.Rdata")
write.xlsx(diet_data,"diet_data.xlsx",rowNames=T,colNames=T)


#2不删除Low列
load("batch_list.Rdata")

result=list()
for (i in 1:length(list)) {
  list[[i]]%>%
    #筛选列
    select(contains("Accession")|contains("Abundance Ratio: "))%>%
    as.data.frame()->batch
  rownames(batch)=batch[,1]
  batch<-batch[,-1]
  #提取文件名中的batch名
  type_i<-str_replace(list.files("diet_protein_alldata/",pattern = ".xlsx")[i],".xlsx","")
  #修改列名
  colnames(batch)<-c(paste0(type_i,str_replace(colnames(batch),"Abundance Ratio: ","AR"),sep=""))
  #行名作为第一列，有助于下一步合并
  Accession=rownames(batch)
  batch=cbind(Accession,batch)
  #输出结果
  result_name=paste0(type_i,"FDRnotsel.xlsx",sep="")
  result[[i]]<-write.xlsx(batch,result_name,rowNames=F,colNames=T) }

#重新合并
list.files("diet_protein_alldata_FDRnotsel/",pattern = ".xlsx",full.names = TRUE)%>%
  lapply(readxl::read_excel)%>%
  reduce(full_join,by="Accession") -> diet_data

options(max.print=1000000)
diet_data<-as.data.frame(diet_data)
rownames(diet_data)=diet_data[,1]
diet_data<-diet_data[,-1]

#列名处理
colnames(diet_data)<-gsub('[(]', '.', colnames(diet_data))
colnames(diet_data)<-gsub('[)/]', '', colnames(diet_data))
colnames(diet_data)<-gsub(' ', '', colnames(diet_data))

#去除异常行"Q5T0Z8"
options(max.print=1000000)
diet_data<-diet_data[-grep("Q5T0Z8",rownames(diet_data)),]


save(diet_data,file="diet_FDRnotsel.Rdata")






  





