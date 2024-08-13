rm(list = ls())#清空环境变量

library("readxl")
library("tidyverse")
library("stringr")
library("vioplot")


list.files("diet_protein_alldata/",pattern = ".xlsx",full.names = TRUE)%>%
  lapply(readxl::read_excel)->list

save(list,file="batch_list.Rdata")

load("batch_list.Rdata")

TR_list<-read.csv("TR_list.csv")


for (i in 1:length(list)) {
  list[[i]]%>%
  select(contains("Accession")|contains("Abundances (Grouped): 133C")|contains(TR_list[which(TR_list$batch==str_replace(list.files("diet_protein_alldata/",pattern = ".xlsx")[i],".xlsx","")),2]))%>%
  as.data.frame()->batch

  rownames(batch)=batch[,1]
  batch<-batch[,-1]
  batch<-log2(batch)
  
  
  head(batch)
  dim(batch)
  sdbatch <- apply(batch, 1, function(sd){sd(sd,na.rm = T)})
  mebatch <- apply(batch, 1, function(me){mean(me,na.rm = T)})
  cvbatch <- sdbatch/mebatch
  
  
  type_i<-str_replace(list.files("diet_protein_alldata/",pattern = ".xlsx")[i],".xlsx","")
  dfvp <- data.frame(type=paste(type_i,"technical repetition",sep ="  "),cv=cvbatch)
  
  result_name=paste0(type_i,"TR.pdf",sep="")
  pdf(result_name,width=5,height=5)
  plot=vioplot(cv~type,data = dfvp,
               main = "coefficient of variation of QC samples",
               col=c("#003C67FF"),
               ylim=c(0,1))
  print(plot)
  dev.off()}

