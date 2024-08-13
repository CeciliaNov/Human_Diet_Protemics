
library("readxl")
library("tidyverse")
library("stringr")
library("vioplot")

read_excel("diet_protein_alldata/batch17.xlsx")%>%
select(contains("Accession")|contains("Abundances (Grouped): 133C")|contains("Abundances (Grouped): 128C"))%>%
as.data.frame()->batch
  
  rownames(batch)=batch[,1]
  batch<-batch[,-1]
  batch<-log2(batch)
  #可选
  na.omit(batch)
  pool<-pool[-grep("Q5T0Z8",rownames(pool)),]
  
  head(batch)
  dim(batch)
  sdbatch <- apply(batch, 1, function(sd){sd(sd,na.rm = T)})
  mebatch <- apply(batch, 1, function(me){mean(me,na.rm = T)})
  cvbatch <- sdbatch/mebatch
  
  dfvp <- data.frame(type="batch17 technical repetition_2",cv=cvbatch)
  
  result_name=paste0("batch17TR_2.pdf",sep="")
  pdf(result_name,width=5,height=5)
  plot=vioplot(cv~type,data = dfvp,
               main = "coefficient of variation of QC samples",
               col=c("#003C67FF"),
               ylim=c(0,1))
  dev.off()

