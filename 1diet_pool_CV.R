library("readxl")
library("tidyverse")
library("stringr")
library("vioplot")


list.files("diet_protein_alldata/",pattern = ".xlsx",full.names = TRUE)%>%
  lapply(readxl::read_excel)%>%
  reduce(full_join,by="Accession")%>%
  arrange(Accession)%>%
  select(contains("Accession")|contains("Abundances (Grouped): 126"))->df

list.files("diet_protein_alldata/",pattern = ".xlsx",full.names = FALSE)
colnames(df)<-c("Accession",str_replace(list.files("diet_protein_alldata/",pattern = ".xlsx"),".xlsx",""))

save(df,file = "pool_merge.Rdata")

load("pool_merge.Rdata")

#没有这一步会报错，可能是"tidyverse"的问题？
pool<-as.data.frame(df)

rownames(pool)=pool[,1]
pool<-pool[,-1]

#可省略：去除异常行"Q5T0Z8"
options(max.print=1000000)
pool<-pool[-grep("Q5T0Z8",rownames(pool)),]

pool<-log2(pool)
sdpool <- apply(pool, 1, function(sd){sd(sd,na.rm = T)})
mepool <- apply(pool, 1, function(me){mean(me,na.rm = T)})
cvpool <- sdpool/mepool

dfvp <- data.frame(type="pool", cv=cvpool)

vioplot(cv~type,data = dfvp,
        main = "coefficient of variation of QC samples",
        col=c("#003C67FF"),
        ylim=c(0,1))

dev.off()

