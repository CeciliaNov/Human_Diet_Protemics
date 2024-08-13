rm(list = ls())

library("openxlsx") 
library(sva)
library(bladderbatch)
library("stringr")
library("tidyverse")


###临床表格前处理（with36和without36版本）----
{clinical<-read.csv("sample_batch_list.csv",check.names=F,row.names=1)
save(clinical,file="sample_batch_list.Rdata")

clinical<-clinical[-grep("batch36",rownames(clinical)),]

save(clinical,file="sample_batch_list_remove36.Rdata")}


###test1:未删除36，NA80，用0替代----

load("diet_FDRnotsel_NA80_rep0.Rdata")
load("sample_batch_list.Rdata")


batch=clinical$batch
combat_diet_data <- ComBat(dat = diet_data, 
                           batch =batch)

combat<-as.data.frame(combat_diet_data)

save(combat,file="combat_FDRnotsel_NA80_rep0.Rdata")

#可视化1(PCA的偷懒画法)，最后决定删batch36
{ 
  #校正前PCA
  colpalettes<-unique(c(pal_igv("default")(51)))
  
  load("sample_batch_list.Rdata")
  Group=clinical$batch
  Group=factor(Group,levels=c(paste0("batch",c(1:36))))
  
  
  load("diet_FDRnotsel_NA80_rep0.Rdata")
  dat=diet_data
 
  pca.plot = draw_pca(dat,Group,color = colpalettes,
                      addEllipses = TRUE,
                      style = "default",
                      color.label = "Group")
  pca.plot
  dev.off()
  
  #校正后PCA
  colpalettes<-unique(c(pal_igv("default")(51)))
  
  load("sample_batch_list.Rdata")
  Group=clinical$batch
  Group=factor(Group,levels=c(paste0("batch",c(1:36))))
  
  load("combat_FDRnotsel_NA80_rep0.Rdata")
  dat=combat

  pca.plot = draw_pca(dat,Group,color = colpalettes,
                      addEllipses = TRUE,
                      style = "default",
                      color.label = "Group")
  pca.plot
  dev.off()}
#可视化2(PCA的基础画法：只画PCA_after)
{load("diet_data_afterNA.Rdata")
  load("combat_diet_data.Rdata")
  sample_batch_list<-read.csv("sample_batch_list.csv",check.names=F,row.names=1)
  Group=sample_batch_list$batch
  Group=factor(Group,levels=c(paste0("batch",c(1:36))))
  
  exp=combat_diet_data
  dat=as.data.frame(t(exp))
  library(FactoMineR)
  library(factoextra) 
  dat.pca <- PCA(dat, graph = FALSE)
  pca_plot <- fviz_pca_ind(dat.pca,
                           geom.ind = c("point", "text") , # show points only (nbut not "text")
                           col.ind = Group, # color by groups
                           palette = colpalettes,
                           addEllipses = TRUE, # Concentration ellipses
                           legend.title = "Groups",
                           repel = TRUE)
  pca_plot}








###test2:删除36，NA0----
load("diet_FDRnotsel_remove36_NA0.Rdata")
load("sample_batch_list_remove36.Rdata")

batch=clinical$batch
combat_diet_data <- ComBat(dat = diet_data, 
                           batch =batch,
                           par.prior = FALSE)

combat<-as.data.frame(combat_diet_data)

save(combat,file="combat_FDRnotsel_remove36_NA0.Rdata")




###test3:删除36，NA80,NAguide impseqrob插补----
load("diet_FDRnotsel_remove36_NAguide_NA80_imp.Rdata")
load("sample_batch_list_remove36.Rdata")

batch=clinical$batch
combat_diet_data <- ComBat(dat = diet_data, 
                           batch =batch,
                           par.prior = FALSE)

combat<-as.data.frame(combat_diet_data)

save(combat,file="combat_FDRnotsel_remove36_NA80_imp.Rdata")



###test4:删除36，NA70,NAguide sek插补----
load("diet_FDRnotsel_remove36_NAguide_NA70_sek.Rdata")
load("sample_batch_list_remove36.Rdata")

batch=clinical$batch
combat_diet_data <- ComBat(dat = diet_data, 
                           batch =batch,
                           par.prior = FALSE)

combat<-as.data.frame(combat_diet_data)

save(combat,file="combat_FDRnotsel_remove36_NA70_sek.Rdata")



###test5:删除36，NA80,NAguide sek插补----
load("diet_FDRnotsel_remove36_NAguide_NA80_sek.Rdata")
load("sample_batch_list_remove36.Rdata")

batch=clinical$batch
combat_diet_data <- ComBat(dat = diet_data, 
                           batch =batch,
                           par.prior = FALSE)

combat<-as.data.frame(combat_diet_data)

save(combat,file="combat_FDRnotsel_remove36_NA80_sek.Rdata")

