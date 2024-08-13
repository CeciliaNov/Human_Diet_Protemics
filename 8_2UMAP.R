rm(list=ls())

library(umap)
library(ggsci)
library(ggplot2)
#设定颜色差异较大的颜色
colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))



###test5:删除36，NA80,NAguide sek插补----
load("sample_batch_list_remove36.Rdata")
clinical$timediet=paste0(clinical$diet,clinical$time)
Participant=as.character(clinical$bianhao)
Time=clinical$time
Diet=clinical$diet
Diet_time=clinical$timediet
Gender=clinical$gender
Weight=clinical$weight

##after
load("combat_FDRnotsel_remove36_NA80_sek.Rdata")
myd=as.data.frame(t(combat))

#共有
myumap<-umap(myd)
mydf<-data.frame(myumap$layout)

p<-ggplot(mydf,aes(x=X1, y=X2, colour=Diet)) +
  scale_color_manual(values=colpalettes)+
  geom_point(size=1)+
  #geom_jitter(width = 0.5,height = 0.5)+
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5),
        panel.background = element_blank())
p

#after保存
ggsave(p,filename = "umap_NA80_sek_after.pdf",width = 10,height = 10)  
