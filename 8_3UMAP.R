rm(list=ls())
load("sample_batch_list_remove36.Rdata")
load("combat_FDRnotsel_remove36_NA80_sek.Rdata")

library(umap)
library(ggsci)
library(ggplot2)
library(Rtsne)
#设定颜色差异较大的颜色
colpalettes<-unique(c(pal_npg("nrc")(10),pal_aaas("default")(10),pal_nejm("default")(8),pal_lancet("lanonc")(9),
                      pal_jama("default")(7),pal_jco("default")(10),pal_ucscgb("default")(26),pal_d3("category10")(10),
                      pal_locuszoom("default")(7),pal_igv("default")(51),
                      pal_uchicago("default")(9),pal_startrek("uniform")(7),
                      pal_tron("legacy")(7),pal_futurama("planetexpress")(12),pal_rickandmorty("schwifty")(12),
                      pal_simpsons("springfield")(16),pal_gsea("default")(12)))


dat=as.data.frame(t(combat))
dat$ID=row.names(dat)

dat2=merge(dat,clinical)
set.seed(42)
tsne_out = Rtsne(dat,perplexity = 30, check_duplicates = FALSE)

clinical$timediet=paste0(clinical$diet,clinical$time)
Participant=as.character(clinical$bianhao)
Time=clinical$time
Diet=clinical$diet
Diet_time=clinical$timediet
Gender=clinical$gender
Weight=clinical$weight

plot(tsne_out$Y,
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2",
     col=rainbow(33)[as.numeric(as.factor(Weight))])
