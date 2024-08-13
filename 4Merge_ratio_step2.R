rm(list = ls())

####1:删除36，NA80-------

#此处选包含batch36数据
load("diet_FDRnotsel.Rdata")


#去除含batch36的列，保存
options(max.print=1000000)
diet_data<-diet_data[,-grep("batch36",colnames(diet_data))]

#去除缺失值方法1：可筛选NA率：80%？
#计算行NA率
#del<-function(x){sum(is.na(x))/ncol(diet_data)*100} 
#apply(diet_data,1,del)
#diet_data<-diet_data[apply(diet_data,1,del)<80,]

save(diet_data,file="diet_FDRnotsel_without36.Rdata")



###2:删除36，在NAguide中筛选NA80%，sek填补（group=饮食+时间点），然后和test3比较结果是否一致------
#smple list整理
load("C:/Users/Cycy/Desktop/diet_R/sample_batch_list_remove36.Rdata")
clinical$diettime=paste0(clinical$diet,clinical$time)
sample_naguide_2=clinical[,c(1,9)]
write.csv(sample_naguide_2,file="sample_naguide_2.csv")

#这是group=timediet的结果2,log,normolnize=T------
diet_data2<-read.csv("diet_NAguide_data/diet_naguide_result_2.csv",row.names=1)
new_diet_data2 <- diet_data2[order(row.names(diet_data2)), order(colnames(diet_data2))]


#这是group=timediet的结果3,log,normolnize=F------
diet_data3<-read.csv("diet_NAguide_data/diet_naguide_result_3.csv",row.names=1)
new_diet_data3 <- diet_data3[order(row.names(diet_data3)), order(colnames(diet_data3))]


#这是group=batch的结果------
load("diet_FDRnotsel_remove36_NAguide_NA80_sek.Rdata")
new_diet_data <- diet_data[order(row.names(diet_data)), order(colnames(diet_data))]

#比较#############
all(new_diet_data==new_diet_data3)

library(daff)
d=diff_data(new_diet_data,new_diet_data3)

