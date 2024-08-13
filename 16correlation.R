library(readxl)
# library(corrplot)
# library(ComplexHeatmap)
# library(circlize)
library(pheatmap)
# library(ggplot2)
library(foreign)
setwd("D:/Cycy1/R_Proj/Diet/cor_R_plasma and urine")

# Plasma_A
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^A',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
physiological = read_xls('检验data_A.xls')
physiological$Serialnumber = sprintf('%02d',physiological$Serialnumber)
physiological$Time = sprintf('%02d',physiological$Time)
physiological$BH = paste0('ABH',physiological$Serialnumber,'SJ',physiological$Time)
physiological = physiological[match(colnames(df),physiological$BH),]
physiological = t(physiological[,c('GLU','insulin','TG','FFA')])
colnames(physiological) = colnames(df)

df = t(df)
physiological = t(physiological)

rho_df = data.frame()
p_df = data.frame()

for (i in colnames(df)) {
  for (j in colnames(physiological)) {
    result = cor.test(df[,i],physiological[,j],method='spearman')
    rho_df[i,j] = result$estimate["rho"]
    p_df[i,j] = result$p.value
  }
}
rho_df = as.matrix(rho_df)
p_df = as.matrix(p_df)
# p_df = ifelse(p_df<0.05,round(p_df,3),'')
for (i in 1:length(p_df)) {
  if (as.numeric(p_df[i]) < 0.001) {
    p_df[i] = '***'
  } else if (as.numeric(p_df[i]) < 0.01) {
    p_df[i] = '**'
  } else if (as.numeric(p_df[i]) < 0.05) {
    p_df[i] = '*'
  } else {
    p_df[i] = ''
  }
}
# n = 1
# df_heatmap = data.frame()
# for (i in colnames(p_df)) {
#   for (j in rownames(p_df)) {
#     df_heatmap[n,'x'] = i
#     df_heatmap[n,'y'] = j
#     df_heatmap[n,'spearman'] = rho_df[j,i]
#     df_heatmap[n,'p'] = p_df[j,i]
#     n = n+1
#   }
# }

pdf("plasma_A.pdf")
# col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# Heatmap(rho_df, cluster_rows=F, cluster_columns=F,
#         rect_gp = gpar(col = "black"),
#         col=col_fun,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(p_df[i,j], x, y, vjust='centre', gp = gpar(fontsize = 12,col="black"))},
#         column_names_rot=45, row_names_side="left",
#         name='Spearman R', border='black'
# )

pheatmap(rho_df, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma A', display_numbers=p_df, fontsize_number=10, fontsize=12
)

# ggplot(df_heatmap, aes(x,y)) +
#   geom_tile(aes(fill=spearman),color='black') +
#   geom_text(aes(label=p), color="black", size=4, vjust=0.7) + 
#   scale_fill_gradient2(low='blue', high='red',mid = 'white', limit=c(-1,1), name="Correlation") + 
#   labs(x=NULL,y=NULL) + # 去掉横纵坐标标题
#   theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
#         axis.text.y = element_text(size=8,color = "black"),
#         axis.ticks.x = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.background=element_blank()) # 做一些简单的美化

dev.off()

# Plasma_B
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^B',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
physiological = read.spss('clinical_data_B.sav',to.data.frame=T)
physiological$SerialNumber = sprintf('%02d',as.numeric(physiological$SerialNumber))
physiological$Time = sprintf('%02d',physiological$Time)
physiological$BH = paste0('BBH',physiological$SerialNumber,'SJ',physiological$Time)
physiological = physiological[match(colnames(df),physiological$BH),]
physiological = t(physiological[,c('GLU','insulin','TG','FFA')])
colnames(physiological) = colnames(df)

df = t(df)
physiological = t(physiological)

rho_df = data.frame()
p_df = data.frame()

for (i in colnames(df)) {
  for (j in colnames(physiological)) {
    result = cor.test(df[,i],physiological[,j],method='spearman')
    rho_df[i,j] = result$estimate["rho"]
    p_df[i,j] = result$p.value
  }
}
rho_df = as.matrix(rho_df)
p_df = as.matrix(p_df)
for (i in 1:length(p_df)) {
  if (as.numeric(p_df[i]) < 0.001) {
    p_df[i] = '***'
  } else if (as.numeric(p_df[i]) < 0.01) {
    p_df[i] = '**'
  } else if (as.numeric(p_df[i]) < 0.05) {
    p_df[i] = '*'
  } else {
    p_df[i] = ''
  }
}

pdf("plasma_B.pdf")
pheatmap(rho_df, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma B', display_numbers=p_df, fontsize_number=10, fontsize=12
)
dev.off()

# Plasma_C
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^C',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
physiological = read.spss('clinical_data_C.sav',to.data.frame=T)
physiological$SerialNumber = sprintf('%02d',as.numeric(physiological$SerialNumber))
physiological$Time = sprintf('%02d',physiological$Time)
physiological$BH = paste0('CBH',physiological$SerialNumber,'SJ',physiological$Time)
physiological = physiological[match(colnames(df),physiological$BH),]
physiological = t(physiological[,c('GLU','insulin','TG','FFA')])
colnames(physiological) = colnames(df)

df = t(df)
physiological = t(physiological)

rho_df = data.frame()
p_df = data.frame()

for (i in colnames(df)) {
  for (j in colnames(physiological)) {
    result = cor.test(df[,i],physiological[,j],method='spearman')
    rho_df[i,j] = result$estimate["rho"]
    p_df[i,j] = result$p.value
  }
}
rho_df = as.matrix(rho_df)
p_df = as.matrix(p_df)
for (i in 1:length(p_df)) {
  if (as.numeric(p_df[i]) < 0.001) {
    p_df[i] = '***'
  } else if (as.numeric(p_df[i]) < 0.01) {
    p_df[i] = '**'
  } else if (as.numeric(p_df[i]) < 0.05) {
    p_df[i] = '*'
  } else {
    p_df[i] = ''
  }
}

pdf("plasma_C.pdf")
pheatmap(rho_df, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma C', display_numbers=p_df, fontsize_number=10, fontsize=12
)
dev.off()

# Plasma_D
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^D',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
physiological = read.spss('clinical_data_D.sav',to.data.frame=T)
physiological$serialnumber = sprintf('%02d',as.numeric(physiological$serialnumber))
physiological$Time = sprintf('%02d',physiological$Time)
physiological$BH = paste0('DBH',physiological$serialnumber,'SJ',physiological$Time)
physiological = physiological[match(colnames(df),physiological$BH),]
physiological = t(physiological[,c('GLU','insulin','TG','FFA')])
colnames(physiological) = colnames(df)

df = t(df)
physiological = t(physiological)

rho_df = data.frame()
p_df = data.frame()

for (i in colnames(df)) {
  for (j in colnames(physiological)) {
    result = cor.test(df[,i],physiological[,j],method='spearman')
    rho_df[i,j] = result$estimate["rho"]
    p_df[i,j] = result$p.value
  }
}
rho_df = as.matrix(rho_df)
p_df = as.matrix(p_df)
for (i in 1:length(p_df)) {
  if (as.numeric(p_df[i]) < 0.001) {
    p_df[i] = '***'
  } else if (as.numeric(p_df[i]) < 0.01) {
    p_df[i] = '**'
  } else if (as.numeric(p_df[i]) < 0.05) {
    p_df[i] = '*'
  } else {
    p_df[i] = ''
  }
}

pdf("plasma_D.pdf")
pheatmap(rho_df, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma D', display_numbers=p_df, fontsize_number=10, fontsize=12
)
dev.off()

# Plasma_E
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^E',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
physiological = read.spss('clinical_data_E.sav',to.data.frame=T)
physiological$Serialnumber = sprintf('%02d',as.numeric(physiological$Serialnumber))
physiological$Time = sprintf('%02d',physiological$Time)
physiological$BH = paste0('EBH',physiological$Serialnumber,'SJ',physiological$Time)
physiological = physiological[match(colnames(df),physiological$BH),]
physiological = t(physiological[,c('GLU','insulin','TG','FFA')])
colnames(physiological) = colnames(df)

df = t(df)
physiological = t(physiological)

rho_df = data.frame()
p_df = data.frame()

for (i in colnames(df)) {
  for (j in colnames(physiological)) {
    result = cor.test(df[,i],physiological[,j],method='spearman')
    rho_df[i,j] = result$estimate["rho"]
    p_df[i,j] = result$p.value
  }
}
rho_df = as.matrix(rho_df)
p_df = as.matrix(p_df)
for (i in 1:length(p_df)) {
  if (as.numeric(p_df[i]) < 0.001) {
    p_df[i] = '***'
    #print(paste0(i,"  ",p_df[i]))
  } else if (as.numeric(p_df[i]) < 0.01) {
    p_df[i] = '**'
    #print(paste0(i,"  ",p_df[i]))
  } else if (as.numeric(p_df[i]) < 0.05) {
    p_df[i] = '*'
    #print(paste0(i,"  ",p_df[i]))
  } else {
    p_df[i] = ''
    #print(paste0(i,"  ",p_df[i]))
  }
}

pdf("plasma_E.pdf")
pheatmap(rho_df, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma E', display_numbers=p_df, fontsize_number=10, fontsize=12
)
dev.off()




### Plasma Overweight vs Normal ###
rm(list=ls())
# load('combat_serum_FDRnotsel_NA20_seqknn.Rdata')
# load('diet_protein_alldata/data/sample2.Rdata')
# colnames(combat) = sample2$sample[match(colnames(combat), sample2$ID)]
# save(combat,file='combat_serum_FDRnotsel_NA20_seqknn_sample.Rdata')

# Plasma_A
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^A',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
phy = read_xls('检验data_A.xls')
phy$Serialnumber = sprintf('%02d',phy$Serialnumber)
phy$Time = sprintf('%02d',phy$Time)
phy$BH = paste0('ABH',phy$Serialnumber,'SJ',phy$Time)
samples = intersect(colnames(df),phy$BH)
df = df[,samples]
phy = phy[phy$BH %in% samples,]
phy = phy[match(colnames(df),phy$BH),]

phy_fat = phy[phy$Overweight==1,]
phy_nor = phy[phy$Overweight==0,]

phy_fat = t(phy_fat[,c('GLU','insulin','TG','FFA','BH')])
phy_nor = t(phy_nor[,c('GLU','insulin','TG','FFA','BH')])

colnames(phy_fat) = phy_fat['BH',]
colnames(phy_nor) = phy_nor['BH',]
df_fat = df[,colnames(df) %in% colnames(phy_fat)]
df_nor = df[,colnames(df) %in% colnames(phy_nor)]

phy_fat = t(phy_fat[-nrow(phy_fat),])
phy_nor = t(phy_nor[-nrow(phy_nor),])
df_fat = t(df_fat)
df_nor = t(df_nor)

rho_df_fat = data.frame()
rho_df_nor = data.frame()
p_df_fat = data.frame()
p_df_nor = data.frame()

for (i in colnames(df_fat)) {
  for (j in colnames(phy_fat)) {
    result = cor.test(as.numeric(df_fat[,i]),as.numeric(phy_fat[,j]),method='spearman')
    rho_df_fat[i,j] = result$estimate["rho"]
    p_df_fat[i,j] = result$p.value
  }
}
rho_df_fat = as.matrix(rho_df_fat)
p_df_fat = as.matrix(p_df_fat)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_fat[i]) < 0.001) {
    p_df_fat[i] = '***'
  } else if (as.numeric(p_df_fat[i]) < 0.01) {
    p_df_fat[i] = '**'
  } else if (as.numeric(p_df_fat[i]) < 0.05) {
    p_df_fat[i] = '*'
  } else {
    p_df_fat[i] = ''
  }
}

for (i in colnames(df_nor)) {
  for (j in colnames(phy_nor)) {
    result = cor.test(as.numeric(df_nor[,i]),as.numeric(phy_nor[,j]),method='spearman')
    rho_df_nor[i,j] = result$estimate["rho"]
    p_df_nor[i,j] = result$p.value
  }
}
rho_df_nor = as.matrix(rho_df_nor)
p_df_nor = as.matrix(p_df_nor)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_nor[i]) < 0.001) {
    p_df_nor[i] = '***'
  } else if (as.numeric(p_df_nor[i]) < 0.01) {
    p_df_nor[i] = '**'
  } else if (as.numeric(p_df_nor[i]) < 0.05) {
    p_df_nor[i] = '*'
  } else {
    p_df_nor[i] = ''
  }
}


# n = 1
# df_heatmap = data.frame()
# for (i in colnames(p_df)) {
#   for (j in rownames(p_df)) {
#     df_heatmap[n,'x'] = i
#     df_heatmap[n,'y'] = j
#     df_heatmap[n,'spearman'] = rho_df[j,i]
#     df_heatmap[n,'p'] = p_df[j,i]
#     n = n+1
#   }
# }

pdf("plasma_A_fat.pdf")
# col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# Heatmap(rho_df, cluster_rows=F, cluster_columns=F,
#         rect_gp = gpar(col = "black"),
#         col=col_fun,
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(p_df[i,j], x, y, vjust='centre', gp = gpar(fontsize = 12,col="black"))},
#         column_names_rot=45, row_names_side="left",
#         name='Spearman R', border='black'
# )

pheatmap(rho_df_fat, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma A Fat', display_numbers=p_df_fat, fontsize_number=10, fontsize=12
)
dev.off()

pdf("plasma_A_nor.pdf")
pheatmap(rho_df_nor, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, show_rownames = FALSE,
         cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma A Normal', display_numbers=p_df_nor, fontsize_number=10, fontsize=12
)
dev.off()

# ggplot(df_heatmap, aes(x,y)) +
#   geom_tile(aes(fill=spearman),color='black') +
#   geom_text(aes(label=p), color="black", size=4, vjust=0.7) + 
#   scale_fill_gradient2(low='blue', high='red',mid = 'white', limit=c(-1,1), name="Correlation") + 
#   labs(x=NULL,y=NULL) + # 去掉横纵坐标标题
#   theme(axis.text.x = element_text(size=8,angle = 30,hjust = 1,color = "black"),
#         axis.text.y = element_text(size=8,color = "black"),
#         axis.ticks.x = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.background=element_blank()) # 做一些简单的美化

# dev.off()

# Plasma_B
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^B',colnames(df))]
df = df[,-grep('^BBH13SJ',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
phy = read.spss('clinical_data_B.sav',to.data.frame=T)
phy$SerialNumber = sprintf('%02d',as.numeric(phy$SerialNumber))
phy$Time = sprintf('%02d',phy$Time)
phy$BH = paste0('BBH',phy$SerialNumber,'SJ',phy$Time)
samples = intersect(colnames(df),phy$BH)
df = df[,samples]
phy = phy[phy$BH %in% samples,]
phy = phy[match(colnames(df),phy$BH),]

phy_fat = phy[phy$Overweight==1,]
phy_nor = phy[phy$Overweight==0,]

phy_fat = t(phy_fat[,c('GLU','insulin','TG','FFA','BH')])
phy_nor = t(phy_nor[,c('GLU','insulin','TG','FFA','BH')])

colnames(phy_fat) = phy_fat['BH',]
colnames(phy_nor) = phy_nor['BH',]
df_fat = df[,colnames(df) %in% colnames(phy_fat)]
df_nor = df[,colnames(df) %in% colnames(phy_nor)]

phy_fat = t(phy_fat[-nrow(phy_fat),])
phy_nor = t(phy_nor[-nrow(phy_nor),])
df_fat = t(df_fat)
df_nor = t(df_nor)

rho_df_fat = data.frame()
rho_df_nor = data.frame()
p_df_fat = data.frame()
p_df_nor = data.frame()

for (i in colnames(df_fat)) {
  for (j in colnames(phy_fat)) {
    result = cor.test(as.numeric(df_fat[,i]),as.numeric(phy_fat[,j]),method='spearman')
    rho_df_fat[i,j] = result$estimate["rho"]
    p_df_fat[i,j] = result$p.value
  }
}
rho_df_fat = as.matrix(rho_df_fat)
p_df_fat = as.matrix(p_df_fat)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_fat[i]) < 0.001) {
    p_df_fat[i] = '***'
  } else if (as.numeric(p_df_fat[i]) < 0.01) {
    p_df_fat[i] = '**'
  } else if (as.numeric(p_df_fat[i]) < 0.05) {
    p_df_fat[i] = '*'
  } else {
    p_df_fat[i] = ''
  }
}

for (i in colnames(df_nor)) {
  for (j in colnames(phy_nor)) {
    result = cor.test(as.numeric(df_nor[,i]),as.numeric(phy_nor[,j]),method='spearman')
    rho_df_nor[i,j] = result$estimate["rho"]
    p_df_nor[i,j] = result$p.value
  }
}
rho_df_nor = as.matrix(rho_df_nor)
p_df_nor = as.matrix(p_df_nor)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_nor[i]) < 0.001) {
    p_df_nor[i] = '***'
  } else if (as.numeric(p_df_nor[i]) < 0.01) {
    p_df_nor[i] = '**'
  } else if (as.numeric(p_df_nor[i]) < 0.05) {
    p_df_nor[i] = '*'
  } else {
    p_df_nor[i] = ''
  }
}

pdf("plasma_B_fat.pdf")
pheatmap(rho_df_fat, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma B Fat', display_numbers=p_df_fat, fontsize_number=10, fontsize=12
)
dev.off()
pdf("plasma_B_nor.pdf")
pheatmap(rho_df_nor, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none', show_rownames = FALSE,
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma B Normal', display_numbers=p_df_nor, fontsize_number=10, fontsize=12
)
dev.off()


# Plasma_C
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')
trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^C',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
phy = read.spss('clinical_data_C.sav',to.data.frame=T)
phy$SerialNumber = sprintf('%02d',as.numeric(phy$SerialNumber))
phy$Time = sprintf('%02d',phy$Time)
phy$BH = paste0('CBH',phy$SerialNumber,'SJ',phy$Time)
samples = intersect(colnames(df),phy$BH)
df = df[,samples]
phy = phy[phy$BH %in% samples,]
phy = phy[match(colnames(df),phy$BH),]

phy_fat = phy[phy$Overweight==1,]
phy_nor = phy[phy$Overweight==0,]

phy_fat = t(phy_fat[,c('GLU','insulin','TG','FFA','BH')])
phy_nor = t(phy_nor[,c('GLU','insulin','TG','FFA','BH')])

colnames(phy_fat) = phy_fat['BH',]
colnames(phy_nor) = phy_nor['BH',]
df_fat = df[,colnames(df) %in% colnames(phy_fat)]
df_nor = df[,colnames(df) %in% colnames(phy_nor)]

phy_fat = t(phy_fat[-nrow(phy_fat),])
phy_nor = t(phy_nor[-nrow(phy_nor),])
df_fat = t(df_fat)
df_nor = t(df_nor)

rho_df_fat = data.frame()
rho_df_nor = data.frame()
p_df_fat = data.frame()
p_df_nor = data.frame()

for (i in colnames(df_fat)) {
  for (j in colnames(phy_fat)) {
    result = cor.test(as.numeric(df_fat[,i]),as.numeric(phy_fat[,j]),method='spearman')
    rho_df_fat[i,j] = result$estimate["rho"]
    p_df_fat[i,j] = result$p.value
  }
}
rho_df_fat = as.matrix(rho_df_fat)
p_df_fat = as.matrix(p_df_fat)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_fat[i]) < 0.001) {
    p_df_fat[i] = '***'
  } else if (as.numeric(p_df_fat[i]) < 0.01) {
    p_df_fat[i] = '**'
  } else if (as.numeric(p_df_fat[i]) < 0.05) {
    p_df_fat[i] = '*'
  } else {
    p_df_fat[i] = ''
  }
}

for (i in colnames(df_nor)) {
  for (j in colnames(phy_nor)) {
    result = cor.test(as.numeric(df_nor[,i]),as.numeric(phy_nor[,j]),method='spearman')
    rho_df_nor[i,j] = result$estimate["rho"]
    p_df_nor[i,j] = result$p.value
  }
}
rho_df_nor = as.matrix(rho_df_nor)
p_df_nor = as.matrix(p_df_nor)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_nor[i]) < 0.001) {
    p_df_nor[i] = '***'
  } else if (as.numeric(p_df_nor[i]) < 0.01) {
    p_df_nor[i] = '**'
  } else if (as.numeric(p_df_nor[i]) < 0.05) {
    p_df_nor[i] = '*'
  } else {
    p_df_nor[i] = ''
  }
}

pdf("plasma_C_fat.pdf")
pheatmap(rho_df_fat, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma C Fat', display_numbers=p_df_fat, fontsize_number=10, fontsize=12
)
dev.off()
pdf("plasma_C_nor.pdf")
pheatmap(rho_df_nor, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none', show_rownames = FALSE,
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma C Normal', display_numbers=p_df_nor, fontsize_number=10, fontsize=12
)
dev.off()


# Plasma_D
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')

trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^D',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
phy = read.spss('clinical_data_D.sav',to.data.frame=T)
phy$serialnumber = sprintf('%02d',as.numeric(phy$serialnumber))
phy$Time = sprintf('%02d',phy$Time)
phy$BH = paste0('DBH',phy$serialnumber,'SJ',phy$Time)
samples = intersect(colnames(df),phy$BH)
df = df[,samples]
phy = phy[phy$BH %in% samples,]
phy = phy[match(colnames(df),phy$BH),]

phy_fat = phy[phy$Overweight==1,]
phy_nor = phy[phy$Overweight==0,]

phy_fat = t(phy_fat[,c('GLU','insulin','TG','FFA','BH')])
phy_nor = t(phy_nor[,c('GLU','insulin','TG','FFA','BH')])

colnames(phy_fat) = phy_fat['BH',]
colnames(phy_nor) = phy_nor['BH',]
df_fat = df[,colnames(df) %in% colnames(phy_fat)]
df_nor = df[,colnames(df) %in% colnames(phy_nor)]

phy_fat = t(phy_fat[-nrow(phy_fat),])
phy_nor = t(phy_nor[-nrow(phy_nor),])
df_fat = t(df_fat)
df_nor = t(df_nor)

rho_df_fat = data.frame()
rho_df_nor = data.frame()
p_df_fat = data.frame()
p_df_nor = data.frame()

for (i in colnames(df_fat)) {
  for (j in colnames(phy_fat)) {
    result = cor.test(as.numeric(df_fat[,i]),as.numeric(phy_fat[,j]),method='spearman')
    rho_df_fat[i,j] = result$estimate["rho"]
    p_df_fat[i,j] = result$p.value
  }
}
rho_df_fat = as.matrix(rho_df_fat)
p_df_fat = as.matrix(p_df_fat)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_fat[i]) < 0.001) {
    p_df_fat[i] = '***'
  } else if (as.numeric(p_df_fat[i]) < 0.01) {
    p_df_fat[i] = '**'
  } else if (as.numeric(p_df_fat[i]) < 0.05) {
    p_df_fat[i] = '*'
  } else {
    p_df_fat[i] = ''
  }
}

for (i in colnames(df_nor)) {
  for (j in colnames(phy_nor)) {
    result = cor.test(as.numeric(df_nor[,i]),as.numeric(phy_nor[,j]),method='spearman')
    rho_df_nor[i,j] = result$estimate["rho"]
    p_df_nor[i,j] = result$p.value
  }
}
rho_df_nor = as.matrix(rho_df_nor)
p_df_nor = as.matrix(p_df_nor)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_nor[i]) < 0.001) {
    p_df_nor[i] = '***'
  } else if (as.numeric(p_df_nor[i]) < 0.01) {
    p_df_nor[i] = '**'
  } else if (as.numeric(p_df_nor[i]) < 0.05) {
    p_df_nor[i] = '*'
  } else {
    p_df_nor[i] = ''
  }
}

pdf("plasma_D_fat.pdf")
pheatmap(rho_df_fat, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma D Fat', display_numbers=p_df_fat, fontsize_number=10, fontsize=12
)
dev.off()
pdf("plasma_D_nor.pdf")
pheatmap(rho_df_nor, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none', show_rownames = FALSE,
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma D Normal', display_numbers=p_df_nor, fontsize_number=10, fontsize=12
)
dev.off()


# Plasma_E
rm(list=ls())
load("combat_FDRnotsel_remove36_NA80_sek_sample.Rdata")
trans = read_xlsx("idmapping_2023_11_09.xlsx")
genes = c('FTH1','PF4V1','RAP1B','LMNA',
          'ARPC2','FTL','CA2','PPBP',
          'SRGN','THBS1','APOA5','SOD1',
          'TXN','IDH1','CAT','PRDX1',
          'PRDX2','S100A7','ORM1','MIF',
          'S100A4','CRP','CA1','FLNA')

trans2 = trans[trans$Gene %in% genes,]
df = combat2[rownames(combat2) %in% trans2$From,]
df = df[,grep('^E',colnames(df))]
rownames(df) = trans$Gene[match(rownames(df),trans$From)]
df = df[order(rownames(df)),]
phy = read.spss('clinical_data_E.sav',to.data.frame=T)
phy$Serialnumber = sprintf('%02d',as.numeric(phy$Serialnumber))
phy$Time = sprintf('%02d',phy$Time)
phy$BH = paste0('EBH',phy$Serialnumber,'SJ',phy$Time)
samples = intersect(colnames(df),phy$BH)
df = df[,samples]
phy = phy[phy$BH %in% samples,]
phy = phy[match(colnames(df),phy$BH),]

phy_fat = phy[phy$Overweight==1,]
phy_nor = phy[phy$Overweight==0,]

phy_fat = t(phy_fat[,c('GLU','insulin','TG','FFA','BH')])
phy_nor = t(phy_nor[,c('GLU','insulin','TG','FFA','BH')])

colnames(phy_fat) = phy_fat['BH',]
colnames(phy_nor) = phy_nor['BH',]
df_fat = df[,colnames(df) %in% colnames(phy_fat)]
df_nor = df[,colnames(df) %in% colnames(phy_nor)]

phy_fat = t(phy_fat[-nrow(phy_fat),])
phy_nor = t(phy_nor[-nrow(phy_nor),])
df_fat = t(df_fat)
df_nor = t(df_nor)

rho_df_fat = data.frame()
rho_df_nor = data.frame()
p_df_fat = data.frame()
p_df_nor = data.frame()

for (i in colnames(df_fat)) {
  for (j in colnames(phy_fat)) {
    result = cor.test(as.numeric(df_fat[,i]),as.numeric(phy_fat[,j]),method='spearman')
    rho_df_fat[i,j] = result$estimate["rho"]
    p_df_fat[i,j] = result$p.value
  }
}
rho_df_fat = as.matrix(rho_df_fat)
p_df_fat = as.matrix(p_df_fat)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_fat[i]) < 0.001) {
    p_df_fat[i] = '***'
  } else if (as.numeric(p_df_fat[i]) < 0.01) {
    p_df_fat[i] = '**'
  } else if (as.numeric(p_df_fat[i]) < 0.05) {
    p_df_fat[i] = '*'
  } else {
    p_df_fat[i] = ''
  }
}

for (i in colnames(df_nor)) {
  for (j in colnames(phy_nor)) {
    result = cor.test(as.numeric(df_nor[,i]),as.numeric(phy_nor[,j]),method='spearman')
    rho_df_nor[i,j] = result$estimate["rho"]
    p_df_nor[i,j] = result$p.value
  }
}
rho_df_nor = as.matrix(rho_df_nor)
p_df_nor = as.matrix(p_df_nor)
for (i in 1:length(p_df_fat)) {
  if (as.numeric(p_df_nor[i]) < 0.001) {
    p_df_nor[i] = '***'
  } else if (as.numeric(p_df_nor[i]) < 0.01) {
    p_df_nor[i] = '**'
  } else if (as.numeric(p_df_nor[i]) < 0.05) {
    p_df_nor[i] = '*'
  } else {
    p_df_nor[i] = ''
  }
}

pdf("plasma_E_fat.pdf")
pheatmap(rho_df_fat, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma E Fat', display_numbers=p_df_fat, fontsize_number=10, fontsize=12
)
dev.off()
pdf("plasma_E_nor.pdf")
pheatmap(rho_df_nor, color=colorRampPalette(c('blue','white','red'))(200),
         cellwidth=15, cellheight=15, cluster_rows=F, scale='none', show_rownames = FALSE,
         cluster_cols=F, breaks=seq(-1,1,0.01), legend_breaks=c(-1,0,1),
         main='Plasma E Normal', display_numbers=p_df_nor, fontsize_number=10, fontsize=12
)
dev.off()


##########################################################################################
library(readxl)
library(pheatmap)

# 定义处理A数据的函数
process_ABCED_data <- function(abc_data_path, trans_file, physiological_file, output_pdf) {
  # 加载数据
  load(abc_data_path)
  trans = read_xlsx(trans_file)
  genes = c('FTH1','PF4V1','RAP1B','LMNA',
            'ARPC2','FTL','CA2','PPBP',
            'SRGN','THBS1','APOA5','SOD1',
            'TXN','IDH1','CAT','PRDX1',
            'PRDX2','S100A7','ORM1','MIF',
            'S100A4','CRP','CA1','FLNA')
  
  # 处理trans数据
  trans2 = trans[trans$Gene %in% genes,]
  df = common_combat_urine[rownames(common_combat_urine) %in% trans2$From,]
  df = df[,grep('^A',colnames(df))]
  rownames(df) = trans$Gene[match(rownames(df),trans$From)]
  
  # 处理physiological数据
  physiological = read_xls(physiological_file)
  physiological$Serialnumber = sprintf('%02d', physiological$Serialnumber)
  physiological$Time = sprintf('%02d', physiological$Time)
  physiological$BH = paste0('ABH', physiological$Serialnumber, 'SJ', physiological$Time)
  physiological = physiological[match(colnames(df), physiological$BH),]
  physiological = t(physiological[,c('GLU','insulin','TG','FFA')])
  colnames(physiological) = colnames(df)
  
  # 转置数据
  df = t(df)
  physiological = t(physiological)
  
  # 计算Spearman相关系数
  rho_df = data.frame()
  p_df = data.frame()
  
  for (i in colnames(df)) {
    for (j in colnames(physiological)) {
      result = cor.test(df[,i], physiological[,j], method='spearman')
      rho_df[i,j] = result$estimate["rho"]
      p_df[i,j] = result$p.value
    }
  }
  rho_df = as.matrix(rho_df)
  p_df = as.matrix(p_df)
  p_df = ifelse(p_df < 0.05, '*', '')
  
  # 按字母从A-Z排序行名
  rownames(rho_df) = sort(rownames(rho_df))
  rownames(p_df) = sort(rownames(p_df))
  
  # 确保行名排序的一致性
  rho_df = rho_df[rownames(rho_df), ]
  p_df = p_df[rownames(p_df), ]
  
  # 绘制热图
  pdf(output_pdf)
  pheatmap(rho_df, color=colorRampPalette(c('blue','white','red'))(200),
           cellwidth=15, cellheight=15, cluster_rows=F, scale='none',
           cluster_cols=F, breaks=seq(-1, 1, 0.01), legend_breaks=c(-1,0,1),
           main='Urine A', display_numbers=p_df, fontsize_number=10, fontsize=12)
  dev.off()
}

# 使用示例
process_ABCED_data("combat_urine_FDRnotsel_NA80_omitb36_seqknn_sample.Rdata", 
                   "ID_trans.xlsx", 
                   'clinical_data_A.xls', 
                   "urine_A.pdf")

#############################################
library(readxl)
library(pheatmap)
library(foreign)  # 用于读取.sav文件

# 定义处理BCDE数据的函数
process_BCDE_data <- function(abc_data_path, trans_file, physiological_file, output_pdf) {
  # 加载数据
  load(abc_data_path)
  trans = read_xlsx(trans_file)
  genes = c('FTH1','PF4V1','RAP1B','LMNA',
            'ARPC2','FTL','CA2','PPBP',
            'SRGN','THBS1','APOA5','SOD1',
            'TXN','IDH1','CAT','PRDX1',
            'PRDX2','S100A7','ORM1','MIF',
            'S100A4','CRP','CA1','FLNA')
  
  # 处理trans数据
  trans2 = trans[trans$Gene %in% genes,]
  df = common_combat_urine[rownames(common_combat_urine) %in% trans2$From,]
  df = df[, grep('^B', colnames(df))]
  rownames(df) = trans$Gene[match(rownames(df), trans$From)]
  
  # 处理physiological数据
  physiological = read.spss(physiological_file, to.data.frame = TRUE)
  physiological$SerialNumber = sprintf('%02d', as.numeric(physiological$SerialNumber))
  physiological$Time = sprintf('%02d', physiological$Time)
  physiological$BH = paste0('BBH', physiological$SerialNumber, 'SJ', physiological$Time)
  physiological = physiological[match(colnames(df), physiological$BH),]
  physiological = t(physiological[, c('GLU', 'insulin', 'TG', 'FFA')])
  colnames(physiological) = colnames(df)
  
  # 转置数据
  df = t(df)
  physiological = t(physiological)
  
  # 计算Spearman相关系数
  rho_df = data.frame()
  p_df = data.frame()
  
  for (i in colnames(df)) {
    for (j in colnames(physiological)) {
      result = cor.test(df[, i], physiological[, j], method = 'spearman')
      rho_df[i, j] = result$estimate["rho"]
      p_df[i, j] = result$p.value
    }
  }
  rho_df = as.matrix(rho_df)
  p_df = as.matrix(p_df)
  p_df = ifelse(p_df < 0.05, '*', '')
  
  # 按字母从A到E排序行名
  sorted_row_names = sort(rownames(rho_df))
  rho_df = rho_df[sorted_row_names, ]
  p_df = p_df[sorted_row_names, ]
  
  # 绘制热图
  pdf(output_pdf)
  pheatmap(rho_df, color=colorRampPalette(c('blue', 'white', 'red'))(200),
           cellwidth=15, cellheight=15, cluster_rows=FALSE, scale='none',
           cluster_cols=FALSE, breaks=seq(-1, 1, 0.01), legend_breaks=c(-1, 0, 1),
           main='Urine B', display_numbers=p_df, fontsize_number=10, fontsize=12)
  dev.off()
}

# 使用示例
process_BCDE_data("combat_urine_FDRnotsel_NA80_omitb36_seqknn_sample.Rdata", 
                  "ID_trans.xlsx", 
                  'clinical_data_B.sav', 
                  "urine_B.pdf")

# 使用示例
process_BCDE_data("combat_urine_FDRnotsel_NA80_omitb36_seqknn_sample.Rdata", 
                  "ID_trans.xlsx", 
                  'clinical_data_C.sav', 
                  "urine_C.pdf")

library(readxl)
library(pheatmap)
library(foreign)  # 用于读取.sav文件

# 定义处理D和E数据的函数
process_DE_data <- function(abc_data_path, trans_file, physiological_file, output_pdf, dataset_letter) {
  # 加载数据
  load(abc_data_path)
  trans = read_xlsx(trans_file)
  genes = c('FTH1','PF4V1','RAP1B','LMNA',
            'ARPC2','FTL','CA2','PPBP',
            'SRGN','THBS1','APOA5','SOD1',
            'TXN','IDH1','CAT','PRDX1',
            'PRDX2','S100A7','ORM1','MIF',
            'S100A4','CRP','CA1','FLNA')
  
  # 处理trans数据
  trans2 = trans[trans$Gene %in% genes,]
  df = common_combat_urine[rownames(common_combat_urine) %in% trans2$From,]
  df = df[, grep(paste0('^', dataset_letter), colnames(df))]
  rownames(df) = trans$Gene[match(rownames(df), trans$From)]
  
  # 处理physiological数据
  physiological = read.spss(physiological_file, to.data.frame = TRUE)
  physiological$serialnumber = sprintf('%02d', as.numeric(physiological$serialnumber))
  physiological$Time = sprintf('%02d', physiological$Time)
  physiological$BH = paste0(dataset_letter, 'BH', physiological$serialnumber, 'SJ', physiological$Time)
  physiological = physiological[match(colnames(df), physiological$BH),]
  physiological = t(physiological[, c('GLU', 'insulin', 'TG', 'FFA')])
  colnames(physiological) = colnames(df)
  
  # 转置数据
  df = t(df)
  physiological = t(physiological)
  
  # 计算Spearman相关系数
  rho_df = data.frame()
  p_df = data.frame()
  
  for (i in colnames(df)) {
    for (j in colnames(physiological)) {
      # 删除含有NA的行
      valid_data = na.omit(data.frame(df[, i], physiological[, j]))
      if (nrow(valid_data) >= 2) {  # 至少需要2个有效值来计算相关系数
        result = cor.test(valid_data[, 1], valid_data[, 2], method = 'spearman')
        rho_df[i, j] = result$estimate["rho"]
        p_df[i, j] = result$p.value
      } else {
        rho_df[i, j] = NA
        p_df[i, j] = NA
      }
    }
  }
  
  rho_df = as.matrix(rho_df)
  p_df = as.matrix(p_df)
  p_df = ifelse(is.na(p_df), '', ifelse(p_df < 0.05, '*', ''))
  
  # 按字母从A到E排序行名
  sorted_row_names = sort(rownames(rho_df))
  rho_df = rho_df[sorted_row_names, ]
  p_df = p_df[sorted_row_names, ]
  
  # 绘制热图
  pdf(output_pdf)
  pheatmap(rho_df, color=colorRampPalette(c('blue', 'white', 'red'))(200),
           cellwidth=15, cellheight=15, cluster_rows=FALSE, scale='none',
           cluster_cols=FALSE, breaks=seq(-1, 1, 0.01), legend_breaks=c(-1, 0, 1),
           main=paste('Urine', dataset_letter), display_numbers=p_df, fontsize_number=10, fontsize=12)
  dev.off()
}

# 使用示例
process_DE_data("combat_urine_FDRnotsel_NA80_omitb36_seqknn_sample.Rdata", 
                "ID_trans.xlsx", 
                'clinical_data_D.sav', 
                "urine_D.pdf", 
                "D")

library(readxl)
library(pheatmap)
library(foreign)  # 用于读取.sav文件

# 定义处理E数据的函数
process_E_data <- function(abc_data_path, trans_file, physiological_file, output_pdf) {
  # 加载数据
  load(abc_data_path)
  trans = read_xlsx(trans_file)
  genes = c('FTH1','PF4V1','RAP1B','LMNA',
            'ARPC2','FTL','CA2','PPBP',
            'SRGN','THBS1','APOA5','SOD1',
            'TXN','IDH1','CAT','PRDX1',
            'PRDX2','S100A7','ORM1','MIF',
            'S100A4','CRP','CA1','FLNA')
  
  # 处理trans数据
  trans2 = trans[trans$Gene %in% genes,]
  df = common_combat_urine[rownames(common_combat_urine) %in% trans2$From,]
  df = df[, grep('^E', colnames(df))]
  rownames(df) = trans$Gene[match(rownames(df), trans$From)]
  
  # 处理physiological数据
  physiological = read.spss(physiological_file, to.data.frame = TRUE)
  physiological$Serialnumber = sprintf('%02d', as.numeric(physiological$Serialnumber))
  physiological$Time = sprintf('%02d', physiological$Time)
  physiological$BH = paste0('EBH', physiological$Serialnumber, 'SJ', physiological$Time)
  
  # 打印检查数据
  print("physiological BH column:")
  print(head(physiological$BH))
  
  # 确保数据匹配
  physiological = physiological[match(colnames(df), physiological$BH),]
  
  # 确保数据完整性
  if (any(is.na(physiological))) {
    print("Warning: NA values detected in physiological data")
  }
  
  physiological = t(physiological[, c('GLU', 'insulin', 'TG', 'FFA')])
  colnames(physiological) = colnames(df)
  
  # 转置数据
  df = t(df)
  physiological = t(physiological)
  
  # 计算Spearman相关系数
  rho_df = data.frame()
  p_df = data.frame()
  
  for (i in colnames(df)) {
    for (j in colnames(physiological)) {
      # 删除含有NA的行
      valid_data = na.omit(data.frame(df[, i], physiological[, j]))
      if (nrow(valid_data) >= 2) {  # 至少需要2个有效值来计算相关系数
        result = cor.test(valid_data[, 1], valid_data[, 2], method = 'spearman')
        rho_df[i, j] = result$estimate["rho"]
        p_df[i, j] = result$p.value
      } else {
        rho_df[i, j] = NA
        p_df[i, j] = NA
      }
    }
  }
  
  rho_df = as.matrix(rho_df)
  p_df = as.matrix(p_df)
  p_df = ifelse(is.na(p_df), '', ifelse(p_df < 0.05, '*', ''))
  
  # 按字母从A到E排序行名
  sorted_row_names = sort(rownames(rho_df))
  rho_df = rho_df[sorted_row_names, ]
  p_df = p_df[sorted_row_names, ]
  
  # 绘制热图
  pdf(output_pdf)
  pheatmap(rho_df, color=colorRampPalette(c('blue', 'white', 'red'))(200),
           cellwidth=15, cellheight=15, cluster_rows=FALSE, scale='none',
           cluster_cols=FALSE, breaks=seq(-1, 1, 0.01), legend_breaks=c(-1, 0, 1),
           main='Urine E', display_numbers=p_df, fontsize_number=10, fontsize=12)
  dev.off()
}

# 使用示例
process_E_data("combat_urine_FDRnotsel_NA80_omitb36_seqknn_sample.Rdata", 
               "ID_trans.xlsx", 
               'clinical_data_E.sav', 
               "urine_E.pdf")

