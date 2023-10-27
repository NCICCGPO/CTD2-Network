# This code replicate ORF screening plots from 
#  "Funtional screenings of amplification outlier oncogenes in organoid models of early tumorigenesis"
# Salahudeen et al.
# Fig 3
# Sup Fig S5, S6, S7
# joseaseoane@vhio.net


library("xlsx")
library(ggpubr)
library(ggrastr)


# esophageal

barcode_counts <- read.delim("../data/Barcode_counts_esophageal.txt",as.is=TRUE)
names(barcode_counts) <- gsub("."," ",names(barcode_counts),fixed=TRUE)
sns <- names(barcode_counts)[-(1:5)]
barcode_counts2 <- as.matrix(barcode_counts[,sns])


write.xlsx(barcode_counts,file="counts_propts_v3.xlsx",sheetName = "ESO_counts",append = F)

barcode_prop <- barcode_counts2[-1,] / rep(colSums(barcode_counts2[-1,]),each=nrow(barcode_counts2)-1)
rownames(barcode_prop) <- barcode_counts$Gene[-1]

write.xlsx(barcode_prop,file="counts_propts_v3.xlsx",sheetName = "ESO_props",append = T)


barcode_table <- matrix(0,ncol=ncol(barcode_counts2),nrow=3)
colnames(barcode_table) <- colnames(barcode_counts2)
rownames(barcode_table) <- c("Other barcodes","eGFP","No barcode")
barcode_table[1,] <- colSums(barcode_counts2[!barcode_counts$Gene %in% c("","eGFP"),])
barcode_table[2,] <- barcode_counts2[barcode_counts$Gene == "eGFP",]
barcode_table[3,] <- barcode_counts2[1,] - barcode_table[1,] - barcode_table[2,]

barcode_table2 <- matrix(0,nrow=2,ncol=ncol(barcode_counts2))
colnames(barcode_table2) <- colnames(barcode_counts2)
rownames(barcode_table2) <- c("Barcodes >= 0.1% reads","Barcodes >= 1 read")
barcode_table2[1,] <- colSums(barcode_prop>=0.001)
barcode_table2[2,] <- colSums(barcode_counts2[-1,]>0)

barcode_counts <- barcode_counts[-1,]
sns_m <- matrix(sns,ncol=3)

barcode_counts2 <- barcode_counts[,sns]
barcode_total2 <- colSums(barcode_counts2)

barcode_prop <- as.matrix(barcode_counts2/rep(barcode_total2,each=nrow(barcode_counts2)))


# threshold to remove ORFs with a proportion of ORF/tota llibrary
thrshold_prop = 0.0003



idx1 <- apply(barcode_prop[,sns_m[,1]] >= thrshold_prop ,1,all)

barcode_counts_filtered_esca = barcode_counts[idx1,]
save(barcode_counts_filtered_esca,file="barcode_counts_filtered_esca_dic.RData")

#### esophagel F


toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns)])
toPlot$sample=gsub("[^A-D]","",toPlot$variable,fixed = F)
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns_m[,3]),"T55",
                     ifelse(toPlot$variable %in% c(sns_m[,1]),"T0",NA))


toPlot2d = expand.grid(A1=sns_m[,3],A2=sns_m[,3])

toPlot3d = merge(toPlot2d,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4d = merge(toPlot3d,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)

toPlot4d$A1 = droplevels(toPlot4d$A1)
toPlot4d$A2 = droplevels(toPlot4d$A2)

toPlot4d$filt = toPlot4d$CloneID %in% barcode_counts_filtered_esca$CloneID

ggplot(na.omit(toPlot4d[which(toPlot4d$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+
  facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_esophageal_v6_F.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4d[which(toPlot4d$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+yscale("log10", .format = TRUE)+xscale("log10", .format = TRUE)+
  facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_esophageal_v6_F_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns)])
toPlot$sample=gsub("[^A-D]","",toPlot$variable,fixed = F)
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns_m[,3]),"T55",
                     ifelse(toPlot$variable %in% c(sns_m[,1]),"T0",NA))

ggpubr::ggpaired(na.omit(toPlot[which(toPlot$Gene=="FGF3"),]),
                 x="comp",y="logvalue",id = "sample",facet.by = "CloneID")+ggtitle("FGF3 F")+ylab("log10 counts")
ggsave("paired_esohpageal_v6_F.pdf",width = 7,height = 7)


toPlot2 = expand.grid(T0=sns_m[,1],T55=sns_m[,3])

toPlot3 = merge(toPlot2,toPlot,by = c("T0"),by.y=c("variable"),
                all = T)

toPlot4 = merge(toPlot3,toPlot,by = c("CloneID","Gene","T55"),by.y=c("CloneID","Gene","variable"),
                all = T)

toPlot4$filt = toPlot4$CloneID %in% barcode_counts_filtered_esca$CloneID

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T55)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_esophageal_v6_T0_F.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+yscale("log10", .format = TRUE)+xscale("log10", .format = TRUE)+
  facet_grid(T0~T55)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_esophageal_v6_T0_F_log.pdf",width = 10,height = 10, useDingbats=FALSE)


# boxplots

control1 <- barcode_counts$Control == "Yes" & idx1
barcode_counts3 <- barcode_counts2[idx1,]
barcode_prop3 <- barcode_prop[idx1,]

barcode_test_en = rep(0,nrow(barcode_prop3))
for(i in 1:nrow(barcode_prop3)){
  barcode_test_en[i]= t.test(log(barcode_prop3[i,sns_m[,2]]+1E-8),log(barcode_prop3[i,sns_m[,1]]+1E-8),alternative="g",paired=T)$p.value
}


barcode_test_F = rep(0,nrow(barcode_prop3))
for(i in 1:nrow(barcode_prop3)){
  barcode_test_F[i]= t.test(log(barcode_prop3[i,sns_m[,3]]+1E-8),log(barcode_prop3[i,sns_m[,1]]+1E-8),alternative="g",paired = T)$p.value
}
names(barcode_test_en)=names(barcode_test_F)=barcode_counts$Gene[idx1]

pval_df = data.frame(geneName=make.names(barcode_counts$Gene[idx1],unique = T),barcode_test_en=barcode_test_en,barcode_test_F=barcode_test_F,stringsAsFactors = F)
pval_df$signE = pval_df$barcode_test_en<0.05
pval_df$signF = pval_df$barcode_test_F<0.05

ratio1 <- matrix(0,nrow=nrow(barcode_prop3),ncol=nrow(sns_m))
rownames(ratio1) <- barcode_counts$Gene[idx1]
colnames(ratio1) <- sns_m[,2]
for (i in 1:nrow(sns_m)) {
  ratio1[,i] <- barcode_prop3[,sns_m[i,2]]/ barcode_prop3[,sns_m[i,1]]
}
ratio2 <- matrix(0,nrow=nrow(barcode_prop3),ncol=nrow(sns_m))
rownames(ratio2) <- barcode_counts$Gene[idx1]
colnames(ratio2) <- sns_m[,3]
for (i in 1:nrow(sns_m)) {
  ratio2[,i] <- barcode_prop3[,sns_m[i,3]]/ barcode_prop3[,sns_m[i,1]]
}

ratio1[which(is.infinite(ratio1))]=0
ratio1[which(is.nan(ratio1))]=0

ratio2[which(is.infinite(ratio2))]=0
ratio2[which(is.nan(ratio2))]=0

order1 <- order(rownames(ratio1) %in% barcode_counts$Gene[control1],
                apply(ratio1,1,median) , 
                decreasing=TRUE)

ratio1 <- pmax(ratio1, 1e-5)

order2 <- order(rownames(ratio2) %in% barcode_counts$Gene[control1],
                apply(ratio2,1,median) , 
                decreasing=TRUE)

ratio2 <- pmax(ratio2, 1e-5)

toPlot2 = reshape2::melt(data.frame(t(ratio1[order1[1:50],])))
p_eso_E = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,20))+geom_text(data=pval_df,aes(x=geneName,y=20,label=ifelse(signE,signif(barcode_test_en,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Mouse Esophagus TP53 EN T46 vs T0")#+geom_hline(yintercept = 2)

p_eso_E
ggsave("eso_TP53_v5_E.pdf",width = 8.5,height = 4)

toPlot2 = reshape2::melt(as.data.frame(t(ratio2[order2[1:50],])))
p_eso_F =ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,21))+geom_text(data=pval_df,aes(x=geneName,y=21.5,label=ifelse(signF,signif(barcode_test_F,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Mouse Esophagus TP53 F12 T46 vs T0")

p_eso_F
ggsave("eso_TP53_v5_F.pdf",width = 8.5,height = 4)


# oral 
barcode_counts <- read.delim("../data/Barcode_counts_oral_gastric.txt",as.is=TRUE)
names(barcode_counts) <- gsub("."," ",names(barcode_counts),fixed=TRUE)
sns <- names(barcode_counts)[-(1:12)]
barcode_counts2 <- as.matrix(barcode_counts[,sns])

barcode_prop <- barcode_counts2[-1,] / rep(colSums(barcode_counts2[-1,]),each=nrow(barcode_counts2)-1)
rownames(barcode_prop) <- barcode_counts$Gene[-1]


write.xlsx(barcode_counts[,c(1:4,13:24)],file="counts_propts_v3.xlsx",sheetName = "ORAL_counts",append = T)
write.xlsx(barcode_prop,file="counts_propts_v3.xlsx",sheetName = "ORAL_props",append = T)


barcode_table <- matrix(0,ncol=ncol(barcode_counts2),nrow=3)
colnames(barcode_table) <- colnames(barcode_counts2)
rownames(barcode_table) <- c("Other barcodes","eGFP","No barcode")
barcode_table[1,] <- colSums(barcode_counts2[!barcode_counts$Gene %in% c("","eGFP"),])
barcode_table[2,] <- barcode_counts2[barcode_counts$Gene == "eGFP",]
barcode_table[3,] <- barcode_counts2[1,] - barcode_table[1,] - barcode_table[2,]

barcode_table2 <- matrix(0,nrow=2,ncol=ncol(barcode_counts2))
colnames(barcode_table2) <- colnames(barcode_counts2)
rownames(barcode_table2) <- c("Barcodes >= 0.1% reads","Barcodes >= 1 read")
barcode_table2[1,] <- colSums(barcode_prop>=0.001)
barcode_table2[2,] <- colSums(barcode_counts2[-1,]>0)

barcode_counts <- barcode_counts[-1,]
barcode_total <- colSums(barcode_counts[,sns])

barcode_prop <- as.matrix(barcode_counts[,sns]/rep(barcode_total,each=nrow(barcode_counts)))
sns_m <- matrix(sns,ncol=3)


idx2 <- apply(barcode_prop[,sns_m[,1]] >= thrshold_prop,1,all)

barcode_counts_filtered_oral = barcode_counts[idx2,]
save(barcode_counts_filtered_oral,file="barcode_counts_filtered_oral_dic.RData")


# oral pairwise plots K
toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns)])
toPlot$sample=sapply(as.character(toPlot$variable), function(x) substr(x,nchar(x),nchar(x)))
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns_m[,2]),"T52",
                     ifelse(toPlot$variable %in% c(sns_m[,1]),"T0",NA))


ggpubr::ggpaired(toPlot[which(toPlot$Gene=="DYRK2" & toPlot$CloneID=="TRCN0000489007"& !is.na(toPlot$comp)),],
                 x="comp",y="logvalue",id = "sample")+ggtitle("DYRK2 Oral E")+ylab("log10 counts")
ggsave("paired_oral_DYRK2_TRCN489007_v6_EN.pdf",width = 7,height = 7)


toPlot2 = expand.grid(T0=sns_m[,1],T52=sns_m[,2])

toPlot3 = merge(toPlot2,toPlot,by = c("T0"),by.y=c("variable"),
                all = T)

toPlot4 = merge(toPlot3,toPlot,by = c("CloneID","Gene","T52"),by.y=c("CloneID","Gene","variable"),
                all = T)

toPlot4$filt = toPlot4$CloneID %in% barcode_counts_filtered_oral$CloneID

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T52)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_T0_T52.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T52)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_T0_T52_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2b = expand.grid(A1=sns_m[,1],A2=sns_m[,1])

toPlot3b = merge(toPlot2b,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4b = merge(toPlot3b,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)

toPlot4b$filt = toPlot4b$CloneID %in% barcode_counts_filtered_oral$CloneID

ggplot(na.omit(toPlot4b),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_T0.pdf",width = 10,height = 10, useDingbats=FALSE)


ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_T0_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2c = expand.grid(A1=sns_m[,2],A2=sns_m[,2])

toPlot3c = merge(toPlot2c,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4c = merge(toPlot3c,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)
toPlot4c$filt = toPlot4c$CloneID %in% barcode_counts_filtered_oral$CloneID

ggplot(na.omit(toPlot4c),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_EN_T55.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_EN_T55_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


# F
toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns)])
toPlot$sample=sapply(as.character(toPlot$variable), function(x) substr(x,nchar(x),nchar(x)))
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns_m[,3]),"T52",
                     ifelse(toPlot$variable %in% c(sns_m[,1]),"T0",NA))


toPlot2 = expand.grid(T0=sns_m[,1],T52=sns_m[,3])

toPlot3 = merge(toPlot2,toPlot,by = c("T0"),by.y=c("variable"),
                all = T)

toPlot4 = merge(toPlot3,toPlot,by = c("CloneID","Gene","T52"),by.y=c("CloneID","Gene","variable"),
                all = T)
toPlot4$filt = toPlot4$CloneID %in% barcode_counts_filtered_oral$CloneID

ggplot(na.omit(toPlot4),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T52)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v3_K_T0_T52_F12.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T52)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_T0_T52_F12_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2b = expand.grid(A1=sns_m[,3],A2=sns_m[,3])

toPlot3b = merge(toPlot2b,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4b = merge(toPlot3b,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)
toPlot4b$filt = toPlot4b$CloneID %in% barcode_counts_filtered_oral$CloneID

ggplot(na.omit(toPlot4b),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v3_K_T55F.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_oral_v6_K_T55F_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)



control1 <- barcode_counts$Control == "Yes" & idx2
barcode_counts2 <- barcode_counts[idx2,c(colnames(barcode_counts)[1:4],as.vector(sns_m))]
barcode_prop2 <- barcode_prop[idx2,]

ratio2a <- ratio2b <- matrix(0,nrow=nrow(barcode_prop2),ncol=nrow(sns_m))
rownames(ratio2a) <- rownames(ratio2b) <- barcode_counts$Gene[idx2]
colnames(ratio2a) <- sns_m[,2]
colnames(ratio2b) <- sns_m[,3]
for (i in 1:nrow(sns_m)) {
  ratio2a[,i] <- barcode_prop2[,sns_m[i,2]]/ barcode_prop2[,sns_m[i,1]]
  ratio2b[,i] <- barcode_prop2[,sns_m[i,3]]/ barcode_prop2[,sns_m[i,1]]
}


barcode_test_en = rep(0,nrow(barcode_prop2))
for(i in 1:nrow(barcode_prop2)){
  barcode_test_en[i]= t.test(log(barcode_prop2[i,sns_m[,2]]+1E-8),log(barcode_prop2[i,sns_m[,1]]+1E-8),alternative="g",paired=T)$p.value
}


barcode_test_F = rep(0,nrow(barcode_prop2))
for(i in 1:nrow(barcode_prop2)){
  barcode_test_F[i]= t.test(log(barcode_prop2[i,sns_m[,3]]+1E-8),log(barcode_prop2[i,sns_m[,1]]+1E-8),alternative="g",paired=T)$p.value
}
names(barcode_test_en)=names(barcode_test_F)=rownames(ratio2b)

pval_df = data.frame(geneName=make.names(rownames(ratio2b),unique = T),barcode_test_en=barcode_test_en,barcode_test_F=barcode_test_F,stringsAsFactors = F)
pval_df$signE = pval_df$barcode_test_en<0.05
pval_df$signF = pval_df$barcode_test_F<0.05
pval_df$yposE=apply(ratio2a,1,function(x) max(x[is.finite(x)]))
pval_df$yposF=apply(ratio2b,1,function(x) max(x[is.finite(x)]))
rownames(pval_df)=pval_df$geneName
order2a <- order(rownames(ratio2a) %in% barcode_counts$Gene[control1],
                 apply(ratio2a,1,median),decreasing=TRUE)

order2b <- order(rownames(ratio2b) %in% barcode_counts$Gene[control1],
                 apply(ratio2b,1,median),decreasing=TRUE)

ratio2a <- pmax(ratio2a, 1e-5)
ratio2b <- pmax(ratio2b, 1e-5)



toPlot2 = reshape2::melt(data.frame(t(ratio2a[order2a,])))
p_oral_2a =ggpubr::ggboxplot(toPlot2,x="variable",y="value",add = "jitter")+
  coord_cartesian(ylim = c(0,50))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Mouse Oral mucosa TP53 KY0 vs KY52 EN")

p_oral_2a
ggsave("oral_TP53_v6_ky_EN.pdf",width = 8.5,height = 4)

p_oral_2a =ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,51),xlim = c(0,43))+
  geom_text(data=pval_df,aes(x=geneName,y=yposE+2,label=ifelse(signE,signif(barcode_test_en,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Mouse Oral mucosa TP53 KY0 vs KY52 EN")


p_oral_2a
ggsave("oral_TP53_v6_ky_EN_v2.pdf",width = 8.5,height = 4)


# lung

barcode_counts <- read.delim("../data/Barcode_counts_lung_colon.txt",as.is=TRUE)
names(barcode_counts) <- gsub("."," ",names(barcode_counts),fixed=TRUE)
sns1 <- matrix(c("B07","B08","B09","B10","A06","B11","A11","B12"),
               ncol=2,byrow=TRUE)

barcode_counts2 <- as.matrix(barcode_counts[,sns1])

barcode_prop <- barcode_counts2[-1,] / rep(colSums(barcode_counts2[-1,]),each=nrow(barcode_counts2)-1)
rownames(barcode_prop) <- barcode_counts$Gene[-1]


write.xlsx(cbind(barcode_counts[,1:6],barcode_counts2),file="counts_propts_v3.xlsx",sheetName = "LUNG_counts",append = T)
write.xlsx(barcode_prop,file="counts_propts_v3.xlsx",sheetName = "LUNG_props",append = T)

barcode_table <- matrix(0,ncol=ncol(barcode_counts2),nrow=3)
colnames(barcode_table) <- colnames(barcode_counts2)
rownames(barcode_table) <- c("Other barcodes","eGFP","No barcode")
barcode_table[1,] <- colSums(barcode_counts2[!barcode_counts$Gene %in% c("","eGFP"),])
barcode_table[2,] <- barcode_counts2[barcode_counts$Gene == "eGFP",]
barcode_table[3,] <- barcode_counts2[1,] - barcode_table[1,] - barcode_table[2,]


idx1 <- barcode_counts$DOD19_20[-c(1)] == "Yes" &
  apply(barcode_prop[,sns1[,1]] >= thrshold_prop,1,all)

barcode_counts_filtered_lung = barcode_counts[idx1,]
save(barcode_counts_filtered_lung,file="barcode_counts_filtered_lung_dic.RData")

# pairwise stuff
toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns1)])
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns1[,2]),"T4",
                     ifelse(toPlot$variable %in% c(sns1[,1]),"T0",NA))


toPlot2 = expand.grid(T0=sns1[,1],T4=sns1[,2])

toPlot3 = merge(toPlot2,toPlot,by = c("T0"),by.y=c("variable"),
                all = T)

toPlot4 = merge(toPlot3,toPlot,by = c("CloneID","Gene","T4"),by.y=c("CloneID","Gene","variable"),
                all = T)

toPlot4$filt = toPlot4$CloneID %in% barcode_counts_filtered_lung$CloneID

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T4)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_mouseLung_v6_T0_T4.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T4)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_mouseLung_v6_T0_T4_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2b = expand.grid(A1=sns1[,1],A2=sns1[,1])

toPlot3b = merge(toPlot2b,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4b = merge(toPlot3b,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)

toPlot4b$filt = toPlot4b$CloneID %in% barcode_counts_filtered_lung$CloneID

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_mouseLung_v3_T0.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_mouseLung_v6_T0_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)

toPlot2c = expand.grid(A1=sns1[,2],A2=sns1[,2])

toPlot3c = merge(toPlot2c,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4c = merge(toPlot3c,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)
toPlot4c$filt = toPlot4c$CloneID %in% barcode_counts_filtered_lung$CloneID

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_mouseLung_v3_T4.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_mouseLung_v6_T4_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


control1 <- barcode_counts$Control[-c(1)] == "Yes" & idx1
ratio_control1 <- colSums(barcode_counts2[control1,as.vector(sns1)])/colSums(barcode_counts2[,as.vector(sns1)])


ratio1 <- matrix(0,nrow=sum(idx1),ncol=nrow(sns1))
rownames(ratio1) <- barcode_counts$Gene[-c(1)][idx1]
for (i in 1:nrow(sns1)) {
  ratio1[,i] <- barcode_prop[idx1,sns1[i,2]]/ barcode_prop[idx1,sns1[i,1]]
}


ratio2 <- matrix(0,nrow=sum(idx1),ncol=nrow(sns1))
rownames(ratio2) <- barcode_counts$Gene[-c(1)][idx1]
for (i in 1:nrow(sns1)) {
  ratio2[,i] <- (barcode_prop[idx1,sns1[i,2]] / ratio_control1[sns1[i,2]]) /
    (barcode_prop[idx1,sns1[i,1]] / ratio_control1[sns1[i,1]])
}


barcode_prop3 = barcode_prop[idx1,]
barcode_test = rep(0,nrow(barcode_prop3))
for(i in 1:nrow(barcode_prop3)){
  barcode_test[i]= t.test(log(barcode_prop3[i,sns1[,2]]+1E-8),log(barcode_prop3[i,sns1[,1]]+1E-8),alternative="g",paired=T)$p.value
}
names(barcode_test)=rownames(barcode_prop3)

pval_df = data.frame(geneName=make.names(rownames(barcode_prop3),unique = T),barcode_test=barcode_test,stringsAsFactors = F)
pval_df$sign = pval_df$barcode_test<0.05 & pval_df$geneName!="BFP"
pval_df$ypos=apply(ratio1,1,function(x) max(x[is.finite(x)]))



order1 = order(rownames(ratio1) %in% barcode_counts$Gene[-c(1)][control1],apply(ratio1,1,median) ,decreasing=TRUE)

order2 = order(rownames(ratio2) %in% barcode_counts$Gene[-c(1)][control1],apply(ratio2,1,median),decreasing=TRUE)



toPlot2 = reshape2::melt(data.frame(t(ratio1[order1,])))
p_lung_1 = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,8))+geom_text(data=pval_df,aes(x=geneName,y=4,label=ifelse(sign,signif(barcode_test,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Mouse Lung normalized by total barcode read counts TP4 vs T0")

p_lung_1
ggsave("lung_norm_total_v6.pdf",width = 8.5,height = 4)



# pancreas   

barcode_counts <- read.delim("../data/Barcode_counts_pancreas.txt",as.is=TRUE)
sns <- names(barcode_counts)[-(1:5)]
barcode_counts2 <- as.matrix(barcode_counts[,sns])

barcode_prop <- barcode_counts2[-1,] / rep(colSums(barcode_counts2[-1,]),each=nrow(barcode_counts2)-1)
rownames(barcode_prop) <- barcode_counts$Gene[-1]

write.xlsx(barcode_counts,file="counts_propts_v3.xlsx",sheetName = "PANCREAS_counts",append = T)
write.xlsx(barcode_prop,file="counts_propts_v3.xlsx",sheetName = "PANCREAS_props",append = T)

barcode_table <- matrix(0,ncol=ncol(barcode_counts2),nrow=3)
colnames(barcode_table) <- colnames(barcode_counts2)
rownames(barcode_table) <- c("Other barcodes","eGFP","No barcode")
barcode_table[1,] <- colSums(barcode_counts2[!barcode_counts$Gene %in% c("","eGFP"),])
barcode_table[2,] <- barcode_counts2[barcode_counts$Gene == "eGFP",]
barcode_table[3,] <- barcode_counts2[1,] - barcode_table[1,] - barcode_table[2,]

barcode_table2 <- matrix(0,nrow=2,ncol=ncol(barcode_counts2))
colnames(barcode_table2) <- colnames(barcode_counts2)
rownames(barcode_table2) <- c("Barcodes >= 0.1% total reads","Barcodes >= 0.01% total reads")
barcode_table2[1,] <- colSums(barcode_prop>=0.001)
barcode_table2[2,] <- colSums(barcode_prop>=0.0001)

sns_m <- matrix(sns[c(4,5,6,1,2,3)],ncol=2)

barcode_counts <- barcode_counts[-1,]
barcode_total <- colSums(barcode_counts[,sns])
barcode_counts <- barcode_counts[barcode_counts$PAAD == "Yes",]

barcode_prop <- as.matrix(barcode_counts[,sns]/rep(barcode_total,each=nrow(barcode_counts)))
rownames(barcode_prop) <- barcode_counts$Gene


idx1 <- apply(barcode_prop[,sns_m[,1]] >= thrshold_prop,1,all)

barcode_counts_filtered_pancreas = barcode_counts[idx1,]
save(barcode_counts_filtered_pancreas,file="barcode_counts_filtered_pancreas_dic.RData")

# pairwise stuff
toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns)])
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns_m[,2]),"T33",
                     ifelse(toPlot$variable %in% c(sns_m[,1]),"T0",NA))


toPlot2 = expand.grid(T0=sns_m[,1],T33=sns_m[,2])

toPlot3 = merge(toPlot2,toPlot,by = c("T0"),by.y=c("variable"),
                all = T)

toPlot4 = merge(toPlot3,toPlot,by = c("CloneID","Gene","T33"),by.y=c("CloneID","Gene","variable"),
                all = T)
toPlot4$filt = toPlot4$CloneID %in% barcode_counts_filtered_pancreas$CloneID

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T33)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_pancreas_v6_T0_T33.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T33)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_pancreas_v6_T0_T33_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2b = expand.grid(A1=sns_m[,1],A2=sns_m[,1])

toPlot3b = merge(toPlot2b,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4b = merge(toPlot3b,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)
toPlot4b$filt = toPlot4b$CloneID %in% barcode_counts_filtered_pancreas$CloneID

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_pancreas_v6_T0.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_pancreas_v6_T0_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2c = expand.grid(A1=sns_m[,2],A2=sns_m[,2])

toPlot3c = merge(toPlot2c,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4c = merge(toPlot3c,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)

toPlot4c$filt = toPlot4c$CloneID %in% barcode_counts_filtered_pancreas$CloneID

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_pancreas_v6_T33.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_pancreas_v6_T33_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


control1 <- barcode_counts$Control == "Yes" & idx1
barcode_counts1 <- barcode_counts[idx1,c(colnames(barcode_counts)[1:5],as.vector(sns_m[,1:2]))]
barcode_prop1 <- barcode_prop[idx1,]


ratio1 <- matrix(0,nrow=nrow(barcode_prop1),ncol=nrow(sns_m))
rownames(ratio1) <- barcode_counts$Gene[idx1] #rownames(barcode_prop1)
colnames(ratio1) <- sns_m[,2]
for (i in 1:nrow(sns_m)) {
  ratio1[,i] <- barcode_prop1[,sns_m[i,2]]/ barcode_prop1[,sns_m[i,1]]
}

barcode_test = rep(0,nrow(barcode_prop1))
for(i in 1:nrow(barcode_prop1)){
  barcode_test[i]= t.test(log(barcode_prop1[i,sns_m[,2]]+1E-8),log(barcode_prop1[i,sns_m[,1]]+1E-8),alternative="g",paired=T)$p.value
}
names(barcode_test)=rownames(barcode_prop1)

pval_df = data.frame(geneName=make.names(rownames(ratio1),unique = T),barcode_test=barcode_test,stringsAsFactors = F)
pval_df$sign = pval_df$barcode_test<0.05
pval_df$ypos=apply(ratio1,1,function(x) max(x[is.finite(x)]))



order1 <- order(rownames(ratio1) %in% barcode_counts$Gene[control1],
                apply(ratio1,1,median), ##+ apply(ratio2,1,median), 
                decreasing=TRUE)

toPlot2 = reshape2::melt(data.frame(t(ratio1[order1,])))
p_pancreas = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,10))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Mouse Pancreast KRAS  T0 vs T33")

p_pancreas
ggsave("pancreas_v3_KRAS.pdf",width = 8.5,height = 4)


p_pancreas = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,10))+geom_text(data=pval_df,aes(x=geneName,y=4,label=ifelse(sign,signif(barcode_test,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Mouse Pancreast KRAS  T0 vs T33")

p_pancreas
ggsave("pancreas_v5_KRAS.pdf",width = 8.5,height = 4)



# Human gastric p53 R175h
# 20170511
barcode_counts <- read.delim("../data/Barcode_counts_oral_gastric.txt",as.is=TRUE)
sns <- names(barcode_counts)[5:12]
barcode_counts2 <- as.matrix(barcode_counts[,sns])

barcode_prop <- barcode_counts2[-1,] / rep(colSums(barcode_counts2[-1,]),each=nrow(barcode_counts2)-1)
rownames(barcode_prop) <- barcode_counts$Gene[-1]

write.xlsx(cbind(barcode_counts[,1:4],barcode_counts2),file="counts_propts_v3.xlsx",sheetName = "GASTRIC_R175H_counts",append = T)
write.xlsx(barcode_prop,file="counts_propts_v3.xlsx",sheetName = "GASTRIC_R175H_props",append = T)



barcode_table <- matrix(0,ncol=ncol(barcode_counts2),nrow=3)
colnames(barcode_table) <- colnames(barcode_counts2)
rownames(barcode_table) <- c("Other barcodes","eGFP","No barcode")
barcode_table[1,] <- colSums(barcode_counts2[!barcode_counts$Gene %in% c("","eGFP"),])
barcode_table[2,] <- barcode_counts2[barcode_counts$Gene == "eGFP",]
barcode_table[3,] <- barcode_counts2[1,] - barcode_table[1,] - barcode_table[2,]

barcode_table2 <- matrix(0,nrow=2,ncol=ncol(barcode_counts2))
colnames(barcode_table2) <- colnames(barcode_counts2)
rownames(barcode_table2) <- c("Barcodes >= 0.1% total reads","Barcodes >= 0.01% total reads")
barcode_table2[1,] <- colSums(barcode_prop>=0.001)
barcode_table2[2,] <- colSums(barcode_prop>=0.0001)

sns_m <- matrix(sns[c(1,2,3,4,5,6,7,8)],ncol=2)


idx1 <- apply(barcode_prop[,sns_m[,1]] >= thrshold_prop,1,all)

barcode_counts_filtered_gastric = barcode_counts[idx1,]
save(barcode_counts_filtered_gastric,file="barcode_counts_filtered_gastric_dic.RData")


# GASTRIC pairwise plots Amanda
toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns)])
toPlot$sample=sapply(as.character(toPlot$variable), function(x) substr(x,nchar(x),nchar(x)))
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns_m[,2]),"T55",
                     ifelse(toPlot$variable %in% c(sns_m[,1]),"T0",NA))

ggpubr::ggpaired(na.omit(toPlot[which(toPlot$Gene=="CDK6"),]),
                 x="comp",y="logvalue",id = "sample",facet.by = "CloneID")+ggtitle("CDK6 Gastric")+ylab("log10 counts")
ggsave("paired_gastric_v3_cdk6.pdf",width = 7,height = 7)


toPlot2 = expand.grid(T0=sns_m[,1],T55=sns_m[,2])

toPlot3 = merge(toPlot2,toPlot,by = c("T0"),by.y=c("variable"),
                all = T)

toPlot4 = merge(toPlot3,toPlot,by = c("CloneID","Gene","T55"),by.y=c("CloneID","Gene","variable"),
                all = T)

toPlot4$filt = toPlot4$CloneID %in% barcode_counts_filtered_gastric$CloneID

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T55)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_gastric_v6_T0_T55.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4[which(toPlot4$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T55)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_gastric_v6_T0_T55_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)

toPlot2b = expand.grid(A1=sns_m[,1],A2=sns_m[,1])

toPlot3b = merge(toPlot2b,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4b = merge(toPlot3b,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)

toPlot4b$filt = toPlot4b$CloneID %in% barcode_counts_filtered_gastric$CloneID

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_gastric_v6_AT0.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+stat_cor()
ggsave("pairwise_gastric_v6_AT0_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2c = expand.grid(A1=sns_m[,2],A2=sns_m[,2])

toPlot3c = merge(toPlot2c,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4c = merge(toPlot3c,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)
toPlot4c$filt = toPlot4c$CloneID %in% barcode_counts_filtered_gastric$CloneID

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_gastric_v6_AT55.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_gastric_v6_AT55_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)

control1 <- barcode_counts$Control[-c(1)] == "Yes" & idx1
barcode_prop1 <- barcode_prop[idx1,]


ratio1 <- matrix(0,nrow=nrow(barcode_prop1),ncol=nrow(sns_m))
rownames(ratio1) <- rownames(barcode_prop1)
colnames(ratio1) <- sns_m[,2]
for (i in 1:nrow(sns_m)) {
  ratio1[,i] <- barcode_prop1[,sns_m[i,2]]/ barcode_prop1[,sns_m[i,1]]
}

barcode_test = rep(0,nrow(barcode_prop1))
for(i in 1:nrow(barcode_prop1)){
  barcode_test[i]= t.test(log(barcode_prop1[i,sns_m[,2]]+1E-8),log(barcode_prop1[i,sns_m[,1]]+1E-8),alternative="g",paired=T)$p.value
}
names(barcode_test)=rownames(barcode_prop1)

pval_df = data.frame(geneName=make.names(rownames(barcode_prop1),unique = T),barcode_test=barcode_test,stringsAsFactors = F)
pval_df$sign = pval_df$barcode_test<0.05
pval_df$ypos=apply(ratio1,1,function(x) max(x[is.finite(x)]))




order1 <- order(rownames(ratio1) %in% barcode_counts$Gene[-c(1)][control1],
                apply(ratio1,1,median), ##+ apply(ratio2,1,median), 
                decreasing=TRUE)

toPlot2 = reshape2::melt(as.data.frame(t(ratio1[order1[1:50],])))
p_gastric_p53_R175 = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,10))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Human Gastric p53 R175H  T0 vs T55")

p_gastric_p53_R175
ggsave("gastric_p53_R175_v3.pdf",width = 8.5,height = 4)

toPlot2 = reshape2::melt(data.frame(t(ratio1[order1[1:50],])))
p_gastric_p53_R175 = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,10))+geom_text(data=pval_df[which(pval_df$geneName %in% toPlot2$variable),],aes(x=geneName,y=ypos+1,label=ifelse(sign,signif(barcode_test,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Human Gastric p53 R175H  T0 vs T55")

p_gastric_p53_R175
ggsave("gastric_p53_R175_v5.pdf",width = 8.5,height = 4)


#colon 
#c("A09","B05","A09","B03","A10","B04","A10","B06")

barcode_counts <- read.delim("../data/Barcode_counts_lung_colon.txt",as.is=TRUE)

sns1 <- matrix(c("A09","B05","B03","A10","B04","B06"),
               ncol=3,byrow=TRUE)

barcode_counts2 <- as.matrix(barcode_counts[,sns1])

barcode_prop <- barcode_counts2[-1,] / rep(colSums(barcode_counts2[-1,]),each=nrow(barcode_counts2)-1)
rownames(barcode_prop) <- barcode_counts$Gene[-1]

barcode_table <- matrix(0,ncol=ncol(barcode_counts2),nrow=3)
colnames(barcode_table) <- colnames(barcode_counts2)
rownames(barcode_table) <- c("Other barcodes","eGFP","No barcode")
barcode_table[1,] <- colSums(barcode_counts2[!barcode_counts$Gene %in% c("","eGFP"),])
barcode_table[2,] <- barcode_counts2[barcode_counts$Gene == "eGFP",]
barcode_table[3,] <- barcode_counts2[1,] - barcode_table[1,] - barcode_table[2,]

write.xlsx(cbind(barcode_counts[,1:6],barcode_counts2),file="counts_propts_v3.xlsx",sheetName = "COLON_counts",append = T)
write.xlsx(barcode_prop,file="counts_propts_v3.xlsx",sheetName = "COLON_props",append = T)

idx1 <- barcode_counts$COAD[-c(1)] == "Yes" &
  apply(barcode_prop[,sns1[,1]] >= 0.00003,1,all)

barcode_counts_filtered_colon = barcode_counts[idx1,]
save(barcode_counts_filtered_colon,file="barcode_counts_filtered_colon_dic.RData")

# plot pairwise

toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns1)])
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns1[,2]),"T1",
                     ifelse(toPlot$variable %in% c(sns1[,1]),"T0",NA))

toPlot = reshape2::melt(barcode_counts[,c("CloneID","Gene",sns1)])
toPlot$logvalue= log10(toPlot$value)
toPlot$comp = ifelse(toPlot$variable %in% c(sns1[,3]),"T1",
                     ifelse(toPlot$variable %in% c(sns1[,1]),"T0",NA))

ggpubr::ggpaired(na.omit(toPlot[which(toPlot$Gene=="KRAS"),]),
                 x="comp",y="logvalue",facet.by = "CloneID")+ggtitle("KRAS")+ylab("log10 counts")
ggsave("paired_COLON_v3_KRAS_T0_T2.pdf",width = 7,height = 7)


#pairwise
toPlot2 = expand.grid(T0=sns1[,1],T1=sns1[,3])

toPlot3 = merge(toPlot2,toPlot,by = c("T0"),by.y=c("variable"),
                all = T)

toPlot4 = merge(toPlot3,toPlot,by = c("CloneID","Gene","T1"),by.y=c("CloneID","Gene","variable"),
                all = T)
toPlot4$filt = toPlot4$CloneID %in% barcode_counts_filtered_colon$CloneID

ggplot(na.omit(toPlot4[which(toPlot4$filt ),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T1)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_colon_v6_T0_T2.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4[which(toPlot4$filt ),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(T0~T1)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_colon_v6_T0_T2_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


toPlot2b = expand.grid(A1=sns1[,1],A2=sns1[,1])

toPlot3b = merge(toPlot2b,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4b = merge(toPlot3b,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)

toPlot4b$filt = toPlot4b$CloneID %in% barcode_counts_filtered_colon$CloneID
ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_colon_v6_AT0.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4b[which(toPlot4b$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_colon_v6_AT0_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)

toPlot2c = expand.grid(A1=sns1[,3],A2=sns1[,3])

toPlot3c = merge(toPlot2c,toPlot,by = c("A1"),by.y=c("variable"),
                 all = T)

toPlot4c = merge(toPlot3c,toPlot,by = c("CloneID","Gene","A2"),by.y=c("CloneID","Gene","variable"),
                 all = T)

toPlot4c$A1 = droplevels(toPlot4c$A1)
toPlot4c$A2 = droplevels(toPlot4c$A2)

toPlot4c$filt = toPlot4c$CloneID %in% barcode_counts_filtered_colon$CloneID

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_colon_v6_AT2.pdf",width = 10,height = 10, useDingbats=FALSE)

ggplot(na.omit(toPlot4c[which(toPlot4c$filt),]),aes(x=value.x,y=value.y))+geom_point_rast()+facet_grid(A1~A2)+xscale("log10", .format = TRUE)+yscale("log10", .format = TRUE)+
  theme_bw()+xlab("counts")+ylab("counts")+ggpubr::stat_cor()
ggsave("pairwise_colon_v6_AT2_vlog.pdf",width = 10,height = 10, useDingbats=FALSE)


control1 <- barcode_counts$Control[-c(1)] == "Yes" & idx1
ratio_control1 <- colSums(barcode_counts2[control1,unique(as.vector(sns1)),drop=FALSE])/colSums(barcode_counts2[,unique(as.vector(sns1))])

ratio1 <- matrix(0,nrow=sum(idx1),ncol=nrow(sns1))
rownames(ratio1) <- barcode_counts$Gene[-c(1)][idx1]
for (i in 1:nrow(sns1)) {
  ratio1[,i] <- barcode_prop[idx1,sns1[i,2]]/ barcode_prop[idx1,sns1[i,1]]
}
order1 = order(rownames(ratio1) %in% barcode_counts$Gene[-c(1)][control1],apply(ratio1,1,median),decreasing=TRUE)

ratio2 <- matrix(0,nrow=sum(idx1),ncol=nrow(sns1))
rownames(ratio2) <- barcode_counts$Gene[-c(1)][idx1]
for (i in 1:nrow(sns1)) {
  ratio2[,i] <- barcode_prop[idx1,sns1[i,3]]/ barcode_prop[idx1,sns1[i,1]]
}
order2 = order(rownames(ratio2) %in% barcode_counts$Gene[-c(1)][control1],apply(ratio2,1,median),decreasing=TRUE)

barcode_prop3 = barcode_prop[idx1,]
barcode_test_t1 = rep(0,nrow(barcode_prop3))
for(i in 1:nrow(barcode_prop3)){
  barcode_test_t1[i]= t.test(log(barcode_prop3[i,sns1[,2]]+1E-8),log(barcode_prop3[i,sns1[,1]]+1E-8),alternative="g",paired=T)$p.value
}


barcode_test_t2 = rep(0,nrow(barcode_prop3))
for(i in 1:nrow(barcode_prop3)){
  barcode_test_t2[i]= t.test(log(barcode_prop3[i,sns1[,3]]+1E-8),log(barcode_prop3[i,sns1[,1]]+1E-8),alternative="g",paired=T)$p.value
}
names(barcode_test_t1)=names(barcode_test_t2)=rownames(barcode_prop3)

pval_df = data.frame(geneName=make.names(rownames(barcode_prop3),unique = T),barcode_test_t1=barcode_test_t1,barcode_test_t2=barcode_test_t2,stringsAsFactors = F)
pval_df$sign1 = pval_df$barcode_test_t1<0.05
pval_df$sign2 = pval_df$barcode_test_t2<0.05
pval_df$ypos1=apply(ratio1,1,function(x) max(x[is.finite(x)]))
pval_df$ypos2=apply(ratio2,1,function(x) max(x[is.finite(x)]))





toPlot2 = reshape2::melt(data.frame(t(ratio1[order1[1:50],])))
p_colon_1 = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,10))+geom_text(data=pval_df[which(pval_df$geneName %in% toPlot2$variable),],aes(x=geneName,y=10,label=ifelse(sign1,signif(barcode_test_t1,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Human Colon APC null T0 vs T1")

p_colon_1
ggsave("colon_v5_APC_T0_1.pdf",width = 8.5,height = 4)


toPlot2 = reshape2::melt(data.frame(t(ratio2[order2[1:50],])))
p_colon_2 = ggpubr::ggboxplot(toPlot2,x="variable",y="value")+
  coord_cartesian(ylim = c(0,10))+geom_text(data=pval_df[which(pval_df$geneName %in% toPlot2$variable),],aes(x=geneName,y=ypos2+3,label=ifelse(sign2,signif(barcode_test_t2,2),"")))+
  theme(axis.text.x = element_text(angle = 55, hjust = 1))+ylab("Ratio")+xlab("Gene")+
  ggtitle("Human Colon APC null T0 vs T2")

p_colon_2
ggsave("colon_v5_APC_T0_2.pdf",width = 8.5,height = 4)


