# This code replicate ORF screening lorenz plots from 
#  "Funtional screenings of amplification outlier oncogenes in organoid models of early tumorigenesis"
# Salahudeen et al.
# joseaseoane@vhio.net


library(gglorenz)
library(ineq)
library(dplyr)

# eso

load("barcode_counts_filtered_esca_dic.RData")


eso_counts_l = reshape2::melt(barcode_counts_filtered_esca)

eso_counts_l_2 = eso_counts_l[which(!eso_counts_l$variable %in% c("EA","EB","EC","ED")),]
m1 = eso_counts_l_2$value[which(eso_counts_l_2$variable=="A0")]
Gini(m1)
compute_stats= function(x){
  lc1 = Lc(x)
  g1 = Gini(x)
  auc = pracma::trapz(lc1$L,lc1$p)
  list(gini=g1,auc=auc)
}

eso_stats = do.call(rbind.data.frame,eso_counts_l_2%>% 
                      group_by(variable) %>%
                      group_map(~compute_stats(.x$value))) 

eso_stats$sample = unique(eso_counts_l_2$variable)
eso_stats$time = c(rep("T0",4),rep("T46",4))
eso_stats$string = paste(eso_stats$sample," gini:",signif(eso_stats$gini,2)," auc:" ,signif(eso_stats$auc,2),sep="")
eso_counts_l_2$sample = eso_stats$string[match(eso_counts_l_2$variable,eso_stats$sample)]
eso_counts_l_2$time = eso_stats$time[match(eso_counts_l_2$variable,eso_stats$sample)]
eso_counts_l_2$sample_0 = gsub("[^A-D]","",eso_counts_l_2$variable)

p_eso = ggplot(eso_counts_l_2[eso_counts_l_2$variable %in% c("A0","B0","C0","D0"),],aes(value,col=sample))+
  stat_lorenz(desc = T)+
  coord_fixed()+
  geom_abline(linetype="dashed")+
  theme_minimal()+theme(legend.position="bottom")+
  labs(x ="ORF ranked by abundance",
       y= "cumulative fraction of ORF reads",
       title="Esophageal")

gghistogram(eso_counts_l_2,x="value",fill="time",facet.by = "sample_0",bins=20,alpha =0.5  )

p_dis_eso=ggplot(eso_counts_l_2[which(eso_counts_l_2$Control=='No'),],aes(x=reorder(CloneID,value),y=value,col=time))+
  geom_point()+facet_grid(time~sample_0)+
  theme_minimal()+  theme(axis.text.x=element_blank())+labs(x="ORF sorted by total ammount",y="counts",title="Esophageal")


# oral
load("barcode_counts_filtered_oral_dic.RData")

oral_counts_l = reshape2::melt(barcode_counts_filtered_oral)
oral_counts_l_2 = oral_counts_l[which(!oral_counts_l$variable %in% c("AM T0 A", "AM T0 B", "AM T0 C","AM T0 D",
                                                                     "AM T55 A","AM T55 B","AM T55 C","AM T55 D",
                                                                     "KY52 F12 A", "KY52 F12 B", "KY52 F12 C","KY52 F12 D")),]

oral_stats = do.call(rbind.data.frame,oral_counts_l_2 %>% 
                       group_by(variable) %>%
                       group_map(~compute_stats(.x$value))) 

oral_stats$sample = unique(oral_counts_l_2$variable)
oral_stats$time = c(rep("KY0",4),rep("KY52",4))
oral_stats$string = paste(oral_stats$sample," gini:",signif(oral_stats$gini,2)," auc:" ,signif(oral_stats$auc,2),sep="")
oral_counts_l_2$sample = oral_stats$string[match(oral_counts_l_2$variable,oral_stats$sample)]
oral_counts_l_2$time = oral_stats$time[match(oral_counts_l_2$variable,oral_stats$sample)]
oral_counts_l_2$sample_0 = gsub("[^A-D]","",oral_counts_l_2$variable)

p_oral=ggplot(oral_counts_l_2[which(oral_counts_l_2$variable %in% c("KY0 A", "KY0 B", "KY0 C", "KY0 D" )),],aes(value,col=sample))+
  stat_lorenz(desc = T)+
  coord_fixed()+
  geom_abline(linetype="dashed")+
  theme_minimal()+theme(legend.position="bottom")+
  labs(x ="ORF ranked by abundance",
       y= "cumulative fraction of ORF reads",
       title="Oral")

p_dis_oral=ggplot(oral_counts_l_2[which(oral_counts_l_2$Control=='No'),],aes(x=reorder(CloneID,value),y=value,col=time))+
  geom_point()+facet_grid(time~sample_0)+
  theme_minimal()+  theme(axis.text.x=element_blank())+labs(x="ORF sorted by total ammount",y="counts",title="Oral")



#lung

load("barcode_counts_filtered_lung_dic.RData")
lung_counts_l = reshape2::melt(barcode_counts_filtered_lung[-c(1),])

lung_counts_l_2 = lung_counts_l[which(lung_counts_l$variable %in% c("B07","B08","B09","B10","A06","B11","A11","B12")),]

lung_stats = do.call(rbind.data.frame,lung_counts_l_2 %>% 
                       group_by(variable) %>%
                       group_map(~compute_stats(.x$value))) 

lung_stats$sample = unique(lung_counts_l_2$variable)
lung_stats$time = c("T0","T0","T0","TP4","T0","TP4","TP4","TP4")
lung_stats$string = paste(lung_stats$sample," gini:",signif(lung_stats$gini,2)," auc:" ,signif(lung_stats$auc,2),sep="")
lung_stats$sample_0 = c("A06","A11","B07","A06","B09","A11","B07","B09")
lung_counts_l_2$sample = lung_stats$string[match(lung_counts_l_2$variable,lung_stats$sample)]
lung_counts_l_2$time = lung_stats$time[match(lung_counts_l_2$variable,lung_stats$sample)]
lung_counts_l_2$sample_0 = lung_stats$sample_0[match(lung_counts_l_2$variable,lung_stats$sample)]



p_lung=ggplot(lung_counts_l_2[lung_counts_l_2$variable %in% c("A06","A11","B07","B08"),],aes(value,col=sample))+
  stat_lorenz(desc = T)+
  coord_fixed()+
  geom_abline(linetype="dashed")+
  theme_minimal()+theme(legend.position="bottom")+
  labs(x ="ORF ranked by abundance",
       y= "cumulative fraction of ORF reads",
       title="Lung")

p_dis_lung=ggplot(lung_counts_l_2[which(lung_counts_l_2$Control=='No'),],aes(x=reorder(CloneID,value),y=value,col=time))+
  geom_point()+facet_grid(time~sample_0)+
  theme_minimal()+  theme(axis.text.x=element_blank())+labs(x="ORF sorted by total ammount",y="counts",title="Lung")


#pancreas
load("barcode_counts_filtered_pancreas_dic.RData")
pancreas_counts_l = reshape2::melt(barcode_counts_filtered_pancreas)

pancreas_counts_l_2 = pancreas_counts_l[which(!pancreas_counts_l$variable =="mPanKBO4_D0"),]

pancreas_stats = do.call(rbind.data.frame,pancreas_counts_l_2 %>% 
                           group_by(variable) %>%
                           group_map(~compute_stats(.x$value))) 

pancreas_stats$sample = unique(pancreas_counts_l_2$variable)
pancreas_stats$time = c(rep("T33",3),rep("T0",3))
pancreas_stats$sample_0 = c("01","02","03","01","02","03")
pancreas_stats$string = paste(pancreas_stats$sample," gini:",signif(pancreas_stats$gini,2)," auc:" ,signif(pancreas_stats$auc,2),sep="")
pancreas_counts_l_2$sample = pancreas_stats$string[match(pancreas_counts_l_2$variable,pancreas_stats$sample)]
pancreas_counts_l_2$time = pancreas_stats$time[match(pancreas_counts_l_2$variable,pancreas_stats$sample)]
pancreas_counts_l_2$sample_0 = pancreas_stats$sample_0[match(pancreas_counts_l_2$variable,pancreas_stats$sample)]

p_pancreas=ggplot(pancreas_counts_l_2[which(pancreas_counts_l_2$variable %in% c("mPanKBO1_D0","mPanKBO2_D0","mPanKBO3_D0")),],aes(value,col=sample))+
  stat_lorenz(desc = T)+
  coord_fixed()+
  geom_abline(linetype="dashed")+
  theme_minimal()+theme(legend.position="bottom")+
  labs(x ="ORF ranked by abundance",
       y= "cumulative fraction of ORF reads",
       title="Pancreas")

p_dis_paad=ggplot(pancreas_counts_l_2[which(pancreas_counts_l_2$Control=='No'),],aes(x=reorder(CloneID,value),y=value,col=time))+
  geom_point()+facet_grid(time~sample_0)+
  theme_minimal()+  theme(axis.text.x=element_blank())+labs(x="ORF sorted by total ammount",y="counts",title="Pancreas")


#gastric
load("barcode_counts_filtered_gastric_dic.RData")
gastric_counts_l = reshape2::melt(barcode_counts_filtered_gastric[-c(1),])

gastric_counts_l_2 = gastric_counts_l[which(gastric_counts_l$variable %in% colnames(barcode_counts_filtered_gastric)[5:12] ),]

gastric_stats = do.call(rbind.data.frame,gastric_counts_l_2 %>% 
                          group_by(variable) %>%
                          group_map(~compute_stats(.x$value))) 

gastric_stats$sample = unique(gastric_counts_l_2$variable)
gastric_stats$time = c(rep("T0",4),rep("T55",4))
gastric_stats$sample_0 = rep(c("A","B","C","D"),2)
gastric_stats$string = paste(gastric_stats$sample," gini:",signif(gastric_stats$gini,2)," auc:" ,signif(gastric_stats$auc,2),sep="")
gastric_counts_l_2$sample = gastric_stats$string[match(gastric_counts_l_2$variable,gastric_stats$sample)]
gastric_counts_l_2$time = gastric_stats$time[match(gastric_counts_l_2$variable,gastric_stats$sample)]
gastric_counts_l_2$sample_0 = gastric_stats$sample_0[match(gastric_counts_l_2$variable,gastric_stats$sample)]


p_gastric=ggplot(gastric_counts_l_2[which(gastric_counts_l_2$variable %in% c("AM.T0.A","AM.T0.B","AM.T0.C","AM.T0.D")),],aes(value,col=sample))+
  stat_lorenz(desc = T)+
  coord_fixed()+
  geom_abline(linetype="dashed")+
  theme_minimal()+theme(legend.position="bottom")+
  labs(x ="ORF ranked by abundance",
       y= "cumulative fraction of ORF reads",
       title="Gastric")


p_dis_gastric=ggplot(gastric_counts_l_2[which(gastric_counts_l_2$Control=='No'),],aes(x=reorder(CloneID,value),y=value,col=time))+
  geom_point()+facet_grid(time~sample_0)+
  theme_minimal()+  theme(axis.text.x=element_blank())+labs(x="ORF sorted by total ammount",y="counts",title="Gastric")



#colon
load("barcode_counts_filtered_colon_dic.RData")
colon_counts_l = reshape2::melt(barcode_counts_filtered_colon[-c(1),])

colon_counts_l_2 = colon_counts_l[which(colon_counts_l$variable %in% c("A09","B03","A10","B06")),]

colon_stats = do.call(rbind.data.frame,colon_counts_l_2 %>% 
                        group_by(variable) %>%
                        group_map(~compute_stats(.x$value))) 

colon_stats$sample = unique(colon_counts_l_2$variable)
colon_stats$time = c(rep("T0",2),rep("T2",2))
colon_stats$string = paste(colon_stats$sample," gini:",signif(colon_stats$gini,2)," auc:" ,signif(colon_stats$auc,2),sep="")
colon_stats$sample_0 = c("A09","A10","A09","A10")
colon_counts_l_2$sample = colon_stats$string[match(colon_counts_l_2$variable,colon_stats$sample)]
colon_counts_l_2$time = colon_stats$time[match(colon_counts_l_2$variable,colon_stats$sample)]
colon_counts_l_2$sample_0 = colon_stats$sample_0[match(colon_counts_l_2$variable,colon_stats$sample)]


p_colon=ggplot(colon_counts_l_2[which(colon_counts_l_2$variable %in% c("A09","A10")),],aes(value,col=sample))+
  stat_lorenz(desc = T)+
  coord_fixed()+
  geom_abline(linetype="dashed")+
  theme_minimal()+ theme(legend.position="bottom")+
  labs(x ="ORF ranked by abundance",
       y= "cumulative fraction of ORF reads",
       title="Colon")

p_dis_colon=ggplot(colon_counts_l_2[which(colon_counts_l_2$Control=='No'),],aes(x=reorder(CloneID,value),y=value,col=time))+
  geom_point()+facet_grid(time~sample_0)+
  theme_minimal()+  theme(axis.text.x=element_blank())+labs(x="ORF sorted by total ammount",y="counts",title="Colon")


library(ggpubr)
ggarrange(p_eso,p_oral,p_lung,p_pancreas,p_gastric,p_colon,ncol = 2,nrow=3)
ggsave("lorenz_plots_thresholded_v2.pdf",width = 20,height = 18)

ggarrange(p_dis_eso,p_dis_oral,p_dis_lung,p_dis_paad,p_dis_gastric,p_dis_colon,ncol = 2,nrow=3)
ggsave("dist_plots_thresholded_v2.pdf",width = 20,height = 18)
