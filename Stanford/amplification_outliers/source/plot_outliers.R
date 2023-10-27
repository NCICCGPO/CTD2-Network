# This code replicate genome wide outlier plots from
#  "Funtional screenings of amplification outlier oncogenes in organoid models of early tumorigenesis"
# Salahudeen et al.
# Fig 1
# joseaseoane@vhio.net


theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}


library(GenomicRanges)
# PAAD
ourliers = read.table("PAAD/outliers_.RData_all.tsv",sep="\t",header=T)

load("PAAD/CNA_samples_annotationTCGA_PAAD.RData")
ourliers = cbind(ourliers,r2[match(ourliers$EntrezId,r2$hgnc_symbol),4:6])

ourliers = ourliers[with(ourliers,order(chrom,start)),]
ourliers$index = 1:nrow(ourliers)
library(plyr)
kk2 = dlply(ourliers,.(chrom),function(miniOut) smooth.spline(x=1:nrow(miniOut),y=miniOut$AmpUpCount,spar=0.5 )$y )
smoothedChrom = do.call(c,kk2)
ourliers$AmpUpCountS = ifelse(smoothedChrom>0,2*smoothedChrom,0)

gr.genes.all = with(ourliers,GRanges(paste("chr",chrom,sep=""),IRanges(start,end),score=AmpUpCount,scorei=AmpUpCountS))
gr.genes.all = keepSeqlevels(gr.genes.all,paste("chr",1:22,sep=""),pruning.mode = "coarse")
gr.genes.all = sortSeqlevels(gr.genes.all)


p.paad = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = score),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="PAAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))
p.paad.i = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = scorei),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="PAAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))


# COAD
ourliers = read.table("COAD/outliers_.RData_all.tsv",sep="\t",header=T)

load("COAD/CNA_samples_annotationTCGA_COAD.RData")
ourliers = cbind(ourliers,r2[match(ourliers$EntrezId,r2$hgnc_symbol),4:6])

ourliers = ourliers[with(ourliers,order(chrom,start)),]
ourliers$index = 1:nrow(ourliers)
library(plyr)
kk2 = dlply(ourliers,.(chrom),function(miniOut) smooth.spline(x=1:nrow(miniOut),y=miniOut$AmpUpCount,spar=0.5 )$y )
smoothedChrom = do.call(c,kk2)
ourliers$AmpUpCountS = ifelse(smoothedChrom>0,2.5*smoothedChrom,0)

gr.genes.all = with(ourliers,GRanges(paste("chr",chrom,sep=""),IRanges(start,end),score=AmpUpCount,scorei=AmpUpCountS))
gr.genes.all = keepSeqlevels(gr.genes.all,paste("chr",1:22,sep=""),pruning.mode = "coarse")
gr.genes.all = sortSeqlevels(gr.genes.all)

p.coad = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = score),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="COAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))
p.coad.i = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = scorei),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="COAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))

# LUAD
ourliers = read.table("LUAD/outliers_.RData_all.tsv",sep="\t",header=T)

load("LUAD/CNA_samples_annotationTCGA_LUAD.RData")
ourliers = cbind(ourliers,r2[match(ourliers$EntrezId,r2$hgnc_symbol),4:6])

ourliers = ourliers[with(ourliers,order(chrom,start)),]
ourliers$index = 1:nrow(ourliers)
library(plyr)
kk2 = dlply(ourliers,.(chrom),function(miniOut) smooth.spline(x=1:nrow(miniOut),y=miniOut$AmpUpCount,spar=0.5 )$y )
smoothedChrom = do.call(c,kk2)
ourliers$AmpUpCountS = ifelse(smoothedChrom>0,2*smoothedChrom,0)

gr.genes.all = with(ourliers,GRanges(paste("chr",chrom,sep=""),IRanges(start,end),score=AmpUpCount,scorei=AmpUpCountS))
gr.genes.all = keepSeqlevels(gr.genes.all,paste("chr",1:22,sep=""),pruning.mode = "coarse")
gr.genes.all = sortSeqlevels(gr.genes.all)

p.luad = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = score),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="LUAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))

p.luad.i = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = scorei),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="LUAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))


# ESCA_SC
ourliers = read.table("ESCA_SC/outliers_.RData_all.tsv",sep="\t",header=T)

load("ESCA_SC/CNA_samples_annotationTCGA_ESCA.RData")
ourliers = cbind(ourliers,r2[match(ourliers$EntrezId,r2$hgnc_symbol),4:6])

ourliers = ourliers[with(ourliers,order(chrom,start)),]
ourliers$index = 1:nrow(ourliers)
library(plyr)
kk2 = dlply(ourliers,.(chrom),function(miniOut) smooth.spline(x=1:nrow(miniOut),y=miniOut$AmpUpCount,spar=0.5 )$y )
smoothedChrom = do.call(c,kk2)
ourliers$AmpUpCountS = ifelse(smoothedChrom>0,2*smoothedChrom,0)

gr.genes.all = with(ourliers,GRanges(paste("chr",chrom,sep=""),IRanges(start,end),score=AmpUpCount,scorei=AmpUpCountS))
gr.genes.all = keepSeqlevels(gr.genes.all,paste("chr",1:22,sep=""),pruning.mode = "coarse")
gr.genes.all = sortSeqlevels(gr.genes.all)

p.escasc = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = score),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="ESCA_SQ")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))
p.escasc.i = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = scorei),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="ESCA_SQ")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))

# STAD
ourliers = read.table("STAD/outliers_.RData_all.tsv",sep="\t",header=T)

load("STAD/CNA_samples_annotationTCGA_STAD.RData")
ourliers = cbind(ourliers,r2[match(ourliers$EntrezId,r2$hgnc_symbol),4:6])

ourliers = ourliers[with(ourliers,order(chrom,start)),]
ourliers$index = 1:nrow(ourliers)
library(plyr)
kk2 = dlply(ourliers,.(chrom),function(miniOut) smooth.spline(x=1:nrow(miniOut),y=miniOut$AmpUpCount)$y )
smoothedChrom = do.call(c,kk2)
ourliers$AmpUpCountS = ifelse(smoothedChrom>0,2*smoothedChrom,0)

gr.genes.all = with(ourliers,GRanges(paste("chr",chrom,sep=""),IRanges(start,end),score=AmpUpCount,scorei=AmpUpCountS))
gr.genes.all = keepSeqlevels(gr.genes.all,paste("chr",1:22,sep=""),pruning.mode = "coarse")
gr.genes.all = sortSeqlevels(gr.genes.all)

p.stad = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = score),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="STAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))
p.stad.i = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = scorei),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="STAD")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))


# HNSC
ourliers = read.table("HNSC/outliers_.RData_all.tsv",sep="\t",header=T)

load("HNSC/CNA_samples_annotationTCGA_HNSC.RData")
ourliers = cbind(ourliers,r2[match(ourliers$EntrezId,r2$hgnc_symbol),4:6])

ourliers = ourliers[with(ourliers,order(chrom,start)),]
ourliers$index = 1:nrow(ourliers)
library(plyr)
kk2 = dlply(ourliers,.(chrom),function(miniOut) smooth.spline(x=1:nrow(miniOut),y=miniOut$AmpUpCount)$y )
smoothedChrom = do.call(c,kk2)
ourliers$AmpUpCountS = ifelse(smoothedChrom>0,2*smoothedChrom,0)

gr.genes.all = with(ourliers,GRanges(paste("chr",chrom,sep=""),IRanges(start,end),score=AmpUpCount,scorei=AmpUpCountS))
gr.genes.all = keepSeqlevels(gr.genes.all,paste("chr",1:22,sep=""),pruning.mode = "coarse")
gr.genes.all = sortSeqlevels(gr.genes.all)

p.hnsc = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = score),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="HNSC")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))
p.hnsc.i = plotGrandLinear(gr.genes.all, coord = "genome", geom = "bar", aes(y = scorei),color="#386cb0", space.skip = 0.005,legend=F,spaceline=T,main="HNSC")+
  theme_Publication(base_size = 9,base_family = "Helvetica")+
  theme(legend.position = "none",axis.text.x = element_text(angle = 55, hjust = 1))


pdf("p_smoothed.pdf",width = 4.5,height = 3.5)
p.coad.i
p.escasc.i
p.luad.i
p.paad.i
p.stad.i
dev.off()

pdf("p_original.pdf",width = 4.5, height = 3.5)
p.coad
p.escasc
p.luad
p.paad
p.stad
dev.off()

pdf("p_smoothed_v2.pdf",width = 9,height = 7)
p.coad.i
p.escasc.i
p.luad.i
p.paad.i
p.stad.i
p.hnsc.i
dev.off()

pdf("p_original_v2.pdf",width = 9, height = 7)
p.coad
p.escasc
p.luad
p.paad
p.stad
p.hnsc
dev.off()

pdf("p_smoothed_b.pdf",width = 10,height = 8)
p.coad.i
p.escasc.i
p.luad.i
p.paad.i
p.stad.i
dev.off()

pdf("p_original_b.pdf",width = 10, height = 8)
p.coad
p.escasc
p.luad
p.paad
p.stad
dev.off()
