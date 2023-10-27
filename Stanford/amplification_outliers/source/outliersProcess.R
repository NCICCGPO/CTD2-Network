# helper functions for outlier plots from
#  "Funtional screenings of amplification outlier oncogenes in organoid models of early tumorigenesis"
# Salahudeen et al.
# joseaseoane@vhio.net


getOutliersmRNA = function(CNAFile,exp_met,Des,samples,outputName,cores=1,r2){
  
  exp_met = apply(exp_met,1,function(x) log2(x+1))
  exp_met = exp_met[rownames(exp_met) %in% samples,]
  
  load(CNAFile)
  samplesAmp = samplesAmp[samples,]
  samplesGAIN = samplesGAIN[samples,]
  samplesHETD = samplesHETD[samples,]
  samplesHOMD = samplesHOMD[samples,]
  
  
  #Des = as.character(Des)
  
  
  library(MASS)
  library(doMC)
  registerDoMC(cores=cores)
  
  #geneNames = colnames(samplesAmp)
  #geneNames = unique(Des[,1])
  geneNames = unique(as.character(Des))
  
  tempOut  = NULL
  outliersCount = foreach(i =1:length(geneNames),.combine = rbind,.errorhandling="pass")%dopar%{
    #outliersCount = foreach(i =1:100,.combine = rbind)%dopar%{
    
    geneId = geneNames[i]
    #geneName = as.character(annotationMet$Gene_symbol[which(annotationMet$Entrez_Gene_ID_0==geneId)][1])
    
    #idProbes = as.character(annotationMet$Probe_id[which(annotationMet$Entrez_Gene_ID_0==geneId)])
    if(as.character(geneId) %in% colnames(samplesAmp)){
      
    
    est = fitdistr(na.omit(exp_met[,which(Des==geneId)]),"normal")$estimate
    
    qvalues = qnorm(c(0.05,0.95),est[1],est[2])
    leftTail = which(exp_met[,which(Des==geneId)]<qvalues[1])
    rigthTail = which(exp_met[,which(Des==geneId)]>qvalues[2])
    
    leftSamples = rownames(exp_met)[ leftTail]
    rigthSamples = rownames(exp_met)[ rigthTail]
    
    AmpUpCount = sum(samplesAmp[ rownames(samplesAmp) %in% rigthSamples,as.character(geneId)])
    AmpDownCount = sum(samplesAmp[rownames(samplesAmp) %in% leftSamples,as.character(geneId)])
    DelUpCount = sum(samplesHOMD[rownames(samplesHOMD) %in% rigthSamples,as.character(geneId)])
    DelDownCount = sum(samplesHOMD[rownames(samplesHOMD) %in% leftSamples,as.character(geneId)])
    AmpGainUpCount = sum(samplesAmp[rownames(samplesAmp) %in% rigthSamples,as.character(geneId)] | samplesGAIN[rownames(samplesGAIN) %in% rigthSamples,as.character(geneId)])
    HetHomDDownCount = sum(samplesHOMD[rownames(samplesHOMD) %in% leftSamples,as.character(geneId)] | samplesHETD[rownames(samplesHETD) %in% leftSamples,as.character(geneId)])
    totalAmp = sum(samplesAmp[,as.character(geneId)])
    totalGain= sum(samplesGAIN[,as.character(geneId)])
    totalNeut= sum(samplesNEUT[,as.character(geneId)])
    totalHomD= sum(samplesHOMD[,as.character(geneId)])
    totalHetD= sum(samplesHETD[,as.character(geneId)])
    gc()
    
    tempOut =rbind(tempOut,c(geneId,AmpUpCount,AmpDownCount,DelUpCount,DelDownCount,AmpGainUpCount,HetHomDDownCount,totalAmp,totalGain,totalNeut,totalHetD,totalHomD,length(rigthSamples),length(leftSamples)))
    tempOut
    }
  }
  
  colnames(outliersCount)= c("EntrezId","AmpUpCount","AmpDownCount","DelUpCount","DelDownCount","AmpGainUpCount","HetHomDDownCount","totalAmp","totalGain","totalNeut","totalHetD","totalHomD","ExpOutRight","ExpOutLeft")

  write.table(outliersCount,file=paste(outputName,"_all.tsv",sep=""),quote = F,sep = "\t",row.names = F,col.names = T)
  #outliersCount = read.table("~/projects/outliers_pancan/PAAD/outliers_.RData_all.tsv",header = T,sep="\t",stringsAsFactors = F)
  idch = match(outliersCount[,1],r2$hgnc_symbol)
  
  AmpOutliers = NULL
  try({
  AmpOutliers = data.frame(Genename=r2$hgnc_symbol[idch],Chr=r2$chrom[idch],Start=r2$start[idch],End=r2$end[idch],Rank=NA,Outlier=as.numeric(outliersCount[,2]),Amp=outliersCount[,8],Gain=outliersCount[,9],ExpCount = outliersCount[,13] )
  AmpOutliers = AmpOutliers[ AmpOutliers$Outlier>0,]
  idord = order(AmpOutliers$Outlier,decreasing = T)
  AmpOutliers = AmpOutliers[idord,]
  AmpOutliers$Rank=1:nrow(AmpOutliers)
  })
  
  try({
  AmpGainOutliers = data.frame(Genename=r2$hgnc_symbol[idch],Chr=r2$chrom[idch],Start=r2$start[idch],End=r2$end[idch],Rank=NA,Outlier=as.numeric(outliersCount[,6]),Amp=outliersCount[,8],Gain=outliersCount[,9])
  AmpGainOutliers = AmpGainOutliers[ AmpGainOutliers$Outlier>0,]
  idord = order(AmpGainOutliers$Outlier,decreasing = T)
  AmpGainOutliers = AmpGainOutliers[idord,]
  AmpGainOutliers$Rank=1:nrow(AmpGainOutliers)
  })
  
  HomDOutliers = NULL
  try({
  HomDOutliers = data.frame(Genename=r2$hgnc_symbol[idch],Chr=r2$chrom[idch],Start=r2$start[idch],End=r2$end[idch],Rank=NA,Outlier=as.numeric(outliersCount[,5]),HetD=outliersCount[,11],HomD=outliersCount[,12],ExpCount = outliersCount[,14] )
  HomDOutliers = HomDOutliers[ HomDOutliers$Outlier>0,]
  idord = order(HomDOutliers$Outlier,decreasing = T)
  HomDOutliers = HomDOutliers[idord,]
  HomDOutliers$Rank=1:nrow(HomDOutliers)
  })
  
  try({
  HomHetDOutliers = data.frame(Genename=r2$hgnc_symbol[idch],Chr=r2$chrom[idch],Start=r2$start[idch],End=r2$end[idch],Rank=NA,Outlier=as.numeric(outliersCount[,7]),HetD=outliersCount[,11],HomD=outliersCount[,12])
  HomHetDOutliers = HomHetDOutliers[ HomHetDOutliers$Outlier>0,]
  idord = order(HomHetDOutliers$Outlier,decreasing = T)
  HomHetDOutliers = HomHetDOutliers[idord,]
  HomHetDOutliers$Rank=1:nrow(HomHetDOutliers)
  })
  
  write.table(AmpOutliers,file=paste(outputName,"_amp.tsv",sep=""),quote = F,sep = "\t",row.names = F)
  write.table(AmpGainOutliers,file=paste(outputName,"_ampgain.tsv",sep=""),quote = F,sep = "\t",row.names = F)
  write.table(HomDOutliers,file=paste(outputName,"_homd.tsv",sep=""),quote = F,sep = "\t",row.names = F)
  write.table(HomHetDOutliers,file=paste(outputName,"_homhetD.tsv",sep=""),quote = F,sep = "\t",row.names = F)
  
  library(ggbio)
  library(GenomicRanges)
  library(biovizBase)
  
    
  upr = AmpOutliers
  dwn = HomDOutliers
  
  gene.upr = with(upr,GRanges(paste("chr",Chr,sep=""),IRanges(Start,End),score=Outlier))
  gene.dwn = with(dwn,GRanges(paste("chr",Chr,sep=""),IRanges(Start,End),score=Outlier))
  
  
  gr.genes.upr = gene.upr #keepSeqlevels(gene.upr,c(1:22))
  
  gr.genes.dwn = gene.dwn #keepSeqlevels(gene.dwn,c(1:22))
  
  gr.genes.all = gr.genes.upr
  gr.genes.all$var = "Upregulated"
  aux = gr.genes.dwn
  aux$var ="Downregulated"
  aux$score = -aux$score
  gr.genes.all = c(gr.genes.all,aux)
  
  gr.genes.all = sortSeqlevels(gr.genes.all)
  
  
  #pdf(file=paste(outputName,"_Plot.pdf",sep=""),width = 15)
  #autoplot(gr.genes.all, coord = "genome", geom = "bar", aes(y = score,col=var), space.skip = 0.01,legend=F)
  #dev.off()

  plotGrandLinear(sortSeqlevels(gr.genes.all),aes(y = score,col=var),geom="bar")
  ggsave(filename = paste(outputName,"_Plot.pdf",sep=""),width = 15)
  
  
  
}

plotOutliersmRNA = function(outputName,amp,del){
  library(ggbio)
  library(GenomicRanges)
  library(biovizBase)
  
  
  upr = read.table(amp,header = T,sep = "\t",stringsAsFactors = F)
  dwn = read.table(del,header = T,sep = "\t",stringsAsFactors = F)
  
  gene.upr = with(upr,GRanges(paste("chr",Chr,sep=""),IRanges(Start,End),score=Outlier))
  gene.dwn = with(dwn,GRanges(paste("chr",Chr,sep=""),IRanges(Start,End),score=Outlier))
  
  seqLevels = paste("chr",c(1:22,"X","Y"),sep="")
  
  
  gr.genes.upr = gene.upr #keepSeqlevels(gene.upr,c(1:22))
  #seqlevels(gr.genes.upr)=seqLevels#sortSeqlevels(seqLevels)
  gr.genes.dwn = gene.dwn #keepSeqlevels(gene.dwn,c(1:22))
  #seqlevels(gr.genes.dwn)=seqLevels#sortSeqlevels(seqLevels)
  
  gr.genes.all = gr.genes.upr
  gr.genes.all$var = "Upregulated"
  aux = gr.genes.dwn
  aux$var ="Downregulated"
  aux$score = -aux$score
  #gr.genes.all = c(gr.genes.all,aux)
  gr.genes.all.1 = c(gr.genes.all,aux)
 
  
  #gr.genes.all = sortSeqlevels(gr.genes.all)
  #gr.genes.all.srt = sort(gr.genes.all)
  
  #pdf(file=paste(outputName,"_Plot.pdf",sep=""),width = 15)
  #autoplot(gr.genes.all.1, geom = "bar", aes(y = score,col=var), space.skip = 0.01,legend=F)
  #autoplot(gr.genes.all.1, coord = "genome",geom = "bar", aes(y = score,col=var), space.skip = 0.01,legend=F)
    plotGrandLinear(sortSeqlevels(gr.genes.all.1),aes(y = score,col=var),geom="bar")
  ggsave(filename = paste(outputName,"_Plot.pdf",sep=""),width = 15)
  #dev.off()
  
  plotGrandLinear(sort(gene.upr),aes(y=score))
  
}
