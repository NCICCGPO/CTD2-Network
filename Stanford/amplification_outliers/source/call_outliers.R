# This code replicate outlier analysis from
#  "Funtional screenings of amplification outlier oncogenes in organoid models of early tumorigenesis"
# Salahudeen et al.
# Fig 1
# joseaseoane@vhio.net


# call outliers

# Download TCGA expression data and CN profiles
# LUAD
load("LUAD__illuminahiseq_rnaseqv2__rsem.genes.normalized.rda")
## data is a matrix with row genes and colum samples, Des is a dictionary for different gene identifiers
types = substr(colnames(Data),start = 14,stop = 15)
mRNA.Data = Data[,which(types == "01")]
colnames(mRNA.Data)=substr(colnames(mRNA.Data),1,12)
mRNA.Des = Des

# first get the thresholds from raw data (not gene anotated), and the n
CNAraw = read.table("LUAD__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Feb-12-2015.txt",sep="\t",header = T) 
type = substr(CNAraw$Sample,start = 14,stop = 15)
CNAraw = CNAraw[ type=="01",]
CNAraw$Sample = substr(CNAraw$Sample,1,12) 

CNAsamples = unique(CNAraw$Sample)
CNAraw$call = NA
for(i in 1:length(CNAsamples)){
  sample1 = CNAraw[ CNAraw$Sample== as.character(CNAsamples[i]),6]
  qunt <- quantile(sample1, c(0.075, 0.125, 0.15, 0.25, 0.275, 0.315, 0.35, 0.375, 0.5, 0.625, 0.65, 0.685, 0.725, 0.75, 0.85, 0.875, 0.925), na.rm=TRUE)
  median = qunt[9]
  med = median(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  all.sds0.5 = sd(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  ampsc =1* (sample1>(med+2*all.sds0.5)) + 1*  (sample1>(med+6*all.sds0.5)) -1* (sample1<(med-2.5*all.sds0.5)) -1 * (sample1<(med-5*all.sds0.5))
  CNAraw$call[CNAraw$Sample== as.character(CNAsamples[i])]=ampsc
}

library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])


CNAamp = CNAraw[CNAraw$call=="2",]
CNAneu = CNAraw[CNAraw$call=="0",]
CNAgain = CNAraw[CNAraw$call=="1",]
CNAhomd = CNAraw[CNAraw$call=="-2",]
CNAhetd = CNAraw[CNAraw$call=="-1",]

nsamples = length(CNAsamples)
ngenes = length(geneNames)
samplesAmp = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesAmp)=geneNames
rownames(samplesAmp)=CNAsamples


samplesGAIN = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesGAIN)=geneNames
rownames(samplesGAIN)=CNAsamples

samplesHETD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHETD)=geneNames
rownames(samplesHETD)=CNAsamples

samplesHOMD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHOMD)=geneNames
rownames(samplesHOMD)=CNAsamples

samplesNEUT = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesNEUT)=geneNames
rownames(samplesNEUT)=CNAsamples

samplesNUM = matrix(nrow = nsamples,ncol = ngenes,data=NA)
colnames(samplesNUM)=geneNames
rownames(samplesNUM)=CNAsamples

for (i in 1:length(geneNames)){
  idr = which(r2$hgnc_symbol==geneNames[i])
  chr = as.character(r2$chrom[idr])
  ini = as.numeric(r2$start[idr])
  end = as.numeric(r2$end[idr])
  
  idamp = which(CNAamp$Chromosome==chr&( (ini>=CNAamp$Start&end<=CNAamp$End)| #
                                           (ini<=CNAamp$End&ini>=CNAamp$Start)| # 
                                           (ini<=CNAamp$Start& end>=CNAamp$End )| # 
                                           (end>= CNAamp$Start & end<=CNAamp$End ) ) ) 
  
  idamp = which(CNAamp$Chromosome==chr&( (end>=CNAamp$Start&end<=CNAamp$End)|(ini<=CNAamp$End&ini>=CNAamp$Start) ) )  
  idgain = which(CNAgain$Chromosome==chr&( (ini>=CNAgain$Start&end<=CNAgain$End)| #
                                             (ini<=CNAgain$End&ini>=CNAamp$Start)| # 
                                             (ini<=CNAgain$Start& end>=CNAgain$End )| # 
                                             (end>= CNAgain$Start & end<=CNAgain$End ) ) ) 
  idneu = which(CNAneu$Chromosome==chr&( (ini>=CNAneu$Start&end<=CNAneu$End)| #
                                           (ini<=CNAneu$End&ini>=CNAneu$Start)| # 
                                           (ini<=CNAneu$Start& end>=CNAneu$End )| # 
                                           (end>= CNAneu$Start & end<=CNAneu$End ) ) ) 
  idhomd = which(CNAhomd$Chromosome==chr&( (ini>=CNAhomd$Start&end<=CNAhomd$End)| #
                                             (ini<=CNAhomd$End&ini>=CNAhomd$Start)| # 
                                             (ini<=CNAhomd$Start& end>=CNAhomd$End )| # 
                                             (end>= CNAhomd$Start & end<=CNAhomd$End ) ) ) 
  idhetd = which(CNAhetd$Chromosome==chr&( (ini>=CNAhetd$Start&end<=CNAhetd$End)| #
                                             (ini<=CNAhetd$End&ini>=CNAhetd$Start)| # 
                                             (ini<=CNAhetd$Start& end>=CNAhetd$End )| # 
                                             (end>= CNAhetd$Start & end<=CNAhetd$End ) ) ) 
  
  
  
  samplesNEUT[as.character(CNAneu$Sample[idneu]),i]=T
  samplesAmp[as.character(CNAamp$Sample[idamp]),i]=T
  samplesGAIN[as.character(CNAgain$Sample[idgain]),i]=T
  samplesHETD[as.character(CNAhetd$Sample[idhetd]),i]=T
  samplesHOMD[as.character(CNAhomd$Sample[idhomd]),i]=T
  
}

save(samplesNEUT,samplesAmp,samplesGAIN,samplesHETD,samplesHOMD,samplesNUM,r2,file="LUAD/CNA_samples_annotationTCGA_LUAD.RData")


load("LUAD/CNA_samples_annotationTCGA_LUAD.RData")


library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])

source("outliersProcess.R")

allSamples = rownames(samplesAmp)
getOutliersmRNA(CNAFile = "LUAD/CNA_samples_annotationTCGA_LUAD.RData",exp_met = mRNA.Data,Des = mRNA.Des[,1],outputName = "LUAD/outliers_.RData",samples = allSamples,cores = 10,r2 = r2)


# ESCA_SC
load("ESCA__illuminahiseq_rnaseqv2__gene.quantification.rda")

types = substr(colnames(Data),start = 14,stop = 15)
mRNA.Data = Data[,which(types == "01")]
colnames(mRNA.Data)=substr(colnames(mRNA.Data),1,12)
mRNA.Des = Des

# first get the thresholds from raw data (not gene anotated), and the n
CNAraw = read.table("ESCA__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Feb-12-2015.txt",sep="\t",header = T) 
type = substr(CNAraw$Sample,start = 14,stop = 15)
CNAraw = CNAraw[ type=="01",]
CNAraw$Sample = substr(CNAraw$Sample,1,12) 

CNAsamples = unique(CNAraw$Sample)
CNAraw$call = NA
for(i in 1:length(CNAsamples)){
  sample1 = CNAraw[ CNAraw$Sample== as.character(CNAsamples[i]),6]
  qunt <- quantile(sample1, c(0.075, 0.125, 0.15, 0.25, 0.275, 0.315, 0.35, 0.375, 0.5, 0.625, 0.65, 0.685, 0.725, 0.75, 0.85, 0.875, 0.925), na.rm=TRUE)
  median = qunt[9]
  med = median(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  all.sds0.5 = sd(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  ampsc =1* (sample1>(med+2*all.sds0.5)) + 1*  (sample1>(med+6*all.sds0.5)) -1* (sample1<(med-2.5*all.sds0.5)) -1 * (sample1<(med-5*all.sds0.5))
  CNAraw$call[CNAraw$Sample== as.character(CNAsamples[i])]=ampsc
}

library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])


CNAamp = CNAraw[CNAraw$call=="2",]
CNAneu = CNAraw[CNAraw$call=="0",]
CNAgain = CNAraw[CNAraw$call=="1",]
CNAhomd = CNAraw[CNAraw$call=="-2",]
CNAhetd = CNAraw[CNAraw$call=="-1",]

nsamples = length(CNAsamples)
ngenes = length(geneNames)
samplesAmp = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesAmp)=geneNames
rownames(samplesAmp)=CNAsamples


samplesGAIN = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesGAIN)=geneNames
rownames(samplesGAIN)=CNAsamples

samplesHETD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHETD)=geneNames
rownames(samplesHETD)=CNAsamples

samplesHOMD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHOMD)=geneNames
rownames(samplesHOMD)=CNAsamples

samplesNEUT = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesNEUT)=geneNames
rownames(samplesNEUT)=CNAsamples

samplesNUM = matrix(nrow = nsamples,ncol = ngenes,data=NA)
colnames(samplesNUM)=geneNames
rownames(samplesNUM)=CNAsamples

for (i in 1:length(geneNames)){
  idr = which(r2$hgnc_symbol==geneNames[i])
  chr = as.character(r2$chrom[idr])
  ini = as.numeric(r2$start[idr])
  end = as.numeric(r2$end[idr])
  
 
  
  idamp = which(CNAamp$Chromosome==chr&( (ini>=CNAamp$Start&end<=CNAamp$End)| #
                                           (ini<=CNAamp$End&ini>=CNAamp$Start)| # 
                                           (ini<=CNAamp$Start& end>=CNAamp$End )| # 
                                           (end>= CNAamp$Start & end<=CNAamp$End ) ) ) 
  
  #idamp = which(CNAamp$Chromosome==chr&( (end>=CNAamp$Start&end<=CNAamp$End)|(ini<=CNAamp$End&ini>=CNAamp$Start) ) )  
  idgain = which(CNAgain$Chromosome==chr&( (ini>=CNAgain$Start&end<=CNAgain$End)| #
                                             (ini<=CNAgain$End&ini>=CNAamp$Start)| # 
                                             (ini<=CNAgain$Start& end>=CNAgain$End )| # 
                                             (end>= CNAgain$Start & end<=CNAgain$End ) ) ) 
  idneu = which(CNAneu$Chromosome==chr&( (ini>=CNAneu$Start&end<=CNAneu$End)| #
                                           (ini<=CNAneu$End&ini>=CNAneu$Start)| # 
                                           (ini<=CNAneu$Start& end>=CNAneu$End )| # 
                                           (end>= CNAneu$Start & end<=CNAneu$End ) ) ) 
  idhomd = which(CNAhomd$Chromosome==chr&( (ini>=CNAhomd$Start&end<=CNAhomd$End)| #
                                             (ini<=CNAhomd$End&ini>=CNAhomd$Start)| # 
                                             (ini<=CNAhomd$Start& end>=CNAhomd$End )| # 
                                             (end>= CNAhomd$Start & end<=CNAhomd$End ) ) ) 
  idhetd = which(CNAhetd$Chromosome==chr&( (ini>=CNAhetd$Start&end<=CNAhetd$End)| #
                                             (ini<=CNAhetd$End&ini>=CNAhetd$Start)| # 
                                             (ini<=CNAhetd$Start& end>=CNAhetd$End )| # 
                                             (end>= CNAhetd$Start & end<=CNAhetd$End ) ) ) 
  
 
  
  samplesNEUT[as.character(CNAneu$Sample[idneu]),i]=T
  samplesAmp[as.character(CNAamp$Sample[idamp]),i]=T
  samplesGAIN[as.character(CNAgain$Sample[idgain]),i]=T
  samplesHETD[as.character(CNAhetd$Sample[idhetd]),i]=T
  samplesHOMD[as.character(CNAhomd$Sample[idhomd]),i]=T
 
  
}

save(samplesNEUT,samplesAmp,samplesGAIN,samplesHETD,samplesHOMD,samplesNUM,r2,file="ESCA_SC/CNA_samples_annotationTCGA_ESCA.RData")


allSamples = rownames(samplesAmp)
clinic = read.table("20151001_ESCA_STAD_Master_Patient_Table.tsv",header = T,sep = "\t",stringsAsFactors = F)

escaSccSamples = intersect(clinic$barcode [ clinic$EC==1],clinic$barcode [ clinic$SCC==1])

getOutliersmRNA(CNAFile = "ESCA_SC/CNA_samples_annotationTCGA_ESCA.RData",exp_met = mRNA.Data,Des = mRNA.Des[,1],outputName = "ESCA_SC/outliers_.RData",samples = escaSccSamples,cores = 10,r2 = r2)

#### HNSQ
load("HNSC_illuminahiseq_rnaseqv2__rsem.genes.rda")

types = substr(colnames(Data),start = 14,stop = 15)
mRNA.Data = Data[,which(types == "01")]
colnames(mRNA.Data)=substr(colnames(mRNA.Data),1,12)
mRNA.Des = Des

# first get the thresholds from raw data (not gene anotated), and the n
CNAraw = read.table("HNSC__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Feb-12-2015.txt",sep="\t",header = T) 
type = substr(CNAraw$Sample,start = 14,stop = 15)
CNAraw = CNAraw[ type=="01",]
CNAraw$Sample = substr(CNAraw$Sample,1,12) 

CNAsamples = unique(CNAraw$Sample)
CNAraw$call = NA
for(i in 1:length(CNAsamples)){
  sample1 = CNAraw[ CNAraw$Sample== as.character(CNAsamples[i]),6]
  qunt <- quantile(sample1, c(0.075, 0.125, 0.15, 0.25, 0.275, 0.315, 0.35, 0.375, 0.5, 0.625, 0.65, 0.685, 0.725, 0.75, 0.85, 0.875, 0.925), na.rm=TRUE)
  median = qunt[9]
  med = median(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  all.sds0.5 = sd(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  ampsc =1* (sample1>(med+2*all.sds0.5)) + 1*  (sample1>(med+6*all.sds0.5)) -1* (sample1<(med-2.5*all.sds0.5)) -1 * (sample1<(med-5*all.sds0.5))
  CNAraw$call[CNAraw$Sample== as.character(CNAsamples[i])]=ampsc
}

library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene_id","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])


CNAamp = CNAraw[CNAraw$call=="2",]
CNAneu = CNAraw[CNAraw$call=="0",]
CNAgain = CNAraw[CNAraw$call=="1",]
CNAhomd = CNAraw[CNAraw$call=="-2",]
CNAhetd = CNAraw[CNAraw$call=="-1",]

nsamples = length(CNAsamples)
ngenes = length(geneNames)
samplesAmp = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesAmp)=geneNames
rownames(samplesAmp)=CNAsamples


samplesGAIN = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesGAIN)=geneNames
rownames(samplesGAIN)=CNAsamples

samplesHETD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHETD)=geneNames
rownames(samplesHETD)=CNAsamples

samplesHOMD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHOMD)=geneNames
rownames(samplesHOMD)=CNAsamples

samplesNEUT = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesNEUT)=geneNames
rownames(samplesNEUT)=CNAsamples

samplesNUM = matrix(nrow = nsamples,ncol = ngenes,data=NA)
colnames(samplesNUM)=geneNames
rownames(samplesNUM)=CNAsamples

for (i in 1:length(geneNames)){
  idr = which(r2$hgnc_symbol==geneNames[i])
  chr = as.character(r2$chrom[idr])
  ini = as.numeric(r2$start[idr])
  end = as.numeric(r2$end[idr])
  
  
  
  idamp = which(CNAamp$Chromosome==chr&( (ini>=CNAamp$Start&end<=CNAamp$End)| #
                                           (ini<=CNAamp$End&ini>=CNAamp$Start)| # 
                                           (ini<=CNAamp$Start& end>=CNAamp$End )| # 
                                           (end>= CNAamp$Start & end<=CNAamp$End ) ) ) 
  
  #idamp = which(CNAamp$Chromosome==chr&( (end>=CNAamp$Start&end<=CNAamp$End)|(ini<=CNAamp$End&ini>=CNAamp$Start) ) )  
  idgain = which(CNAgain$Chromosome==chr&( (ini>=CNAgain$Start&end<=CNAgain$End)| #
                                             (ini<=CNAgain$End&ini>=CNAamp$Start)| # 
                                             (ini<=CNAgain$Start& end>=CNAgain$End )| # 
                                             (end>= CNAgain$Start & end<=CNAgain$End ) ) ) 
  idneu = which(CNAneu$Chromosome==chr&( (ini>=CNAneu$Start&end<=CNAneu$End)| #
                                           (ini<=CNAneu$End&ini>=CNAneu$Start)| # 
                                           (ini<=CNAneu$Start& end>=CNAneu$End )| # 
                                           (end>= CNAneu$Start & end<=CNAneu$End ) ) ) 
  idhomd = which(CNAhomd$Chromosome==chr&( (ini>=CNAhomd$Start&end<=CNAhomd$End)| #
                                             (ini<=CNAhomd$End&ini>=CNAhomd$Start)| # 
                                             (ini<=CNAhomd$Start& end>=CNAhomd$End )| # 
                                             (end>= CNAhomd$Start & end<=CNAhomd$End ) ) ) 
  idhetd = which(CNAhetd$Chromosome==chr&( (ini>=CNAhetd$Start&end<=CNAhetd$End)| #
                                             (ini<=CNAhetd$End&ini>=CNAhetd$Start)| # 
                                             (ini<=CNAhetd$Start& end>=CNAhetd$End )| # 
                                             (end>= CNAhetd$Start & end<=CNAhetd$End ) ) ) 
  
  
  samplesNEUT[as.character(CNAneu$Sample[idneu]),i]=T
  samplesAmp[as.character(CNAamp$Sample[idamp]),i]=T
  samplesGAIN[as.character(CNAgain$Sample[idgain]),i]=T
  samplesHETD[as.character(CNAhetd$Sample[idhetd]),i]=T
  samplesHOMD[as.character(CNAhomd$Sample[idhomd]),i]=T

  
}

save(samplesNEUT,samplesAmp,samplesGAIN,samplesHETD,samplesHOMD,samplesNUM,r2,file="HNSC/CNA_samples_annotationTCGA_HNSC.RData")


allSamples = rownames(samplesAmp)
getOutliersmRNA(CNAFile = "HNSC/CNA_samples_annotationTCGA_HNSC.RData",exp_met = mRNA.Data,Des = mRNA.Des[,1],outputName = "HNSC/outliers_.RData",samples = allSamples,cores = 10,r2 = r2)



# PAAD
load("PAAD__illuminahiseq_rnaseqv2__rsem_norm.genes.rda")

types = substr(colnames(Data),start = 14,stop = 15)
mRNA.Data = Data[,which(types == "01")]
colnames(mRNA.Data)=substr(colnames(mRNA.Data),1,12)
mRNA.Des = Des

# first get the thresholds from raw data (not gene anotated), and the n
CNAraw = read.table("PAAD__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Oct--18:55:19.txt",sep="\t",header = T) 
type = substr(CNAraw$Sample,start = 14,stop = 15)
CNAraw = CNAraw[ type=="01",]
CNAraw$Sample = substr(CNAraw$Sample,1,12) 

CNAsamples = unique(CNAraw$Sample)
CNAraw$call = NA
for(i in 1:length(CNAsamples)){
  sample1 = CNAraw[ CNAraw$Sample== as.character(CNAsamples[i]),6]
  qunt <- quantile(sample1, c(0.075, 0.125, 0.15, 0.25, 0.275, 0.315, 0.35, 0.375, 0.5, 0.625, 0.65, 0.685, 0.725, 0.75, 0.85, 0.875, 0.925), na.rm=TRUE)
  median = qunt[9]
  med = median(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  all.sds0.5 = sd(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  ampsc =1* (sample1>(med+2*all.sds0.5)) + 1*  (sample1>(med+6*all.sds0.5)) -1* (sample1<(med-2.5*all.sds0.5)) -1 * (sample1<(med-5*all.sds0.5))
  CNAraw$call[CNAraw$Sample== as.character(CNAsamples[i])]=ampsc
}

library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])


CNAamp = CNAraw[CNAraw$call=="2",]
CNAneu = CNAraw[CNAraw$call=="0",]
CNAgain = CNAraw[CNAraw$call=="1",]
CNAhomd = CNAraw[CNAraw$call=="-2",]
CNAhetd = CNAraw[CNAraw$call=="-1",]

nsamples = length(CNAsamples)
ngenes = length(geneNames)
samplesAmp = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesAmp)=geneNames
rownames(samplesAmp)=CNAsamples


samplesGAIN = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesGAIN)=geneNames
rownames(samplesGAIN)=CNAsamples

samplesHETD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHETD)=geneNames
rownames(samplesHETD)=CNAsamples

samplesHOMD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHOMD)=geneNames
rownames(samplesHOMD)=CNAsamples

samplesNEUT = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesNEUT)=geneNames
rownames(samplesNEUT)=CNAsamples

samplesNUM = matrix(nrow = nsamples,ncol = ngenes,data=NA)
colnames(samplesNUM)=geneNames
rownames(samplesNUM)=CNAsamples

for (i in 1:length(geneNames)){
  idr = which(r2$hgnc_symbol==geneNames[i])
  chr = as.character(r2$chrom[idr])
  ini = as.numeric(r2$start[idr])
  end = as.numeric(r2$end[idr])
  
  idamp = which(CNAamp$Chromosome==chr&( (ini>=CNAamp$Start&end<=CNAamp$End)| #
                                           (ini<=CNAamp$End&ini>=CNAamp$Start)| # 
                                           (ini<=CNAamp$Start& end>=CNAamp$End )| # 
                                           (end>= CNAamp$Start & end<=CNAamp$End ) ) ) 
  

  idgain = which(CNAgain$Chromosome==chr&( (ini>=CNAgain$Start&end<=CNAgain$End)| #
                                             (ini<=CNAgain$End&ini>=CNAamp$Start)| # 
                                             (ini<=CNAgain$Start& end>=CNAgain$End )| # 
                                             (end>= CNAgain$Start & end<=CNAgain$End ) ) ) 
  idneu = which(CNAneu$Chromosome==chr&( (ini>=CNAneu$Start&end<=CNAneu$End)| #
                                           (ini<=CNAneu$End&ini>=CNAneu$Start)| # 
                                           (ini<=CNAneu$Start& end>=CNAneu$End )| # 
                                           (end>= CNAneu$Start & end<=CNAneu$End ) ) ) 
  idhomd = which(CNAhomd$Chromosome==chr&( (ini>=CNAhomd$Start&end<=CNAhomd$End)| #
                                             (ini<=CNAhomd$End&ini>=CNAhomd$Start)| # 
                                             (ini<=CNAhomd$Start& end>=CNAhomd$End )| # 
                                             (end>= CNAhomd$Start & end<=CNAhomd$End ) ) ) 
  idhetd = which(CNAhetd$Chromosome==chr&( (ini>=CNAhetd$Start&end<=CNAhetd$End)| #
                                             (ini<=CNAhetd$End&ini>=CNAhetd$Start)| # 
                                             (ini<=CNAhetd$Start& end>=CNAhetd$End )| # 
                                             (end>= CNAhetd$Start & end<=CNAhetd$End ) ) ) 
  

  
  samplesNEUT[as.character(CNAneu$Sample[idneu]),i]=T
  samplesAmp[as.character(CNAamp$Sample[idamp]),i]=T
  samplesGAIN[as.character(CNAgain$Sample[idgain]),i]=T
  samplesHETD[as.character(CNAhetd$Sample[idhetd]),i]=T
  samplesHOMD[as.character(CNAhomd$Sample[idhomd]),i]=T
  
  
}

save(samplesNEUT,samplesAmp,samplesGAIN,samplesHETD,samplesHOMD,samplesNUM,r2,file="./PAAD/CNA_samples_annotationTCGA_PAAD.RData")



allSamples = rownames(samplesAmp)
getOutliersmRNA(CNAFile = "PAAD/CNA_samples_annotationTCGA_PAAD.RData",exp_met = mRNA.Data,Des = mRNA.Des[,1],outputName = "PAAD/outliers_.RData",samples = allSamples,cores = 10,r2 = r2)


## COAD
load("COAD__illuminahiseq_rnaseqv2__genes.normalized_results.rda")

types = substr(colnames(Data),start = 14,stop = 15)
mRNA.Data = Data[,which(types == "01")]
colnames(mRNA.Data)=substr(colnames(mRNA.Data),1,12)
mRNA.Des = Des

# first get the thresholds from raw data (not gene anotated), and the n
CNAraw = read.table("COAD__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Feb-12-2015.txt",sep="\t",header = T) 
type = substr(CNAraw$Sample,start = 14,stop = 15)
CNAraw = CNAraw[ type=="01",]
CNAraw$Sample = substr(CNAraw$Sample,1,12) 

CNAsamples = unique(CNAraw$Sample)
CNAraw$call = NA
for(i in 1:length(CNAsamples)){
  sample1 = CNAraw[ CNAraw$Sample== as.character(CNAsamples[i]),6]
  qunt <- quantile(sample1, c(0.075, 0.125, 0.15, 0.25, 0.275, 0.315, 0.35, 0.375, 0.5, 0.625, 0.65, 0.685, 0.725, 0.75, 0.85, 0.875, 0.925), na.rm=TRUE)
  median = qunt[9]
  med = median(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  all.sds0.5 = sd(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  ampsc =1* (sample1>(med+2*all.sds0.5)) + 1*  (sample1>(med+6*all.sds0.5)) -1* (sample1<(med-2.5*all.sds0.5)) -1 * (sample1<(med-5*all.sds0.5))
  CNAraw$call[CNAraw$Sample== as.character(CNAsamples[i])]=ampsc
}

library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])


CNAamp = CNAraw[CNAraw$call=="2",]
CNAneu = CNAraw[CNAraw$call=="0",]
CNAgain = CNAraw[CNAraw$call=="1",]
CNAhomd = CNAraw[CNAraw$call=="-2",]
CNAhetd = CNAraw[CNAraw$call=="-1",]

nsamples = length(CNAsamples)
ngenes = length(geneNames)
samplesAmp = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesAmp)=geneNames
rownames(samplesAmp)=CNAsamples


samplesGAIN = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesGAIN)=geneNames
rownames(samplesGAIN)=CNAsamples

samplesHETD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHETD)=geneNames
rownames(samplesHETD)=CNAsamples

samplesHOMD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHOMD)=geneNames
rownames(samplesHOMD)=CNAsamples

samplesNEUT = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesNEUT)=geneNames
rownames(samplesNEUT)=CNAsamples

samplesNUM = matrix(nrow = nsamples,ncol = ngenes,data=NA)
colnames(samplesNUM)=geneNames
rownames(samplesNUM)=CNAsamples

for (i in 1:length(geneNames)){
  idr = which(r2$hgnc_symbol==geneNames[i])
  chr = as.character(r2$chrom[idr])
  ini = as.numeric(r2$start[idr])
  end = as.numeric(r2$end[idr])
  
  idamp = which(CNAamp$Chromosome==chr&( (ini>=CNAamp$Start&end<=CNAamp$End)| #
                                           (ini<=CNAamp$End&ini>=CNAamp$Start)| # 
                                           (ini<=CNAamp$Start& end>=CNAamp$End )| # 
                                           (end>= CNAamp$Start & end<=CNAamp$End ) ) ) 
  
  #idamp = which(CNAamp$Chromosome==chr&( (end>=CNAamp$Start&end<=CNAamp$End)|(ini<=CNAamp$End&ini>=CNAamp$Start) ) )  
  idgain = which(CNAgain$Chromosome==chr&( (ini>=CNAgain$Start&end<=CNAgain$End)| #
                                             (ini<=CNAgain$End&ini>=CNAamp$Start)| # 
                                             (ini<=CNAgain$Start& end>=CNAgain$End )| # 
                                             (end>= CNAgain$Start & end<=CNAgain$End ) ) ) 
  idneu = which(CNAneu$Chromosome==chr&( (ini>=CNAneu$Start&end<=CNAneu$End)| #
                                           (ini<=CNAneu$End&ini>=CNAneu$Start)| # 
                                           (ini<=CNAneu$Start& end>=CNAneu$End )| # 
                                           (end>= CNAneu$Start & end<=CNAneu$End ) ) ) 
  idhomd = which(CNAhomd$Chromosome==chr&( (ini>=CNAhomd$Start&end<=CNAhomd$End)| #
                                             (ini<=CNAhomd$End&ini>=CNAhomd$Start)| # 
                                             (ini<=CNAhomd$Start& end>=CNAhomd$End )| # 
                                             (end>= CNAhomd$Start & end<=CNAhomd$End ) ) ) 
  idhetd = which(CNAhetd$Chromosome==chr&( (ini>=CNAhetd$Start&end<=CNAhetd$End)| #
                                             (ini<=CNAhetd$End&ini>=CNAhetd$Start)| # 
                                             (ini<=CNAhetd$Start& end>=CNAhetd$End )| # 
                                             (end>= CNAhetd$Start & end<=CNAhetd$End ) ) ) 
  
  # idgen = which(CNAraw$Chromosome==chr&( (ini>=CNAraw$Start&end<=CNAraw$End)| #
  #                                          (ini<=CNAraw$End&ini>=CNAraw$Start)| # 
  #                                          (ini<=CNAraw$Start& end>=CNAraw$End )| # 
  #                                          (end>= CNAraw$Start & end<=CNAraw$End ) ) ) 
  # 
  samplesNEUT[as.character(CNAneu$Sample[idneu]),i]=T
  samplesAmp[as.character(CNAamp$Sample[idamp]),i]=T
  samplesGAIN[as.character(CNAgain$Sample[idgain]),i]=T
  samplesHETD[as.character(CNAhetd$Sample[idhetd]),i]=T
  samplesHOMD[as.character(CNAhomd$Sample[idhomd]),i]=T
  #samplesNUM[ as.character(CNAraw$Sample[idgen]),i]=CNAraw$Segment_Mean[idgen]
  
}

save(samplesNEUT,samplesAmp,samplesGAIN,samplesHETD,samplesHOMD,samplesNUM,r2,file="./COAD/CNA_samples_annotationTCGA_COAD.RData")

allSamples = rownames(samplesAmp)
getOutliersmRNA(CNAFile = "COAD/CNA_samples_annotationTCGA_COAD.RData",exp_met = mRNA.Data,Des = mRNA.Des[,1],outputName = "COAD/outliers_.RData",samples = allSamples,cores = 10,r2 = r2)

#STAD 
load("STAD__illuminahiseq_rnaseqv2__rsem_gene.quantification_normalized.rda")

types = substr(colnames(Data),start = 14,stop = 15)
mRNA.Data = Data[,which(types == "01")]
colnames(mRNA.Data)=substr(colnames(mRNA.Data),1,12)
mRNA.Des = Des

# first get the thresholds from raw data (not gene anotated), and the n
CNAraw = read.table("STAD__broad.mit.edu__genome_wide_snp_6__nocnv_hg19__Feb-12-2015.txt",sep="\t",header = T) 
type = substr(CNAraw$Sample,start = 14,stop = 15)
CNAraw = CNAraw[ type=="01",]
CNAraw$Sample = substr(CNAraw$Sample,1,12) 

CNAsamples = unique(CNAraw$Sample)
CNAraw$call = NA
for(i in 1:length(CNAsamples)){
  sample1 = CNAraw[ CNAraw$Sample== as.character(CNAsamples[i]),6]
  qunt <- quantile(sample1, c(0.075, 0.125, 0.15, 0.25, 0.275, 0.315, 0.35, 0.375, 0.5, 0.625, 0.65, 0.685, 0.725, 0.75, 0.85, 0.875, 0.925), na.rm=TRUE)
  median = qunt[9]
  med = median(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  all.sds0.5 = sd(sample1[sample1>=qunt[4] & sample1 <= qunt[14]], na.rm=TRUE)
  ampsc =1* (sample1>(med+2*all.sds0.5)) + 1*  (sample1>(med+6*all.sds0.5)) -1* (sample1<(med-2.5*all.sds0.5)) -1 * (sample1<(med-5*all.sds0.5))
  CNAraw$call[CNAraw$Sample== as.character(CNAsamples[i])]=ampsc
}

library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])


CNAamp = CNAraw[CNAraw$call=="2",]
CNAneu = CNAraw[CNAraw$call=="0",]
CNAgain = CNAraw[CNAraw$call=="1",]
CNAhomd = CNAraw[CNAraw$call=="-2",]
CNAhetd = CNAraw[CNAraw$call=="-1",]

nsamples = length(CNAsamples)
ngenes = length(geneNames)
samplesAmp = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesAmp)=geneNames
rownames(samplesAmp)=CNAsamples


samplesGAIN = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesGAIN)=geneNames
rownames(samplesGAIN)=CNAsamples

samplesHETD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHETD)=geneNames
rownames(samplesHETD)=CNAsamples

samplesHOMD = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesHOMD)=geneNames
rownames(samplesHOMD)=CNAsamples

samplesNEUT = matrix(nrow = nsamples,ncol = ngenes,data=F)
colnames(samplesNEUT)=geneNames
rownames(samplesNEUT)=CNAsamples

samplesNUM = matrix(nrow = nsamples,ncol = ngenes,data=NA)
colnames(samplesNUM)=geneNames
rownames(samplesNUM)=CNAsamples

for (i in 1:length(geneNames)){
  idr = which(r2$hgnc_symbol==geneNames[i])
  chr = as.character(r2$chrom[idr])
  ini = as.numeric(r2$start[idr])
  end = as.numeric(r2$end[idr])
  
  idamp = which(CNAamp$Chromosome==chr&( (ini>=CNAamp$Start&end<=CNAamp$End)| #
                                           (ini<=CNAamp$End&ini>=CNAamp$Start)| # 
                                           (ini<=CNAamp$Start& end>=CNAamp$End )| # 
                                           (end>= CNAamp$Start & end<=CNAamp$End ) ) ) 
  
  
  idgain = which(CNAgain$Chromosome==chr&( (ini>=CNAgain$Start&end<=CNAgain$End)| #
                                             (ini<=CNAgain$End&ini>=CNAamp$Start)| # 
                                             (ini<=CNAgain$Start& end>=CNAgain$End )| # 
                                             (end>= CNAgain$Start & end<=CNAgain$End ) ) ) 
  idneu = which(CNAneu$Chromosome==chr&( (ini>=CNAneu$Start&end<=CNAneu$End)| #
                                           (ini<=CNAneu$End&ini>=CNAneu$Start)| # 
                                           (ini<=CNAneu$Start& end>=CNAneu$End )| # 
                                           (end>= CNAneu$Start & end<=CNAneu$End ) ) ) 
  idhomd = which(CNAhomd$Chromosome==chr&( (ini>=CNAhomd$Start&end<=CNAhomd$End)| #
                                             (ini<=CNAhomd$End&ini>=CNAhomd$Start)| # 
                                             (ini<=CNAhomd$Start& end>=CNAhomd$End )| # 
                                             (end>= CNAhomd$Start & end<=CNAhomd$End ) ) ) 
  idhetd = which(CNAhetd$Chromosome==chr&( (ini>=CNAhetd$Start&end<=CNAhetd$End)| #
                                             (ini<=CNAhetd$End&ini>=CNAhetd$Start)| # 
                                             (ini<=CNAhetd$Start& end>=CNAhetd$End )| # 
                                             (end>= CNAhetd$Start & end<=CNAhetd$End ) ) ) 
  
  
  samplesNEUT[as.character(CNAneu$Sample[idneu]),i]=T
  samplesAmp[as.character(CNAamp$Sample[idamp]),i]=T
  samplesGAIN[as.character(CNAgain$Sample[idgain]),i]=T
  samplesHETD[as.character(CNAhetd$Sample[idhetd]),i]=T
  samplesHOMD[as.character(CNAhomd$Sample[idhomd]),i]=T
 
  
}

save(samplesNEUT,samplesAmp,samplesGAIN,samplesHETD,samplesHOMD,samplesNUM,r2,file="./STAD/CNA_samples_annotationTCGA_STAD.RData")


## STAD

load("LUAD__illuminahiseq_rnaseqv2__rsem.genes.normalized.rda")

types = substr(colnames(Data),start = 14,stop = 15)
mRNA.Data = Data[,which(types == "01")]
colnames(mRNA.Data)=substr(colnames(mRNA.Data),1,12)
mRNA.Des = Des

load("STAD/CNA_samples_annotationTCGA_STAD.RData")


library(biomaRt)

ensembl = useMart("ENSEMBL_MART_ENSEMBL",
                  dataset="hsapiens_gene_ensembl",
                  host="grch37.ensembl.org",
                  path="/biomart/martservice",
                  archive=FALSE)

r = getBM(attributes=c("hgnc_symbol","entrezgene","ensembl_gene_id","chromosome_name","start_position","end_position"),
          filters="hgnc_symbol",
          values=unique(mRNA.Des[,1]),mart=ensembl)

r2 = r[ r$chromosome_name %in% c(1:22,"X","Y"),]

colnames(r2)=c("hgnc_symbol",     "entrezgene","ensembl_gene_id", "chrom" ,"start",  "end")

idIn = which(mRNA.Des[,1] %in% r2$hgnc_symbol)



mRNA.Data = mRNA.Data[idIn,]
mRNA.Des = mRNA.Des[idIn,]


geneNames = as.character(mRNA.Des[,"GeneSymbol"])


allSamples = rownames(samplesAmp)
getOutliersmRNA(CNAFile = "STAD/CNA_samples_annotationTCGA_STAD.RData",exp_met = mRNA.Data,Des = mRNA.Des[,1],outputName = "STAD/outliers_.RData",samples = allSamples,cores = 10,r2 = r2)

