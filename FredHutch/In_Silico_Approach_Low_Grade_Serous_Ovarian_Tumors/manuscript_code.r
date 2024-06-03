

#2.1. Dataset Selection
####data1 from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18520
####same procedure for data2-9
#BiocManager::install("affy")
#BiocManager::install("Affyhgu133Plus2Expr")
#BiocManager::install("hgu133plus2cdf")
#BiocManager::install("hgu133plus2.db")
#BiocManager::install("annotate")
library(affy)
library(Affyhgu133Plus2Expr)
library(hgu133plus2cdf)
library(hgu133plus2.db)
library(annotate)

#2.2. Microarray Gene-Expression Data Preprocessing
affybatch1 <- ReadAffy ()
sampleNames(affybatch1)
sampleNames(affybatch1) <-c("B1-NOR-1","B1-NOR-2","B1-NOR-3","B1-NOR-4","B1-NOR-5",
                            "B1-NOR-6","B1-NOR-7","B1-NOR-8","B1-NOR-9","B1-NOR-10")
eset1<-rma(affybatch1)
par(mfrow=c(1,2))
# Generates a box plot of un-normalized log intensity values.
boxplot(affybatch1,col="red",main="Raw Data")
# Generates a box plot of normalized log intensity values.
boxplot(data.frame(exprs(eset1)), col="blue", main="Normalized Data") 
write.exprs(eset1, file="mydata1.txt")
probeIDs <- featureNames(eset1)
length(probeIDs)
only_exp1 <- exprs(eset1)
dim(only_exp1)

#2.3. Data Aggregation and Data Integration with Batch-Effect Correction
e1<-as.data.frame(eset1@assayData[["exprs"]])
e2<-as.data.frame(eset2@assayData[["exprs"]])
e3<-as.data.frame(eset3@assayData[["exprs"]])
e4<-as.data.frame(eset4@assayData[["exprs"]])
e5<-as.data.frame(eset5@assayData[["exprs"]])
e6<-as.data.frame(eset6@assayData[["exprs"]])
e7<-as.data.frame(eset7@assayData[["exprs"]])
e8<-as.data.frame(eset8@assayData[["exprs"]])
e9<-as.data.frame(eset9@assayData[["exprs"]])

annotation(affybatch1)
geneAnnot <- data.frame(SYMBOL=sapply(contents(hgu133plus2ENTREZID), paste, collapse=", "))
dforagg <- as.data.frame(cbind(as.vector(geneAnnot),e1,e2,e3,e4,e5,e6,e7,e8,e9))
dataAggAnno<-aggregate(dforagg[, -1],
                       by = list(Gene = dforagg$SYMBOL),
                       FUN = median,
                       na.rm = TRUE)
rownames(dataAggAnno)=as.vector(dataAggAnno$Gene)
dataAgg<-dataAggAnno[-20858, -1]
write.table(dataAgg,"dataAgg.csv")

################3.step combat batch
#BiocManager::install("sva")
library(sva)
pheno <- data.frame("sample" = 1:107, 
                    "batch" = pheno2$batch,
                    "cancer" = pheno2$cancer)
rownames(pheno)<-colnames(dataAgg)
batch = pheno$batch
mod = model.matrix(~as.factor(cancer), data=pheno)
combat_total_exp2 = ComBat(dat=dataAgg, batch=batch, mod=mod, par.prior=TRUE,prior.plots=FALSE)
write.table(combat_total_exp2,"combat2_dataAggBatch107.tsv")

combat_total_exp<-combat_total_exp2[,c(1:10,32:38,49:50,82:93,98:107,
                                       11:18,64:81,94:97,
                                       19:31,39:48,51:63)]
colnames(combat_total_exp)
write.table(combat_total_exp,"combat_dataAggBatch107.tsv")

#2.4. Coexpression Network Analysis
library(xlsx)
library(flashClust)
library(WGCNA)
library(impute)
#BiocManager::install('GO.db')
#BiocManager::install('impute')
allowWGCNAThreads()

phenotype<-data.frame(
  Normal=c(rep(1, 41),rep(0,66)),
  Borderline=c(rep(0, 41),rep(1,30),rep(0,36)),
  Low_Serous=c(rep(0,71),rep(1,36)))
  
data <-read.table(file = 'combat_dataAggBatch107.tsv', sep=" ")
data<-t(data) # rows are samples, columns are genes
powers = c(1:10);
sft=pickSoftThreshold(data,dataIsExpr = TRUE,powerVector = powers,corFnc = cor,corOptions = list(use = 'p'),networkType = "unsigned")

par(mfrow = c(1, 2))
plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", ylab = "SFT, signed R^2", type = "n", main = paste("Scale independence"))
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2], 
     labels = powers, col = "red")
abline(h = 0.8, col = "red")
plot(sft$fitIndices[, 1], sft$fitIndices[, 5], type = "n", xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "red")

A = adjacency(data, power = 8)
dissTOM = TOMdist(A)
geneTree = flashClust(as.dist(dissTOM), method = "average")
moduleLabelsManual1 = cutreeDynamic(dendro = geneTree, distM = dissTOM, method = "hybrid",
                                    deepSplit = 2, pamRespectsDendro = F, minClusterSize = 30)
table(moduleLabelsManual1)
moduleColorsManual1 = labels2colors(moduleLabelsManual1)
plotDendroAndColors(geneTree, moduleColorsManual1, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
MEList = moduleEigengenes(data, colors = moduleColorsManual1)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
windows()
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.25
abline(h=MEDissThres, col = "red")
MET = orderMEs(MEs)
windows()
plotEigengeneNetworks(MET, "", marDendro = c(0, 4, 1, 2),
                      marHeatmap = c(3,4, 1, 2), cex.lab = 0.8, xLabelsAngle = 90)


merge = mergeCloseModules(data, moduleColorsManual1, cutHeight = MEDissThres, verbose = 3)
moduleColorsManual2 = merge$colors;
MEsManual = merge$newMEs;
nGenes = ncol(data)
nSamples = nrow(data)
MEs3 = moduleEigengenes(data, moduleColorsManual2)$eigengenes
MEs_ordered3 = orderMEs(MEs3)
modTraitCor3 = cor(MEs_ordered3, phenotype, use = "p")
modTraitP3 = corPvalueStudent(modTraitCor3, nSamples)
textMatrix3 = paste(signif(modTraitCor3, 2), "\n(", signif(modTraitP3, 1), ")", 
                    sep = "")
dim(textMatrix3) = dim(modTraitCor3)
openxlsx::write.xlsx(textMatrix3, file = "train_all_manualhybrid.xlsx")
windows()
par(mar = c(2, 2, 1, 1));
labeledHeatmap(Matrix = modTraitCor3, xLabels = names(phenotype), yLabels = names(MEs_ordered3), 
               ySymbols = names(MEs_ordered3), colorLabels = T, colors = greenWhiteRed(50), 
               textMatrix = textMatrix3, setStdMargins = FALSE, cex.text = 0.7, zlim = c(-1,1), main = paste("Module-trait relationships for manual hybrid module detection"))
mnames<-as.vector(colnames(MEs_ordered3))
for (i in 1:24)
{
  write.xlsx(colnames(data)[moduleColorsManual2==str_remove_all(mnames[i],"[ME]")],'wgcna_module.xlsx',sheetName =as.character(mnames[i]),append = T,row.names = F)
}

#2.5. Differentially Expressed Genes Analysis
library(limma)
data<-read.table("combat_dataAggBatch107.tsv")
ourData <- data[,c(1:71)]
pheno<-rep(c(' Normal','Borderline'),c(41,30))

design <- model.matrix(~ pheno)
fit <- lmFit(ourData, design)
fit <- eBayes(fit)
analyzed.data<-topTable(fit,n=20857)

nvsb <- analyzed.data[ (abs(analyzed.data[,1]) >= 1.0 & analyzed.data[,5] < 0.05), ]
nvsb$mRNA<-rownames(nvsb)

geneName<-getBM(attributes = c('entrezgene_id','hgnc_symbol' ),             
                filters = 'entrezgene_id', 
                values = nvsb$mRNA,
                mart = mart)
m_all_deg<-merge (nvsb, geneName, by.x = 'mRNA', by.y ='entrezgene_id')
write.xlsx(m_all_deg,file = "limma.xlsx",sheetName = "1nvsb",append=T,row.names = F)

#2.6. Gene-Set Enrichment Analysis
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
geneName<-getBM(attributes = c('entrezgene_id','hgnc_symbol' ),             
                filters = 'entrezgene_id', 
                values = dataDEGsFiltLevel$mRNA,
                mart = mart)
dataDEGsFiltLevel_X1<- dataDEGsFiltLevel[order(dataDEGsFiltLevel$fcS),-c(2:46)]
dataDEGsFiltLevel_X <-merge (dataDEGsFiltLevel_X1, geneName, by.x = 'mRNA', by.y ='entrezgene_id')
dataDEGsFiltLevel_X<-dataDEGsFiltLevel_X[,c(1,5,2:4)]
m_all_deg<- dataDEGsFiltLevel_X[order(dataDEGsFiltLevel_X$fcS),]
m_all_deg$label<-""
m_all_deg$label[m_all_deg$fcS<0]<-"down"
m_all_deg$label[m_all_deg$fcS>0]<-"up"
m_all_deg$group<-"Borderline"
write.xlsx(m_all_deg,file = "1NvsB.xlsx",sheetName = "dataDEGsFiltLevel_X",append=T,row.names = F)

####enrich
library(affy)
library(hgu133plus2cdf)
library(hgu133plus2.db)
library(annotate)
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
combat_total_exp<-read.table(file ="combat_dataAggBatch107.tsv", header = T)

dataAggBatchNB<-combat_total_exp[,c(1:71)]
x <- log2(dataAggBatchNB)
fcS <- rowMeans(dataAggBatchNB[,42:71]) - rowMeans(dataAggBatchNB[, 1:41])
p.val = apply(dataAggBatchNB, 1, function(x) { t.test(x[1:41], x[42:71]) $p.value } )
fdr.pval = p.adjust(p.val, method="fdr")
mRNA<-row.names(combat_total_exp)
analyzed.data <- cbind(mRNA, dataAggBatchNB, fcS, p.val, fdr.pval)
sapply(analyzed.data, class)
analyzed.data <- as.data.frame(apply(analyzed.data, 2, as.numeric))
geneName<-getBM(attributes = c('entrezgene_id','hgnc_symbol' ),             
                filters = 'entrezgene_id', 
                values = analyzed.data$mRNA,
                mart = mart)
m_all_deg<-merge (analyzed.data, geneName, by.x = 'mRNA', by.y ='entrezgene_id')
dataDEGsFiltLevel2 <- analyzed.data[ (abs(analyzed.data[,73]) >= 1.0 & analyzed.data[,75] < 0.05), ]

library(ggplot2)
library(enrichR)
dbs<-listEnrichrDbs()
dbs<-c("KEGG_2019_Human","GO_Biological_Process_2018",
       "MSigDB_Hallmark_2020","HMDB_Metabolites")
eup <- enrichr(m_all_deg$hgnc_symbol, dbs)
###kegg
up<-eup$KEGG_2019_Human[eup$KEGG_2019_Human$Adjusted.P.value<0.05,c(4,1,2,9)]
names(up)[2]<-'Term_Description'
names(up)[1]<-"lowest_p"
up$Up_regulated<-""
up$Down_regulated<-""
Genes<-(strsplit(up$Genes,"\\;"))
for(i in 1:length(Genes))
  {
  a<-subset(m_all_deg, m_all_deg$hgnc_symbol %in% Genes[[i]]) 
  splitlabel<-split (a, a$label)
  up$Up_regulated[i]<-paste0(splitlabel$up$hgnc_symbol, collapse=", ")
  up$Down_regulated[i]<-paste0(splitlabel$down$hgnc_symbol, collapse=", ")
}
up<-up[,-4]
write.xlsx(up,file = "1NvsB.xlsx",sheetName = "kegg",append = T)

#2.7. Drug–Gene Interaction Analysis
d1<-read.xlsx("1NvsBvsM2.xlsx",sheetIndex = 1)
genesUP<-subset(d1,label=="up",select="hgnc_symbol")
genesdown<-subset(d1,label=="down",select="hgnc_symbol")

library(rDGIdb)

resultUP <- queryDGIdb(genesUP$hgnc_symbol,interactionTypes= c("antagonist","antibody", "antisense oligonucleotide","blocker","cleavage",
                                                   "inhibitor" ,"inhibitory allosteric modulator" ,"inverse agonist" ,"negative modulator",
                                                   "partial antagonist", "suppressor"))
write.xlsx(resultUP@detailedResults[,-1],"1NvsBvsM2.xlsx",sheetName = "1NvsBvsM2_drugUP",append = T,row.names = F)
write.xlsx(resultUP@byGene[,c(2,3,5)],"1NvsBvsM2.xlsx",sheetName = "1NvsBvsM2_geneUP",append = T,row.names = F)


resultDOWN <- queryDGIdb(genesdown$hgnc_symbol,
                       interactionTypes= c("activator", "agonist" , "chaperone","cofactor","inducer", "partial agonist",
                                           "positive modulator", "stimulator","vaccine" ))
write.xlsx(resultDOWN@detailedResults[,-1],"1NvsBvsM2.xlsx",sheetName = "1NvsBvsM2_drugDOWN",append = T,row.names = F)
write.xlsx(resultDOWN@byGene[,c(2,3,5)],"1NvsBvsM2.xlsx",sheetName = "1NvsBvsM2_geneDOWN",append = T,row.names = F)

#2.9. Estimation of Tumor-Microenvironment Infiltration
#install.packages("remotes")
remotes::install_github("grst/immunedeconv")

library(dplyr)
library(tidyr)
library(tibble)
library(immunedeconv)
library(xlsx)
library(ggplot2)
library(ggpubr)
library(MASS)
library(reshape2)
library(reshape)

deg2<-read.table('combat_dataAggBatch107.tsv')
res_mcp_counter <- deconvolute(deg2, "mcp_counter")
rownames(res_mcp_counter) <- res_mcp_counter$cell_type
data2<-as.data.frame(t(res_mcp_counter))
data2 <- data2[-1,]
data2 <- as.data.frame(apply(data2, 2, as.numeric))
data2$Phenotype<-rep(c(' Normal','Borderline','LowSerous'),c(41,30,36))
melt_data2 <- melt(data2, id = c("Phenotype"))
colnames(melt_data2)
colnames(melt_data2)<-c("Phenotype", "Cell_type", "Estimated_Proportion")

p <- ggboxplot(melt_data2, x = "Cell_type", y = "Estimated_Proportion",
               color = "Phenotype", palette = "jco",
               add = "jitter")+ theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust = 0.5))#,facet.by = 'variable'
p + stat_compare_means(aes(group = Phenotype), label = "p.signif")
  #stat_compare_means(comparisons = my_comparisons, label =  "p.signif")
p

#####################
# install.packages("ggstatsplot")
library(ggstatsplot)
ggscatterstats(
  data = corrr[,-23],
  x    = ITGB2,
  y    = Monocyte,          
  type = "p")
library(ggpubr)
ggscatter(corrr, 
          x = "ERBB4", 
          y = c("T cell","T cell CD8+","cytotoxicity score" ,"NK cell","B cell",
                                    "Macrophage/Monocyte","Myeloid dendritic cell", "Neutrophil" ,"Endothelial cell", 
                                    "Cancer associated fibroblast"), 
          size = 0.5,
          combine = TRUE, ylab = "Expression",
          color = "Phenotype", palette = "jco",
          add = "reg.line", conf.int = TRUE) +
  stat_cor(aes(color = Phenotype), method = "spearman")