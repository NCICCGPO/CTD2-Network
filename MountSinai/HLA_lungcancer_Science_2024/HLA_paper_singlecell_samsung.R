library(Seurat)
library(readxl)
library(dplyr)
library(cowplot)
library(ggplot2)
dyn.load("/hpc/packages/minerva-centos7/fftw/3.3.9/lib/libfftw3.so.3")
library(metap)
library(msigdbr)
library(singleseqgset)
library(heatmap3)
library(stringr)
library(ggpubr)

#### READ IN FILES and annotate #### 

samsung_raw <- readRDS("samsung/GSE131907_Lung_Cancer_raw_UMI_matrix.rds")
samsung_meta <- read.table("samsung/GSE131907_Lung_Cancer_cell_annotation.txt",header=TRUE,sep="\t")
samsung_info <- read.csv('samsung/samsung_info.csv', header=TRUE, stringsAsFactors = FALSE)

# samples to keep 
# normal samples to use
normal <- c("LUNG_N01","LUNG_N06","LUNG_N08","LUNG_N09","LUNG_N18","LUNG_N19","LUNG_N20","LUNG_N28","LUNG_N30","LUNG_N31", "LUNG_N34")

# tumor samples to use 
tumor <- c("LUNG_T06","LUNG_T08","LUNG_T09","LUNG_T18","LUNG_T19", "LUNG_T20","LUNG_T25","LUNG_T28","LUNG_T30","LUNG_T31","LUNG_T34","EBUS_06","EBUS_28","EBUS_49","BRONCHO_58")

# define smoking category 
non_smoker <- c("BRONCHO_58","LUNG_N01","LUNG_N08","LUNG_N19","LUNG_N30","LUNG_N34","LUNG_T01","LUNG_T08","LUNG_T19","LUNG_T30","LUNG_T34")
smoker <- c("EBUS_06","EBUS_28","EBUS_49","LUNG_N06","LUNG_N09","LUNG_N18","LUNG_N19","LUNG_N20","LUNG_N28","LUNG_N31","LUNG_T06","LUNG_T09","LUNG_T18","LUNG_T19","LUNG_T20","LUNG_T25","LUNG_T28","LUNG_T31")
non_smoker_full <- samsung_info %>% filter(Smoking == "Never") %>% select(Samples)
smoker_full <- samsung_info %>% filter(Smoking == "Ex" | Smoking == 'Cur') %>% select(Samples)
past_smoker_full <- samsung_info %>% filter(Smoking == "Ex") %>% select(Samples)
current_smoker_full <- samsung_info %>% filter(Smoking == "Cur") %>% select(Samples)

#### process raw data into Seurat object ####
# create seurat object
samsung_seurat <- CreateSeuratObject(counts = samsung_raw)
samsung_seurat[["percent.mt"]] <- PercentageFeatureSet(samsung_seurat, pattern = "^MT-")

# check filtering - fits with what was written in the paper  
range(samsung_seurat@meta.data['nFeature_RNA']) # 200 9750
range(samsung_seurat@meta.data['nCount_RNA']) # 1000 148044
range(samsung_seurat@meta.data['percent.mt']) # 0 20
VlnPlot(samsung_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# normalize and scale
samsung_seurat <- NormalizeData(samsung_seurat, normalization.method = "LogNormalize", scale.factor = 10000)
samsung_seurat <- FindVariableFeatures(samsung_seurat, selection.method = "vst", nfeatures = 2000) # need this for PCA 
samsung_seurat <- ScaleData(samsung_seurat,do.scale=FALSE, do.center=TRUE, scale.max=10)

# dim reduction 
samsung_seurat <- RunPCA(samsung_seurat, features = VariableFeatures(object = samsung_seurat))
samsung_seurat <- RunUMAP(samsung_seurat, dims = 1:10)

##### adding additional metadata ####
samsung_seurat@meta.data$Index<- rownames(samsung_seurat@meta.data)
samsung_seurat@meta.data <- merge(samsung_seurat@meta.data,samsung_meta,by='Index')
samsung_seurat@meta.data <- merge(samsung_seurat@meta.data,samsung_info,by.x='Sample',by.y='Samples')

Idents(samsung_seurat) <- samsung_seurat@meta.data[["Cell_type"]]
samsung_seurat$smoke_binary <- ifelse(samsung_seurat@meta.data$Sample %in% non_smoker_full$Samples,'non_smoker','')
samsung_seurat$smoke_binary <- ifelse(samsung_seurat@meta.data$Sample %in% smoker_full$Samples,'smoker',samsung_seurat$smoke_binary)
samsung_seurat$smoke <- ifelse(samsung_seurat@meta.data$Sample %in% past_smoker_full$Samples,'past','')
samsung_seurat$smoke <- ifelse(samsung_seurat@meta.data$Sample %in% current_smoker_full$Samples,'current',samsung_seurat$smoke)
samsung_seurat$smoke <- ifelse(samsung_seurat@meta.data$Sample %in% non_smoker_full$Samples,'never',samsung_seurat$smoke)
samsung_seurat$Sample_type <- ifelse(samsung_seurat@meta.data$Sample_Origin == 'nLung','normal','')
samsung_seurat$Sample_type <- ifelse(samsung_seurat@meta.data$Sample_Origin == 'tLung'|samsung_seurat@meta.data$Sample_Origin == 'tL/B','tumor',samsung_seurat$Sample_type)
samsung_seurat$Sample_type <- ifelse(samsung_seurat@meta.data$Sample_Origin != 'tLung'&samsung_seurat@meta.data$Sample_Origin != 'tL/B'& samsung_seurat@meta.data$Sample_Origin != 'nLung','exclude',samsung_seurat$Sample_type)
samsung_seurat$lung <- ifelse(samsung_seurat@meta.data$Sample_Origin == 'tLung'|samsung_seurat@meta.data$Sample_Origin == 'tL/B'|samsung_seurat@meta.data$Sample_Origin == 'nLung','lung_sample','non_lung')

#### plotting umaps #### 

umap <- as.data.frame(samsung_seurat@reductions$umap@cell.embeddings)
umap$Index <- rownames(umap)
umap <- merge(umap,samsung_meta,by='Index')
umap <- merge(umap,samsung_info,by.x='Sample',by.y='Samples')

# only smoker umap
p <- ggplot(umap %>% filter(lung=='lung'), aes(x=UMAP_1,y= UMAP_2,group=Cell_type))
p <- p + geom_point(aes(color=Cell_type))
p + facet_wrap(~Sample_type)

# only nonsmoker umap
p <- ggplot(umap, aes(x=UMAP_1,y= UMAP_2,group=Cell_type))
p <- p + geom_point(aes(color=Cell_type))

# normal only umap
p <- ggplot(umap, aes(x=UMAP_1,y= UMAP_2))
p <- p + geom_point(aes(color=Cell_type))

# tumor only umap 
p <- ggplot(umap %>% filter(), aes(x=UMAP_1,y= UMAP_2))

# log2 normalized 
lognorm <- GetAssayData(samsung_seurat, 'data')
lognorm_hla <- lognorm[rownames(lognorm) %in% markers.to.plot, ]
lognorm_hla <- as.matrix(lognorm_hla)
lognorm_hla.t <- t(lognorm_hla)
lognorm_hla.t.df <- as.data.frame(lognorm_hla.t)
lognorm_hla.t.df$Index <- rownames(lognorm_hla.t.df) 
rownames(lognorm_hla.t.df) <- NULL 
lognorm_hla.t.df <- merge(Samsung_sc_umap_lognorm_plus_meta_Dec13, lognorm_hla.t.df, by='Index')
#saveRDS(lognorm_hla.t.df, file="Samsung_sc_umap_lognorm_plus_meta_Dec13.rds")

#### feature plots ####
Idents(samsung_seurat) <- samsung_seurat@meta.data[["lung"]]
lung <- FetchData(samsung_seurat, vars = 'ident')
lung_cells <- grep("lung_sample", lung$ident, ignore.case=T)

Idents(samsung_seurat) <- samsung_seurat@meta.data[["Sample_type"]]
tumor <- FetchData(samsung_seurat, vars = 'ident')
tumor_cells <- grep("tumor", tumor$ident, ignore.case=T) # automatically only lung bc otherwise says exclude 
normal_cells <-  grep("normal", tumor$ident, ignore.case=T) # automatically only lung bc otherwise says exclude 
exclude_cells <- grep("exclude", tumor$ident, ignore.case=T) 

markers.to.plot <- c("HLA-A", "HLA-B","HLA-C","HLA-DRB1",'HLA-DQA1','HLA-DQB1','HLA-DPA1','HLA-DPB1')
Idents(samsung_seurat) <- samsung_seurat@meta.data[["Cell_type"]]
FeaturePlot(samsung_seurat, features=f,cells=lung_cells, pt.size=2, split.by = 'smoke_binary',label=TRUE, reduction = "umap")

#### DEGs smokers vs nonsmokers in normal cells #### 
#readRDS('/sc/arion/projects/samstr01a/HLAproject/singlecell/samsung/DEG_mast_cell_lineage_vars.rds')

# idents column incorporating smoking and tumor together - will compare smoker_normal vs non_smoker_normal
samsung_seurat$celllin.smoke.tumor <- paste(samsung_seurat$Cell_subtype, samsung_seurat$smoke_binary, samsung_seurat$Sample_type, sep = "_") # make column with cell type delin 
Idents(samsung_seurat) <- 'celllin.smoke.tumor'

# DEG with latent var correction 
lineages <- c('AT1','AT2','mo-Mac','Alveolar Mac','Monocytes','CD1c+ DCs','Naive CD8+ T','Pleural Mac','Mesothelial cells','Treg','CD4+ Th','Cytotoxic CD8+ T','NK','pDCs')

DEG_mast_celltype_vars_smoker_v_ns_normal <- list()
Idents(samsung_seurat) <- "celltype.smoke.tumor"
for (i in lineages){
  first <- paste(i, 'smoker','normal',sep='_')
  second <- paste(i, 'non_smoker','normal',sep='_')
  print(first)
  print(second)
  DEGs <- FindMarkers(samsung_seurat, ident.1 = first, ident.2 = second, test.use ="MAST",latent.vars=c('Sex','Stage_category','Age'),logfc.threshold = 0)
  DEG_mast_celltype_vars_smoker_v_ns_normal[[i]] <- DEGs
}
#saveRDS(DEG_mast_celltype_vars_smoker_v_ns_normal, file='DEG_mast_celltype_vars_smoker_v_ns_normal.rds')

# plot all genes that came up - not used in paper 
DEG_mast_celltype_vars_smoker_v_ns_tumor[[2]]$gene <- rownames(DEG_mast_celltype_vars_smoker_v_ns_tumor[[2]])
smoker_v_ns_normal <- ggplot(DEG_mast_celltype_vars_smoker_v_ns_tumor[[2]], aes(x=avg_log2FC,y=-log10(p_val_adj))) + geom_point() + ylim(0,100) + ggtitle('DEG in normal epithelial cells- smoker v nonsmoker')
smoker_v_ns_normal <- LabelPoints(plot = smoker_v_ns_normal, points = markers.to.plot,repel = TRUE) + geom_hline(yintercept=-log10(0.001), linetype="dashed", color = "red")
smoker_v_ns_normal







