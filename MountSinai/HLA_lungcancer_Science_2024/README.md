# An immunogenetic basis for lung cancer risk

This repository contains the script to reproduce the single-cell analyses and figures in Krishna et al Science 2024 of the Samsung lung cancer dataset in Kim et al Nat Comm 2020. Supplementary files for the single cell analysis results are also housed here.

To repeat the single-cell analyses from our paper, first download the processed single-cell data provided by Kim et al from NCBI Gene Expression Omnibus database (accession code GSE131907). More details of the dataset can be found in their paper. Then follow the steps in HLA_paper_singlecell_samsung.R to process the Seurat object from the raw data, select the relevant cells, and perform the analyses. The generated results files can be used to reproduce the relevant figures as outlined in the script. 
Extended results from DEG analysis comparing genes differentially expressed in smokers versus nonsmokers in all cell types, can also be downloaded from the supplementary materials of our study (Table S24). Detailed methods of our analysis as shown in Figure 6, can be found in the materials and methods section of our paper.

