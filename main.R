## Set the working directory
setwd("./")

## Read the Gene Expressions of Cancerous and Healthy Tissues in a dataframe
GE_kirc.cancer <- read.table("./Project_Data/kirc-rsem-fpkm-tcga-t_paired.txt", header = T, row.names=1)
GE_kirc.healthy <- read.table("./Project_Data/kirc-rsem-fpkm-tcga_paired.txt", header = T, row.names=1)

## Read the Copy Number Alterations (CNAs) Cancerous Samples in a dataframe
kirc_CNV <- read.table("./Project_Data/kirc_CNV_core.txt", header = T, row.names=1)
lusc_CNV <- read.table("./Project_Data/lusc_CNV_core.txt", header = T, row.names=1)


# ======== will need a for loop for the most 5 genes 

gene_profile_example = GE_lusc.cancer.keyVal[["TOMM70A"]]

gene_CNA = extract_gene_CNA_df(gene_profile_example, lusc_CNV)

View(gene_CNA$df)
printint(gene_CNA$y)

gene_CNA_regression(gene_CNA$df ,gene_CNA$y )





