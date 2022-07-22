## Set the working directory
setwd("./")

## Read the Gene Expressions of Cancerous Tissues in a dataframe
GE_kirc.cancer <- read.table("./Project_Data/kirc-rsem-fpkm-tcga-t_paired.txt", header = T, row.names=1)
GE_kirc.healthy <- read.table("./Project_Data/kirc-rsem-fpkm-tcga_paired.txt", header = T, row.names=1)

GE_lusc.cancer <- read.table("./Project_Data/lusc-rsem-fpkm-tcga-t_paired.txt", header = T, row.names=1)
GE_lusc.healthy <- read.table("./Project_Data/lusc-rsem-fpkm-tcga_paired.txt", header = T, row.names=1)

#1- clean Kirc dataframes

## Remove the Genes that have more than 50% of its samples as zeros either in healthy or Cancerous tissues.
## We have 68 patients ==> 50% of patients = 34 ==> If number of zeros in GE level of Healthy or Cancerous data > =35 ==> remove this Gene
## skip the Gene index if the GE levels for this Gene in Cancerous tissues are more than 50% zeros ==> 35 or more patient
index_cancer = which(apply(GE_kirc.cancer == 0, 1, sum) >=35)
## skip the Gene index if the GE levels for this Gene in Healthy tissues are more than 50% zeros ==> 35 or more patient
index_healthy = which(apply(GE_kirc.healthy == 0, 1, sum) >=35)

## Here are the unique indices for Genes that need to be skipped
skip_index = unique(c(index_cancer, index_healthy))

## Remove Genes from both Healthy and Cancerous tissues GE levels.
GE_kirc.cancer.clean = GE_kirc.cancer[-skip_index,]
GE_kirc.healthy.clean = GE_kirc.healthy[-skip_index,]

#initialize clean gene list of key value pair 
GE_kirc.cancer.keyVal <- list()
GE_kirc.healthy.keyVal<- list()
for(Gene in rownames(GE_kirc.cancer.clean))
{
  GE_kirc.cancer.clean2 = GE_kirc.cancer.clean[Gene,]
  GE_kirc.healthy.clean2 = GE_kirc.healthy.clean[Gene,]
  ## index of patients with GE level of healthy tissues = 0
  patient_index1 = which(GE_kirc.healthy.clean[Gene,]==0)
  ## index of patients with GE level of Cancerous tissues = 0
  patient_index2 = which(GE_kirc.cancer.clean[Gene,]==0)
  
  skip_patient_index = unique(c(patient_index1, patient_index2))
  if (length(skip_patient_index)!=0)
  {
  GE_kirc.cancer.clean2 = GE_kirc.cancer.clean2[-skip_patient_index]
  GE_kirc.healthy.clean2 = GE_kirc.healthy.clean2[-skip_patient_index]
  }
  else
  {
    GE_kirc.cancer.clean2 = GE_kirc.cancer.clean2
    GE_kirc.healthy.clean2 = GE_kirc.healthy.clean2
  }
  
  
  #fill list of key value pair gene :samples
  GE_kirc.cancer.keyVal[[Gene]]<- GE_kirc.cancer.clean2
  GE_kirc.healthy.keyVal[[Gene]]<- GE_kirc.healthy.clean2
}

#2- clean Lusc dataframes 

## Remove the Genes that have more than 50% of its samples as zeros either in healthy or Cancerous tissues.
## We have 50 patients ==> 50% of patients = 25 ==> If number of zeros in GE level of Healthy or Cancerous data > =35 ==> remove this Gene
## skip the Gene index if the GE levels for this Gene in Cancerous tissues are more than 50% zeros ==> 35 or more patient
index_cancer = which(apply(GE_lusc.cancer == 0, 1, sum) >=25)
## skip the Gene index if the GE levels for this Gene in Healthy tissues are more than 50% zeros ==> 35 or more patient
index_healthy = which(apply(GE_lusc.healthy == 0, 1, sum) >=25)

## Here are the unique indices for Genes that need to be skipped
skip_index = unique(c(index_cancer, index_healthy))

## Remove Genes from both Healthy and Cancerous tissues GE levels.
GE_lusc.cancer.clean = GE_lusc.cancer[-skip_index,]
GE_lusc.healthy.clean = GE_lusc.healthy[-skip_index,]

#initialize clean gene list of key value pair 
GE_lusc.cancer.keyVal <- list()
GE_lusc.healthy.keyVal<- list()
for(Gene in rownames(GE_lusc.cancer.clean))
{
  GE_lusc.cancer.clean2 = GE_lusc.cancer.clean[Gene,]
  GE_lusc.healthy.clean2 = GE_lusc.healthy.clean[Gene,]
  ## index of patients with GE level of healthy tissues = 0
  patient_index1 = which(GE_lusc.healthy.clean[Gene,]==0)
  ## index of patients with GE level of Cancerous tissues = 0
  patient_index2 = which(GE_lusc.cancer.clean[Gene,]==0)
  
  skip_patient_index = unique(c(patient_index1, patient_index2))
  if (length(skip_patient_index)!=0)
   { 
    GE_lusc.cancer.clean2 = GE_lusc.cancer.clean2[-skip_patient_index]
    GE_lusc.healthy.clean2 = GE_lusc.healthy.clean2[-skip_patient_index]
    }
  else
    {
      GE_lusc.cancer.clean2 = GE_lusc.cancer.clean2
      GE_lusc.healthy.clean2 = GE_lusc.healthy.clean2
    }
  #fill list of key value pair gene :samples
  GE_lusc.cancer.keyVal[[Gene]]<-GE_lusc.cancer.clean2
  GE_lusc.healthy.keyVal[[Gene]] <-GE_lusc.healthy.clean2
}

#3-hypothesis testing on each gene, expected output (pvalues and adjusted pvalues)

#4-Fold change

#5-Volcano plot

#6-GSEA

