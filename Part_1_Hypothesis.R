## Set the working directory
setwd("./")

## Read the Gene Expressions of Cancerous and Healthy Tissues in a dataframe
GE_kirc.cancer <- read.table("./Project_Data/kirc-rsem-fpkm-tcga-t_paired.txt", header = T, row.names=1)
GE_kirc.healthy <- read.table("./Project_Data/kirc-rsem-fpkm-tcga_paired.txt", header = T, row.names=1)

GE_lusc.cancer <- read.table("./Project_Data/lusc-rsem-fpkm-tcga-t_paired.txt", header = T, row.names=1)
GE_lusc.healthy <- read.table("./Project_Data/lusc-rsem-fpkm-tcga_paired.txt", header = T, row.names=1)

## Read the Copy Number Alterations (CNAs) Cancerous Samples in a dataframe
kirc_CNV <- read.table("./Project_Data/kirc_CNV_core.txt", header = T, row.names=1)
lusc_CNV <- read.table("./Project_Data/lusc_CNV_core.txt", header = T, row.names=1)


#1- clean Kirc dataframes
##Remove the first column from data farmes
GE_kirc.cancer =GE_kirc.cancer [,-1]
GE_kirc.healthy=GE_kirc.healthy[,-1]

## Remove the Genes that have more than 50% of its samples as zeros either in healthy or Cancerous tissues.
## We have 68 patients ==> 50% of patients = 34 ==> If number of zeros in GE level of Healthy or Cancerous data > =34 ==> remove this Gene
## skip the Gene index if the GE levels for this Gene in Cancerous tissues are more than 50% zeros ==> 34 or more patient
index_cancer = which(apply(GE_kirc.cancer == 0, 1, sum) >=34 )
## skip the Gene index if the GE levels for this Gene in Healthy tissues are more than 50% zeros ==> 34 or more patient
index_healthy = which(apply(GE_kirc.healthy == 0, 1, sum) >=34)

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
  
  #combine patient_index1, patient_index2 to skip them
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
##Remove the first column from data farmes
GE_lusc.cancer = GE_lusc.cancer[,-1]
GE_lusc.healthy=GE_lusc.healthy[,-1]

## Remove the Genes that have more than 50% of its samples as zeros either in healthy or Cancerous tissues.
## We have 50 patients ==> 50% of patients = 25 ==> If number of zeros in GE level of Healthy or Cancerous data > =25 ==> remove this Gene
## skip the Gene index if the GE levels for this Gene in Cancerous tissues are more than 50% zeros ==> 25 or more patient
index_cancer = which(apply(GE_lusc.cancer == 0, 1, sum) >=25)
## skip the Gene index if the GE levels for this Gene in Healthy tissues are more than 50% zeros ==> 25 or more patient
index_healthy = which(apply(GE_lusc.healthy == 0, 1, sum) >=25)

## Here are the unique indices for Genes that need to be skipped
skip_index = unique(c(index_cancer, index_healthy))

## Remove Genes from both Healthy and Cancerous tissues GE levels.
GE_lusc.cancer.clean = GE_lusc.cancer[-skip_index,]
GE_lusc.healthy.clean = GE_lusc.healthy[-skip_index,]

# Initialize clean gene list of key value pair 
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
  
  #combine patient_index1, patient_index2 to skip them
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
source("Functions.R")
#A- Paired  Test 
normalityCheck_GE_kirc_Paired=CheckNormality50PercentOfGroups(TRUE ,GE_kirc.cancer.keyVal, GE_kirc.healthy.keyVal,0.05)
normalityCheck_GE_lusc_Paired=CheckNormality50PercentOfGroups(TRUE ,GE_lusc.cancer.keyVal, GE_lusc.healthy.keyVal,0.05)
## since  more than 50% of data are not normally distributed, we  will use Wilcoxon signed rank test. Paired = True
#A kirc test
kirc_Paired_pvalues=list()
for (gene in  names(GE_kirc.cancer.keyVal))
{
 result= wilcox.test( unlist(GE_kirc.cancer.keyVal[[gene]]),unlist(GE_kirc.healthy.keyVal[[gene]]), alternative = 'two.sided', paired = TRUE)
 kirc_Paired_pvalues[gene]=result$p.value
}

#A lusc test
lusc_Paired_pvalues=list()
for (gene in  names(GE_lusc.cancer.keyVal))
{
  result= wilcox.test(unlist(GE_lusc.cancer.keyVal[[gene]]),unlist(GE_lusc.healthy.keyVal[[gene]]), alternative = 'two.sided', paired = TRUE)
  lusc_Paired_pvalues[gene]=result$p.value
}
#sorting pvalues  
kirc_Paired_pvalues.sorted=sort(unlist(lusc_Paired_pvalues), decreasing=FALSE)
lusc_Paired_pvalues.sorted=sort(unlist(lusc_Paired_pvalues), decreasing=FALSE)

#B- independent Test
normalityCheck_GE_kirc_independent=CheckNormality50PercentOfGroups(FALSE ,GE_kirc.cancer.keyVal, GE_kirc.healthy.keyVal,0.05)
normalityCheck_GE_lusc_independent=CheckNormality50PercentOfGroups(FALSE,GE_lusc.cancer.keyVal, GE_lusc.healthy.keyVal,0.05)
#since  more than 50% of data are not normally distributed, we  will us Wilcoxon rank sum test Paired=False
#B kirc test
kirc_independent_pvalues=list()
for (gene in  names(GE_kirc.cancer.keyVal))
{
  result= wilcox.test( unlist(GE_kirc.cancer.keyVal[[gene]]),unlist(GE_kirc.healthy.keyVal[[gene]]), alternative = 'two.sided')
  kirc_independent_pvalues[gene]=result$p.value
}

#B lusc test
lusc_independent_pvalues=list()
for (gene in  names(GE_lusc.cancer.keyVal))
{
  result= wilcox.test(unlist(GE_lusc.cancer.keyVal[[gene]]),unlist(GE_lusc.healthy.keyVal[[gene]]), alternative = 'two.sided')
  lusc_independent_pvalues[gene]=result$p.value
}
  
#sorting pvalues  
kirc_independent_pvalues.sorted=sort(unlist(lusc_independent_pvalues), decreasing=FALSE)
lusc_independent_pvalues.sorted=sort(unlist(lusc_independent_pvalues), decreasing=FALSE)




#4-Fold change
# useful resource: https://bioconductor.org/help/course-materials/2015/Uruguay2015/day5-data_analysis.html
#A- lusc Data
GE_lusc.cancer.keyVal.log=list()
GE_lusc.healthy.keyVal.log=list()
lusc_foldchange=list()
for (gene in  names(GE_lusc.cancer.keyVal) )
{
  ##transform data into log2 base.
  GE_lusc.cancer.keyVal.log= log2(GE_lusc.cancer.keyVal[[gene]])
  GE_lusc.healthy.keyVal.log= log2(GE_lusc.healthy.keyVal[[gene]])
  ##calculate the mean for each gene row per group
  cancer = apply(GE_lusc.cancer.keyVal.log, 1, mean)
  healthy = apply(GE_lusc.healthy.keyVal.log, 1, mean)
  ##we can take the difference between the means.  And this is our log2 Fold Change or log2 Ratio == log2(control / test)
  lusc_foldchange<- append(lusc_foldchange, as.numeric(healthy - cancer))
}
#B- kirc Data
GE_kirc.cancer.keyVal.log=list()
GE_kirc.healthy.keyVal.log=list()
kirc_foldchange=list()
for (gene in  names( GE_kirc.healthy.keyVal) )
{
  ##transform data into log2 base.
  GE_kirc.cancer.keyVal.log= log2( GE_kirc.cancer.keyVal[[gene]])
  GE_kirc.healthy.keyVal.log= log2( GE_kirc.healthy.keyVal[[gene]])
  ##calculate the mean for each gene row per group
  cancer = apply( GE_kirc.cancer.keyVal.log, 1, mean)
  healthy = apply( GE_kirc.healthy.keyVal.log, 1, mean)
  ##we can take the difference between the means.  And this is our log2 Fold Change or log2 Ratio == log2(control / test)
  kirc_foldchange<- append(kirc_foldchange, as.numeric(healthy - cancer))
}


#draw hist for fold change
##note: if "Error in plot.new():figure margins too large" appears solve it by expanding the plotting window
hist(unlist(lusc_foldchange), xlab = "log2 Fold Change (healthy  vs cancer) lusc")
hist(unlist(kirc_foldchange), xlab = "log2 Fold Change (healthy  vs cancer) kirc")




#5-Volcano plot

#6-GSEA

