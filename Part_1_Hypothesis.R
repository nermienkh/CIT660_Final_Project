## Set the working directory
setwd("C:/Users/SaraNoeman/Documents/NU_MSC_Informatics/Courses/CIT660_Statistical_Analysis_Visualization/Project/CIT660_Final_Project/")

source("Cleanse_function.R")

cancer_list = c("lusc")

for (cancer_type in cancer_list)
{
  
  print(cancer_type)
  GE_cancer.filename = paste("./Project_Data/", cancer_type, "-rsem-fpkm-tcga-t_paired.txt", sep = "")
  GE_healthy.filename = paste("./Project_Data/", cancer_type, "-rsem-fpkm-tcga_paired.txt", sep = "")
  
  CNA.filename = paste("./Project_Data/", cancer_type, "_CNV_core.txt", sep = "")
  
  ## Read the Gene Expressions of Cancerous and Healthy Tissues in a dataframe
  GE.cancer <- read.table(GE_cancer.filename, header = T, row.names=1)
  GE.healthy <- read.table(GE_healthy.filename, header = T, row.names=1)
  
  
  ## Read the Copy Number Alterations (CNAs) Cancerous Samples in a dataframe
  CNA <- read.table(CNA.filename, header = T, row.names=1)
  
  ## Run the data cleansing function
  ## cleanse() function is responsible for the following:
  ## 1) Removing any Gene who has more than 50% of its samples with zero values (from both GE.cancer and GE.healthy)
  ## 2) For each Gene, check if there are patient samples with zero values, then remove this patient sample from both GE.cancer and GE.healthy for this Gene only.
  ## The function return a list containing GE.cancer.clean and GE.healthy.clean.
  GE_cleansed_data = cleanse(GE.cancer, GE.healthy)
  GE.cancer.clean = GE_cleansed_data[[1]]
  GE.healthy.clean = GE_cleansed_data[[2]]
  
  
  

  # 3-hypothesis testing on each gene, expected output (pvalues and adjusted pvalues)
  
  # Paired Test:
  Is_paired = TRUE
  pvalues_paired = Run_Hypothesis_Test(GE.cancer.clean, GE.healthy.clean, 0.05, Is_paired)
  pvalues_paired.adjusted = p.adjust(pvalues_paired, method = 'bonferroni')
  # Top5_Genes.paired = names(pvalues_paired.sorted[c(1:5)])

  # Independent Test:
  Is_paired = FALSE
  pvalues_independent = Run_Hypothesis_Test(GE.cancer.clean, GE.healthy.clean, 0.05, Is_paired)
  pvalues_independent.adjusted = p.adjust(pvalues_independent, method = 'bonferroni')
  # Top5_Genes.independent = names(p_values_independent.sorted[c(1:5)])
  

  # 4-Fold change
  # useful resource: https://bioconductor.org/help/course-materials/2015/Uruguay2015/day5-data_analysis.html
  # foldchange = data.frame(row.names = rownames(GE.cancer.clean))
  foldchange = c()
  for (gene in  names(GE.cancer.clean) )
  {
    ##transform data into log2 base.
    # GE.cancer.log = log2(GE.cancer.clean[[gene]])
    # GE.healthy.log = log2(GE.healthy.clean[[gene]])
    ##calculate the mean for each gene row per group
    cancer = log2(apply(data.frame(GE.cancer.clean[[gene]]), 1, mean))
    healthy = log2(apply(data.frame(GE.healthy.clean[[gene]]), 1, mean))

    ##we can take the difference between the means.  And this is our log2 Fold Change or log2 Ratio == log2(control / test)
    # foldchange[gene, "Mean_log2"] = as.numeric(healthy - cancer)
    foldchange = append(foldchange, cancer-healthy)
    #print(foldchange[gene, 1])
  }

  #draw hist for fold change
  ##note: if "Error in plot.new():figure margins too large" appears solve it by expanding the plotting window
  x_label = paste("log2 Fold Change (healthy  vs cancer)", cancer_type)
  log2_FoldChange = foldchange
  hist(log2_FoldChange, xlab = x_label, main = "Histogran of the Fold Change")


  #5-Volcano plot
  # Paired p_value vs FoldChange:
  
  PValue.paired = -log10(as.numeric(pvalues_paired))
  PValue.paired.FDR = -log10(as.numeric(pvalues_paired.adjusted))
  PValue.independent = -log10(as.numeric(pvalues_independent))
  PValue.independent.FDR = -log10(as.numeric(pvalues_independent.adjusted))
  
  
  data.paired = data.frame(Genes=names(GE.cancer.clean), logFC=log2_FoldChange, PValue=PValue.paired, FDR=PValue.paired.FDR)
                 
  data.independent = data.frame(Genes=names(GE.cancer.clean), logFC=log2_FoldChange, PValue=PValue.independent, FDR=PValue.independent.FDR)
  
  ggplot2::ggplot(data=data.independent, ggplot2::aes(x=log2_FoldChange, y=PValue.independent.FDR)) + ggplot2::geom_point(size = 1/5, na.rm = TRUE) + ggplot2::xlim(-4, 4) + ggplot2::ylim(0,20)
  
  ggplot2::ggplot(data.paired, ggplot2::aes(log2_FoldChange, PValue.paired.FDR)) + # -log10 conversion  
    ggplot2::geom_point(size = 2/5, na.rm = TRUE) +
    ggplot2::xlab("log2_FC") + 
    ggplot2::ylab("-log10_FDR")  + ggplot2::xlim(-5, 5) + ggplot2::ylim(0, 9) 
  
  # Independent p_value vs FoldChange: 
#  log10_pvalue.independent = -log10(p_values_independent.sorted)
  

  #6-GSEA


  
}


