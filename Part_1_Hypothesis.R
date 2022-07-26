## Set the working directory
setwd("./")
#setwd("C:/Users/SaraNoeman/Documents/NU_MSC_Informatics/Courses/CIT660_Statistical_Analysis_Visualization/Project/CIT660_Final_Project/")

source("Hypothesis_Test_functions.R")

# Threshold for the Hypothesis test
alpha = 0.05

# Threshold of the FC = 2 ==> Then log2_FC_threshold = log2(2) = 1 
log2_FC_threshold = log2(2)

cancer_list = c("lusc","kirc")
#cancer_list = c("kirc")
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
 
  # 2-Data Cleansing
  
  ## Run the data cleansing function
  ## cleanse() function is responsible for the following:
  ## 1) Removing any Gene which has more than 50% of its samples with zero values (from both GE.cancer and GE.healthy)
  ## 2) For each Gene, check if there are patient samples with zero values, then remove this patient sample from both GE.cancer[Gene] and GE.healthy[Gene] for this Gene only.
  ## The function return a list containing GE_cleansed_data which is a list of two (GE.cancer.clean and GE.healthy.clean).
  GE_cleansed_data = cleanse(GE.cancer, GE.healthy)
  GE.cancer.clean = GE_cleansed_data[[1]]
  GE.healthy.clean = GE_cleansed_data[[2]]
  
  
  

  # 3-hypothesis testing on each gene, expected output (pvalues and adjusted pvalues)
  
  # Paired Test:
  Is_paired = TRUE
  Results_Paired_dataframe = Run_Hypothesis_Test(GE.cancer.clean, GE.healthy.clean, alpha, Is_paired)
  #sort Dataframe by Statistic
  Results_Paired_dataframe.sorted=Results_Paired_dataframe[ order( Results_Paired_dataframe$statisticValue, decreasing = TRUE),]
  #P adjusted moved to   Run_Hypothesis_Test function
  #pvalues_paired.adjusted = p.adjust(pvalues_paired, method = 'fdr')
  Top5_Genes.paired = Results_Paired_dataframe.sorted[1:5,]

  # Independent Test:
  Is_paired = FALSE
  Results_Independent_dataframe = Run_Hypothesis_Test(GE.cancer.clean, GE.healthy.clean, alpha, Is_paired)

  Results_Independent_dataframe.sorted=Results_Independent_dataframe[ order( Results_Independent_dataframe$statisticValue, decreasing = TRUE),]
  
  Top5_Genes.independent = Results_Independent_dataframe.sorted[1:5,]
  

  # 4-Fold change
  # useful resource: https://bioconductor.org/help/course-materials/2015/Uruguay2015/day5-data_analysis.html
  foldchange = list()
  FC = list()
  for (gene in  names(GE.cancer.clean) )
  {
    # Fold Change = ratio between the mean(GE_cancer) and mean(GE_healthy)
    # Fold Change = mean(GE_cancer)/mean(GE_healthy)
    # Then we need to get log2 of the Fold Change 
    # log2_FoldChange = log2(mean(GE_cancer)/mean(GE_healthy))
    #                 = log2(mean(GE_cancer)) - log2(mean(GE_healthy))
    
    # So for each gene, we calculate log2 of the mean of GE_cancer, and log2 of the mean of GE_healthy
    cancer.mean.log2 = log2(apply(data.frame(GE.cancer.clean[[gene]]), 1, mean))
    health.mean.log2 = log2(apply(data.frame(GE.healthy.clean[[gene]]), 1, mean))
    
    #cancer.mean = apply(data.frame(GE.cancer.clean[[gene]]), 1, mean)
    #health.mean = apply(data.frame(GE.healthy.clean[[gene]]), 1, mean)
    
    # we can take the difference between the log2 of means.  
    # And this is our log2 Fold Change == log2(condition/normal) = log2(condition) - log2(normal)
    # foldchange for this gene = cancer.mean.log2 - health.mean.log2
    foldchange[gene] = as.numeric(cancer.mean.log2-health.mean.log2)
    #FC[gene] = as.numeric(cancer.mean/health.mean)
  }

  # Draw hist for fold change
  # Note: if "Error in plot.new():figure margins too large" appears solve it by expanding the plotting window
  x_label = paste("log2 Fold Change (healthy  vs cancer)", cancer_type)
  log2_FoldChange = as.numeric(foldchange)
  
  ## The path of the png file that we will save the plot.
  path.name = paste(cancer_type, "_FoldChange_histogram.png")
  ## Define the png, with path.name to store it, and the height & width of the plot
  png(path.name, width = 2000, height = 1500, res=200) # To save the figure in the PNG format.
  plot_title = paste("Fold Change Histogram - ", cancer_type)
  hist(log2_FoldChange, xlab = x_label, main = plot_title)
  # To close the figure file and save it.
  dev.off()  
  #hist(xx, xlab = x_label, main = "Histogran of the Fold Change")
  print ("Part-1 DONE")
  


  # 5-Volcano plot
  # Paired adjusted p_value vs FoldChange:
  
  PValue.paired = -log10(Results_Paired_dataframe$pValue)
  PValue.paired.FDR = -log10(Results_Paired_dataframe$adjPValue)

  
  
  data.paired = data.frame(Genes=names(GE.cancer.clean), logFC=log2_FoldChange, PValue=PValue.paired, FDR=PValue.paired.FDR)
  
  ## The path of the png file that we will save the plot.
  #path.name = paste(cancer_type, "_Volcano_plot.png")
  #png(path.name, width = 2000, height = 1500, res=200) # To save the figure in the PNG format.
  
  ## Let your code determine the y-axis limits automatically 
  # y.min = floor(min(PValue.paired.FDR) - 1)
  # y.max = ceiling(max(PValue.paired.FDR) + 1)
  x.min = floor(min(log2_FoldChange))
  x.max = ceiling(max(log2_FoldChange))
  
  
  ###################################################################################
  # Create new categorical column ------------------------------------------------ 
  data.paired = dplyr::mutate(data.paired, gene_type = dplyr::case_when(log2_FoldChange >= log2_FC_threshold & PValue.paired.FDR >= -log10(alpha) ~ "up",
                                                              log2_FoldChange <= -log2_FC_threshold & PValue.paired.FDR >= -log10(alpha) ~ "down",
                                 TRUE ~ "ns"))   
  
  # Add colour, size and alpha (transparency) to volcano plot --------------------
  cols <- c("up" = "#ffad73", "down" = "#26b3ff", "ns" = "grey") 
  sizes <- c("up" = 2, "down" = 2, "ns" = 1) 
  alphas <- c("up" = 1, "down" = 1, "ns" = 0.5)
  
  title = paste("Volcano Plot for", cancer_type, "(Paired)" , sep = " ")
  
  p1 <- ggplot2::ggplot(data.paired, ggplot2::aes(x = log2_FoldChange,
               y = PValue.paired.FDR,
               fill = gene_type,    
               size = gene_type,
               alpha = gene_type)) + 
  ggplot2::geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
               colour = "black") + 
  ggplot2::geom_hline(yintercept = -log10(alpha),
               linetype = "dashed") + 
  ggplot2::geom_vline(xintercept = c(-log2_FC_threshold, log2_FC_threshold),
               linetype = "dashed") +
  ggplot2::scale_fill_manual(values = cols) + # Modify point colour
  ggplot2::scale_size_manual(values = sizes) + # Modify point size
  ggplot2::scale_alpha_manual(values = alphas) + # Modify point transparency
  ggplot2::scale_x_continuous(breaks = c(seq(-7, 7, 2)),       
                       limits = c(x.min, x.max)) +
  ggplot2::ggtitle(title) +
  ggplot2::xlab("log2(fold change)") + 
  ggplot2::ylab("-log10 (adj.p-value)")
  ###################################################################################
  volcano_plot.file = paste(cancer_type, "_paired_Volcano_plot.png", sep = "")
  ggplot2::ggsave(volcano_plot.file, plot = p1, 
                  device = png,
                  path = "./",
                  scale = 1,
                  units = c("in", "cm", "mm", "px"),
                  dpi = 300,
                  limitsize = TRUE)
  
  PValue.independent = -log10(Results_Independent_dataframe$pValue)
  PValue.independent.FDR = -log10(Results_Independent_dataframe$adjPValue)
  data.independent = data.frame(Genes=names(GE.cancer.clean), logFC=log2_FoldChange, PValue=PValue.independent, FDR=PValue.independent.FDR)

  ## Let your code determine the y-axis limits automatically 
  # y.min = floor(min(PValue.independent.FDR)-0.5)
  # y.max = ceiling(max(PValue.independent.FDR)+2)
  
  ###################################################################################
  # Create new categorical column ------------------------------------------------ 
  data.independent = dplyr::mutate(data.independent, gene_type = dplyr::case_when(log2_FoldChange >= log2_FC_threshold & PValue.independent.FDR >= -log10(alpha) ~ "up",
                                                                        log2_FoldChange <= -log2_FC_threshold & PValue.independent.FDR >= -log10(alpha) ~ "down",
                                                                        TRUE ~ "ns"))   
  
  title = paste("Volcano Plot for", cancer_type, "(Independent)" , sep = " ")
  p2 <- ggplot2::ggplot(data.independent, ggplot2::aes(x = log2_FoldChange,
                                            y = PValue.independent.FDR,
                                            fill = gene_type,    
                                            size = gene_type,
                                            alpha = gene_type)) + 
    ggplot2::geom_point(shape = 21, # Specify shape and colour as fixed local parameters    
                        colour = "black") + 
    ggplot2::geom_hline(yintercept = -log10(alpha),
                        linetype = "dashed") + 
    ggplot2::geom_vline(xintercept = c(-log2_FC_threshold, log2_FC_threshold),
                        linetype = "dashed") +
    ggplot2::scale_fill_manual(values = cols) + # Modify point colour
    ggplot2::scale_size_manual(values = sizes) + # Modify point size
    ggplot2::scale_alpha_manual(values = alphas) + # Modify point transparency
    ggplot2::scale_x_continuous(breaks = c(seq(-10, 10, 2)),       
                                limits = c(x.min, x.max)) +
    ggplot2::ggtitle(title) +
    ggplot2::xlab("log2(fold change)") + 
    ggplot2::ylab("-log10 (adj.p-value)")
  ###################################################################################
  
  volcano_plot.file = paste(cancer_type, "_independent_Volcano_plot.png", sep = "")
  ggplot2::ggsave(volcano_plot.file, plot = p2, 
                  device = png,
                  path = "./",
                  scale = 1,
                  units = c("in", "cm", "mm", "px"),
                  dpi = 300,
                  limitsize = TRUE)
  # Independent p_value vs FoldChange: 
  
  # Significant Genes based on Hypothesis Test 
  # a) Paired:
  # First Arrange the p-values in acsending order ==> to get the 
  # PValue.paired.FDR.sorted=sort(unlist(PValue.paired.FDR), decreasing=FALSE)
  

  # Report the set of DEGs in the above two pairing cases, and report how different these two sets of genes.
  # 1. Hypothesis testing,
  #    a) Paired Set:
  pvalues_paired.adjusted.sorted = Results_Paired_dataframe[ order( Results_Paired_dataframe$adjPValue, decreasing = FALSE),]
  DEGs.paired = Results_Paired_dataframe[(which(pvalues_paired.adjusted.sorted$adjPValue<=alpha)),1]
  
  #    b) Independent Set:
  pvalues_independent.adjusted.sorted = Results_Independent_dataframe[ order( Results_Independent_dataframe$adjPValue, decreasing = FALSE),]
  DEGs.independent = Results_Independent_dataframe[(which(pvalues_independent.adjusted.sorted$adjPValue<=alpha)),1]  
  
  # 2. Fold change:
  # Significant Genes based on Fold Change 
  DEGs.FC = names(which(c((foldchange>=log2_FC_threshold),(foldchange<=-log2_FC_threshold))))
  
  # 3. Volcano plot using the set of DEGs obtained by the hypothesis that data are paired
  # Significant Genes based on Volcano plot
  DEGs.volcano.paired = data.paired[which(data.paired$gene_type==c("up", "down")),]$Genes
  
  
  
  report.df=data.frame (nrow(GE.cancer), length(GE.cancer.clean),length( DEGs.independent), length(DEGs.paired), length(DEGs.FC), length(DEGs.volcano.paired) )
  names(report.df)=c("Total.Genes","Genes.After.Cleansing","DEGs.Independent","DEGS.Paired", "DEGS.LOG2FoldChange","DEGS.Volcano.paired")
  write.csv(report.df, file = paste(cancer_type,"_Report_DEGS.csv" ), row.names = FALSE)
  
  # 4. Top-5 significant Genes
  # Paired DEGs:
  # Sort the w_statistic of paired Hypothesis test ascendingly:
  # w_statistic_paired.sorted = sort(w_statistic_paired, decreasing = T, na.last = NA)

  # 6-Apply GSEA using the set of DEGs obtained by the hypothesis that data are paired
  # Prepare Files needed for GSEA (using Paired Data)
  Expression.file  = paste(cancer_type, "_Expression_dataset.paired.txt", sep = "")
  Phenotype.file = paste(cancer_type, "_Phenotype_dataset.paired.cls", sep = "")
  
  ## Concatenate the GE.Cancer data and GE.healthy data for the DEGs of paired hypothesis (DEGs.paired)
  GSEA.paired = as.data.frame(c(GE.cancer[DEGs.paired,], GE.healthy[DEGs.paired,-1]), row.names = row.names(GE.cancer[DEGs.paired,]))
  GSEA.cols = names(GSEA.paired)
  GSEA.cols[1]="Description"
  names(GSEA.paired) = GSEA.cols
  GSEA.paired[,1]="NA"
  cat("Name", GSEA.cols, "\n", sep = "\t", file = Expression.file)
  write.table(GSEA.paired, file = Expression.file, sep="\t", row.names=T, col.names = F, append = T)
  
  ## Prepare the .cls phenotype file
  cat(length(GSEA.paired)-1, "2 1\n", file = Phenotype.file, sep = " ")
  cat("# Cancer Healthy\n", file = Phenotype.file, sep = " ", append = T)
  Class_label = c(rep(1,length(GE.cancer[DEGs.paired,-1])), rep(0,length(GE.healthy[DEGs.paired,-1])))
  cat(Class_label, file = Phenotype.file, sep = " ", append = T)

  print("THE end of part 1")
  
  
  #Top5_Genes.paired$gene
  #Top5_Genes.independent$gene
  
  

}

