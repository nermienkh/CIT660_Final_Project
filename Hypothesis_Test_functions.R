
################################################################################################
##            Is_Normal_Distribution() function checks if all the different Gene data         ##
##          for Cancerous and Healthy tissues are following normal distribution or not        ##
##            If the dataset of any Gene for either Healthy or Cancerous tissues              ##
##     is not Normally distributed, then the function returns FALSE, else it returns TRUE     ##
################################################################################################
Is_Normal_Distribution<-function(pairedFlag,listOfKeyValCancerous,listOfKeyValHealthy, alpha)
{
  if(!pairedFlag) 
  {
    for (gene in  names(listOfKeyValCancerous) )
    {
      group1shapiro = shapiro.test(unlist(listOfKeyValCancerous[[gene]]))
      group2shapiro = shapiro.test(unlist(listOfKeyValHealthy[[gene]]))
      
      ## p_value of group1shapiro is < alpha ==> Not a normal distribution ==> We can not use t-test ==> Go with Wilcoxon
      if (group1shapiro$p.value<alpha)
      {
        return (FALSE)
      }
      ## p_value of group2shapiro is < alpha ==> Not a normal distribution ==> We can not use t-test ==> Go with Wilcoxon
      else if (group2shapiro$p.value<alpha)
      {
        return (FALSE)
      }
    }
    
  }
  else #paired
  {
    difference=c()
    ## p_value of difference between Cancerous and Healthy data is < alpha ==> Not a normal distribution ==> We can not use t-test ==> Go with Wilcoxon
    for (gene in  names(listOfKeyValCancerous) )
    {
      difference=unlist(listOfKeyValCancerous[[gene]]) - unlist(listOfKeyValHealthy[[gene]])
      
      diff.shappiro=shapiro.test(difference)
      
      if (diff.shappiro$p.value<alpha)
      {
        return (FALSE)
      }
    }
  }
  
  return (TRUE)
}
################################################################################################



################################################################################################
##                    cleanse() function is responsible for the following                     ##
## 1) Removing any Gene who has more than 50% of its samples with zero values                 ## 
##    (from both GE.cancer and GE.healthy)                                                    ##
## 2) For each Gene, check if there are patient samples with zero values,                     ##
##    then remove this patient sample from both GE.cancer and GE.healthy for this Gene only.  ##
## The function return a list containing GE.cancer.clean and GE.healthy.clean.                ##
################################################################################################
cleanse <- function(GE.cancer, GE.healthy)
{
  #1- clean the dataframes
  ##Remove the first column from data farmes
  GE.cancer =GE.cancer[,-1]
  GE.healthy=GE.healthy[,-1]
  
  ## Remove the Genes that have more than 50% of its samples as zeros either in healthy or Cancerous tissues.
  ## We have 68 patients ==> 50% of patients = 34 ==> If number of zeros in GE level of Healthy or Cancerous data > =34 ==> remove this Gene
  ## skip the Gene index if the GE levels for this Gene in Cancerous tissues are more than 50% zeros ==> 34 or more patient
  max_zeros_threshold = length(colnames(GE.cancer))/2
  index_cancer = which(apply(GE.cancer == 0, 1, sum) >=max_zeros_threshold)
  ## skip the Gene index if the GE levels for this Gene in Healthy tissues are more than 50% zeros ==> 34 or more patient
  index_healthy = which(apply(GE.healthy == 0, 1, sum) >=max_zeros_threshold)
  
  ## Here are the unique indices for Genes that need to be skipped
  skip_index = unique(c(index_cancer, index_healthy))
  
  ## Remove Genes from both Healthy and Cancerous tissues GE levels. 
  GE.cancer.clean = GE.cancer[-skip_index,]
  GE.healthy.clean = GE.healthy[-skip_index,]
  
  # #initialize clean gene list of key value pair 
  GE.cancer.keyVal <- list()
  GE.healthy.keyVal<- list()
  for(Gene in rownames(GE.cancer.clean))
  {
    GE.cancer.clean2 = GE.cancer.clean[Gene,]
    GE.healthy.clean2 = GE.healthy.clean[Gene,]
    ## index of patients with GE level of healthy tissues = 0
    patient_index1 = which(GE.healthy.clean[Gene,]==0)
    ## index of patients with GE level of Cancerous tissues = 0
    patient_index2 = which(GE.cancer.clean[Gene,]==0)
    
    #combine patient_index1, patient_index2 to skip them
    skip_patient_index = unique(c(patient_index1, patient_index2))

    if (length(skip_patient_index)!=0)
    {
      # GE.cancer.clean[Gene,skip_patient_index] = 0
      # GE.healthy.clean[Gene,skip_patient_index] = 0
      GE.cancer.clean2 = GE.cancer.clean2[-skip_patient_index]
      GE.healthy.clean2 = GE.healthy.clean2[-skip_patient_index]
    }
    
    #fill list of key value pair gene :samples
    GE.cancer.keyVal[[Gene]]<- GE.cancer.clean2
    GE.healthy.keyVal[[Gene]]<- GE.healthy.clean2
  }
  return (list(GE.cancer.keyVal, GE.healthy.keyVal))
}



#3-hypothesis testing on each gene, expected output (pvalues and adjusted pvalues)
Run_Hypothesis_Test <- function(GE.cancer.keyVal, GE.healthy.keyVal, alpha, Is_paired)
{
  #A- Paired  Test 
  normalityCheck_GE=Is_Normal_Distribution(Is_paired ,GE.cancer.keyVal, GE.healthy.keyVal, alpha)
  ## we  will use Wilcoxon signed rank test. Paired = True
  
  pvalues=list()
  ## If All Gene data (Cancerous and Healthy) are Normally Distributed
  if(normalityCheck_GE == TRUE)
  {
    for (gene in  names(GE.cancer.keyVal))
    {
      ## ==> Apply t-test
      result= t.test( unlist(GE.cancer.keyVal[[gene]]),unlist(GE.healthy.keyVal[[gene]]), alternative = 'two.sided', paired = Is_paired)
      pvalues[gene]=result$p.value
    }
  }
  ## If any of the Gene data (Cancerous and Healthy for independent OR the difference for paired) are Normally Distributed
  else
  {
    for (gene in  names(GE.cancer.keyVal))
    {
      ## ==> Apply Wilcox-test
      paste("GE.cancer.keyVal[[",gene,"]] = ", length(GE.cancer.keyVal[[gene]]), "\n")
      paste("GE.healthy.keyVal[[",gene,"]] = ", length(GE.healthy.keyVal[[gene]]), "\n\n")
      
      result= wilcox.test( unlist(GE.cancer.keyVal[[gene]]),unlist(GE.healthy.keyVal[[gene]]), alternative = 'two.sided', paired = Is_paired)
      pvalues[gene]=result$p.value
    }
  }
  #sorting pvalues  
  #pvalues.sorted=sort(unlist(pvalues), decreasing=FALSE)
  return (pvalues)
}