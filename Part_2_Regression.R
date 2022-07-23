extract_gene_CNA_df <- function (gene_cancer_profile, CNA_data){
  
  # replace the '.' in the gene_cancer_profile with '-' to match the CNA_data 
  names(gene_cancer_profile) =  gsub("[.]" , "-" , names(gene_cancer_profile))
  
  
  # save the gene profile samples and the CNA data samples
  CNA_data_samples = row.names(CNA_data)
  gene_profile_samples =  names(gene_cancer_profile)
  
  
  # find the intersection between the GE profiles and the CNAs profiles
  intersect_samples = intersect(gene_profile_samples, CNA_data_samples)
  # save the intersect in a variable
  intersect_CNA_data = CNA_data[intersect_samples,]
  
  
  # transpose the gene_cancer_profile to merge it correctly with the CNA_data
  gene_cancer_profile = t(gene_cancer_profile)
  # merge the two gene df with the CNA df by the row.names (samples names)
  gene_CNA_df = merge(x = gene_cancer_profile, y= intersect_CNA_data,
                   by.x = 0, by.y = 0,)
  # view the gene_cancer_profile after transposing
  View(gene_cancer_profile)
  
  
  # make the index of the dataframe is the Row.names which is the samples names
  rownames(gene_CNA_df) <- gene_CNA_df$Row.names
  # delete the samples names column from the dataframe after making it the index
  gene_CNA_df$Row.names <- NULL
  
  
  # return a list of objects required from this function
  return_list <- list("df" = gene_CNA_df , "y" = colnames(gene_cancer_profile))
  return(return_list)
}

