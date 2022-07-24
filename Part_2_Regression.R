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

gene_CNA_regression <- function(gene_CNA_df , gene_name ){
  
  # extract the CNA columns as a matrix
  x = data.matrix(gene_CNA$df[,-1])

  # get the vectors of y
  y = gene_CNA_df[[gene_name]]
  
  # use a feature selection technique to penalize features: using LASSO
  library(glmnet)
  # get the lambda by the cross validation technique
  fit.cv <- cv.glmnet(x, y , family = "gaussian" , alpha = 1, standardize = FALSE , nfolds = 10 )
  lambda <- fit.cv$lambda.min
  cat("lambda = " , lambda , "
      ")
  
  # change lambda because generated lambda is not working well
  lambda = 0.05
  
  # calculate the linear regression relation between y and x1
  model <- glmnet(x,y, family = "gaussian", alpha =1 , lambda = lambda , standardize = FALSE)
  
  # print the linear regression relation summary
  print(summary(model))
  
  
  coef.fit <- coef(model, s = lambda)[2:114]
  print(coef.fit)
  
  features.in <- which(abs(coef.fit) > 0 )
  print(features.in)
  
  x = x[, features.in]
  
  relation <- lm(y~x)
  print(summary(relation))
  
  # plot the vectors points with the regression line
  # plot(x,y , col="blue", main = "Gene and CNA regression", abline(model),
       # cex = 1.3,pch = 16,xlab = "CNA",ylab = "Gene Expression Level")
  
}