CheckNormality50PercentOfGroups<-function(pairedFlag,listOfKeyValCancerous,listOfKeyValHealthy, alpha)
{
  count=0
  if(!pairedFlag) 
  {
    for (gene in  names(listOfKeyValCancerous) )
  {
    group1shapiro= shapiro.test(unlist(listOfKeyValCancerous[[gene]]))
    group2shapiro=shapiro.test(unlist(listOfKeyValHealthy[[gene]]))
    
    if (group1shapiro$p.value>alpha)
    {
      count=count+1
    }
    else if (group2shapiro$p.value>alpha)
    {
      count=count+1
    }
  }
  
  }
  else #paired
  {
    difference=c()
    for (gene in  names(listOfKeyValCancerous) )
    {
      difference=unlist(listOfKeyValCancerous[[gene]]) - unlist(listOfKeyValHealthy[[gene]])
      
      diff.shappiro=shapiro.test(difference)
      
      if (diff.shappiro$p.value>alpha)
      {
        count=count+1
        
      }
    }
  }
  
  #check if %50 of genes  are normal or not
  if (count>=length(listOfKeyValCancerous)*0.5)
  {
    return (TRUE);
  }  
  else
  {
    return (FALSE)
  }
}