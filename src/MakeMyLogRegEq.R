# Code to make formula for logReg model
# updated at April 2023

MakeMyLogRegEq = function(Y,X_tested=NULL, X_X_covariates=NULL)
{
  Y = Y
  X = 0
  
  if(!is.null(X_tested))
  {
    X = X_tested[1]
    
    if(length(X_tested) > 1)
    {
      for(i in 2:length(X_tested)) X = paste0(X,'+',X_tested[i])
    }
    
  }
  
  
  if(!is.null(X_covariates))
  {
    if(is.null(X_tested))
    {
      X = X_covariates[1]
    }
    
    if(!is.null(X_tested))
    {
      X = paste0(X,'+',X_covariates[1])
    }
    
    if(length(X_tested) > 1)
    {
      for(i in 1:length(X_tested)) X = paste0(X,'+',X_tested[i])
    }
  }
  
  myFormula = formula(paste0(Y,'~',X))
  return(myFormula)
}
