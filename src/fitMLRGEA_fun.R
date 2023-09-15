#'##################################################################################################'
#'                                                                                                 #'
#' Project........Environmental GWAS (eGWAS) for PandAnd -- Hackathon                              #'
#' Title..........Multionomial Logistic Regression (MLR) for Genomic-Environment Association (GEA) #'
#' Function....... fitMLRGEA(), version 0.0.1 (Sep 15 2023)                                        #'
#' Created at..... 06-03-2023                                                                      #'
#' Updated at..... 09-15-2023                                                                      #'
#' Author: G.Costa-Neto, germano.cneto<at>gmail.com                                                #'
#' ....... Contributors: P. Bradbury (model testing, hypothesis)                                   #'
#'                                                                                                 #'
#'##################################################################################################'
#'
# packages

fitMLRGEA <- function(home.path = NULL,  # home directory. If null, getwd()
                      output.path = NULL,  # output directory. If null, getwd()
                      geno_Y = NULL, # data.frame for genomic features used as responses (taxa x features)
                      envo_X = NULL, # data.frame for environmental features used as predictors (taxa x features)
                      geno_X = NULL, # data.frame for genomic features used as covariates (taxa x features)
                      #   partial = 'full_model' # partial = c('full_model', 'envo_model','geno_model')
                      n.core = NULL, #  numeric value denoting the numbe of cores. 
                      output_name = NULL, # charcter denoting the name of the model under test
                      parallel = TRUE, # if is TRUE, and n.core is null, then use detectCores() - 1
                      verbose = FALSE #,save.meta=TRUE
)
{
  
  if(isTRUE(verbose))
  {
    message(paste0('========================================================================'))
    message(paste0('Multinomial Logistic Regresion for Genomic-Environment Associations (MLR-GEA)'))
    #  message(paste0('Author: G. Costa-Neto'))
    message(paste0('fitMLRGEA() function, v0.0.1 (Sep 15 2023)'))
    message(paste0('-========================================================================'))
    
  }
  
  #'-------------------------------------------------------------------------------------------------
  # paths and additional codes     #####
  #'-------------------------------------------------------------------------------------------------
  
  
  if(is.null(home.path)) home.path = getwd()
  if(is.null(   output.path))    output.path = getwd()
  if(is.null( output_name)) output_name <-'MLR_GEA_MODEL'
  
  source('https://raw.githubusercontent.com/gcostaneto/envFeatures/main/src/MakeMyLogRegEq.R')
  
  if (!requireNamespace("ordinal", quietly = TRUE)) {
    utils::install.packages("ordinal")
  }
  
  if (!requireNamespace("doParallel", quietly = TRUE)) {
    utils::install.packages("doParallel")
  }
  
  if (!requireNamespace("foreach", quietly = TRUE)) {
    utils::install.packages("foreach")
  }
  
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    utils::install.packages("reshape2")
  }
  
  
  #'-------------------------------------------------------------------------------------------------
  # Building the Models     #####
  #'-------------------------------------------------------------------------------------------------
  
  
  if(!is.null( geno_X))
  {
    common_rows <- 
      base::intersect(rownames(geno_X),rownames( envo_X)) %>% 
      base::intersect(rownames(  geno_Y)) %>%
      sort()
    
    geno_Y <- geno_Y[common_rows,]
    envo_X  <- envo_X[common_rows,]
    geno_X   <- geno_X[common_rows,]
    
    EnvTraits = colnames(  envo_X )
    GCovariates = colnames(  geno_X )
    
    model_inputs = 
      c(
        GCovariates,
        EnvTraits
      )
    
    form_null = MakeMyLogRegEq(Y =  'gFeature',X_tested =    '1' ) 
    form_full = MakeMyLogRegEq(Y =  'gFeature',X_tested =    model_inputs )
    form_geno = MakeMyLogRegEq(Y =  'gFeature',X_tested =    GCovariates ) 
    
  }
  
  if(is.null(geno_X))
  {
    common_rows <- 
      base::intersect(rownames(geno_Y),rownames( envo_X))
    
    geno_Y  <- geno_Y[common_rows,]
    envo_X  <- envo_X[common_rows,]
    geno_X  <- data.frame(dummt=matrix(NA,nrow = nrow(geno_Y)))
    
    EnvTraits = colnames(  envo_X )
    
    
    model_inputs = 
      c(
        EnvTraits
      )
    
    form_null = MakeMyLogRegEq(Y =  'gFeature',X_tested =    '1' ) 
    form_full = MakeMyLogRegEq(Y =  'gFeature',X_tested =    model_inputs )
    form_geno = MakeMyLogRegEq(Y =  'gFeature',X_tested =    '1' ) # geno_model = null_model if geno_X is null
    
  }
  
  
  n.geno.features  = ncol(geno_Y )
  
  #detectCores()
  
  if(is.null( n.core))
  {
    if(isTRUE(parallel))
    {
      cl <- makeCluster(detectCores() - 1)
      registerDoParallel(cl)
    }
  }
  if(!is.null(n.core))
  {
    cl <-  makeCluster(n.core)
    registerDoParallel(cl)
  }
  
  
  #'-------------------------------------------------------------------------------------------------
  # Running the models     #####
  #'-------------------------------------------------------------------------------------------------
  
  
  start = Sys.time()
  message(paste0('started at ......',start,'\n'))
  
  
  output_MLRGEA =
    try(
      foreach::foreach(geno.feature.j = 1:n.geno.features , 
                       .packages = c('plyr','reshape2','data.table','tidyverse','ordinal'),.combine = "rbind") %dopar%
        {
          
          #'---------------------------------------------------------
          # Model Step 1: organizing the models and inputs
          #'---------------------------------------------------------
          Feature_j = data.frame(Feature_j = geno_Y[,geno.feature.j])
          
          df_model = data.frame(cbind(data.frame(
            Feature_j = Feature_j$Feature_j  ,
            gFeature = as.factor(Feature_j$Feature_j )),envo_X, geno_X))
          
          
          if(length(unique(df_model$gFeature)) > 1) # we only run the logReg for Y with more than two classes
          {
            # running ordinal multinomial log regression
            
            null_model      <- try(ordinal::clm(formula = form_null,data = df_model))  # null model, H0
            full_model      <- try(ordinal::clm(formula = form_full,data = df_model)) # full model, H1
            geno_model      <- try(ordinal::clm(formula = form_geno,data = df_model))  # only Genetic Strucutres, H2
            #  envo_model      <- try(ordinal::clm(formula = form_ereduced,data = df_model))  # only Genetic Strucutres, H2
            
            
            if(!is.null(full_model))
            {
              
              table_PA =     table(df_model$gFeature)
              #'---------------------------------------------------------
              # Model Step 2: LRT  test and McFadden's pseudo R-sq
              #'---------------------------------------------------------
              
              H1 = as.data.frame(anova(   null_model,  full_model))
              H2 = as.data.frame(anova(  geno_model , full_model))
              H3 = as.data.frame(anova(  null_model, geno_model))
              #    H4 = as.data.frame(anova(  null_model, envo_model))
              
              #    as.data.frame(anova(  null_model, geno_model,full_model))
              #    as.data.frame(anova(  null_model, envo_model,  full_model))
              
              # general model statistics
              output.stat = data.frame(GenoFeature    = colnames(geno_Y)[geno.feature.j],
                                       n.parameters_full = length(model_inputs),
                                       
                                       # Model AIC
                                       AIC_null     =   H1$AIC[1],
                                       AIC_geno     =   H2$AIC[1],
                                       AIC_full     =   H1$AIC[2],
                                       
                                       # Model p-value
                                       log10p_H1 = -log10( H1$`Pr(>Chisq)`[2]),
                                       log10p_H2 = -log10( H2$`Pr(>Chisq)`[2]),
                                       log10p_H3 = -log10( H3$`Pr(>Chisq)`[2]),
                                       
                                       # McFadden's pseudo R-sq
                                       McF.pR2_H1   = NA,
                                       McF.pR2_H2   = NA,
                                       McF.pR2_H3   = NA)
              
              
              
              # McFadden's pseudo R-sq
              output.stat$ McF.pR2_H1   = try(100*round(1 - full_model$logLik/null_model$logLik,4))
              output.stat$ McF.pR2_H2   = try(100*round(1 - full_model$logLik/geno_model$logLik,4))
              output.stat$ McF.pR2_H3   = try(100*round(1 - geno_model$logLik/null_model$logLik,4))
              
              ## gcostaneto note: H1 = H2 + H3, that is, variability for a given feature = % explained by environment, % explained by not environment
              # some researchers would call this adaptation space = environmental space + genomic space
              
              #'---------------------------------------------------------
              #' Model Step 3: coef and partial p-values for each predictor
              #'---------------------------------------------------------
              # gcostaneto note: here I only pull the partial p and effects for the full model.
              # you can modify this code to pull for the env model, for instance. There is no limit for the imagination ;)
              
              output.partial       = summary(full_model)$coefficients
              output.partial       = data.frame( output.partial[which(rownames( output.partial) %in% model_inputs ),])
              
              output.partial <-  data.frame(Feature = rownames( output.partial), output.partial[,c(1:2,4)] )
              names( output.partial)[2:4] = c('effect','se_effect','p_value')
              output.partial$log10 = -log10( output.partial$p_value)
              
              
              # maximum log10p
              nulllog =  output.partial$log10[which( output.partial$log10 == Inf)] 
              if(!is.null( nulllog ))  output.partial$log10[which( output.partial$log10 == Inf)] = 323
              
              output.partial <- 
                output.partial %>% 
                reshape2::melt() %>%
                mutate(GenoFeature    = colnames(geno_Y)[geno.feature.j]) %>% 
                mutate(variable=paste0('partial_',variable)) %>% 
                reshape2::dcast(GenoFeature~variable+Feature,value.var = 'value')
              
              #'---------------------------------------------------------
              #' Model Step 4: Exporting ;)
              #'---------------------------------------------------------
              
              data.table::fwrite(       output.partial, paste0(output.path,'/model_partial',output_name ,".csv"),append = T)
              data.table::fwrite(       output.stat, paste0(output.path,'/model_statistics',output_name ,".csv"),append = T)
              
              output <-
                output.stat %>% 
                merge(output.partial ,by='GenoFeature')
            }
            
          }
          
          
          return(output)
        }
    )
  
  stopCluster(cl) 
  
  end = Sys.time()
  message(paste0('ended at ......',end,'\n'))
  
  message(paste0('saving files at......',output.path,'\n'))
  data.table::fwrite(  output_MLRGEA, paste0(output.path,'/Full_Results_',output_name ,".csv"))
  return(  output_MLRGEA)
}

