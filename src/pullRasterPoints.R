#'==================================================================================================
#' Title.    : Pulling Point information from digital image files -- raster
#' Author.   : G Costa-Neto
#' Created at: 2023-01-15
#' Updated at: 2023-04-12
#' 
#' Previous versions: EnvRtype::extract_GIS()
#' Based on terra:extract()
#' 
#' pullRasterPoints()
#' Package: envFeatures (in developing, current version 0.0.1,"beta" April 2023)
#'==================================================================================================



pullRasterPoints <- function(digital.raster =NULL,
                             which.raster.number = NULL, 
                             lat =NULL, lng =NULL, env.dataframe=NULL,
                             env.id=NULL,name.covariate = NULL, merge=TRUE){
  
  
  message(paste0('========================================================================'))
  message(paste0('EnvFeatures: democratizing the environmental data collection'))
  message(paste0('Author: G. Costa-Neto'))
  message(paste0('pullRasterPoints() function, v0.0.1'))
  message(paste0('-========================================================================'))
  
  
  if (!requireNamespace("terra", quietly = TRUE)) {
    utils::install.packages("terra")
  }
  if (!requireNamespace("raster", quietly = TRUE)) {
    utils::install.packages("raster")
  }
  
  if (!requireNamespace("sp", quietly = TRUE)) {
    utils::install.packages("sp")
  }
  
  if(length(digital.raster) == 1 )  digital.raster <- digital.raster
  if(!is.null( which.raster.number )) 
  {
    digital.raster <- digital.raster[[ which.raster.number ]]
  }
  
  if(!is.null(name.covariate))  names(digital.raster) = name.covariate
  if(is.null(name.covariate)) name.covariate = names(digital.raster)
  if(is.null(lat)) message('provide the latitude column name for lat in env.dataframe')
  if(is.null(lng)) message('provide the longitude column name for lng in env.dataframe')
  if(is.null(env.id)) message('provide the environmental id column name for env.id in env.dataframe')
  
  if(is.numeric(lat) & is.numeric(lng))
  {
    coords    <- data.frame(x = lng, y = lat)
  }
  else{
    coords    <- data.frame(x = env.dataframe[,lng], y = env.dataframe[,lat])
  }
  
  sp::proj4string(  coords) = sp::CRS("+proj=longlat +datum=WGS84") 
  sp_vector <- terra::vect(sp::SpatialPoints(coords))
  
  if(isFALSE(class(digital.raster)[1] == "SpatRaster")) digital.raster = terra::rast(digital.raster)
  extracted_values <- terra::extract(digital.raster , sp_vector)
  
  
  if(!is.null(env.id))
  {
    if(is.null(env.dataframe))
    {
      extracted_values$ID = env.id
    }
    
    if(!is.null(env.dataframe))
    {
      extracted_values$ID = env.dataframe[,env.id]
      names(extracted_values)[1] = env.id
      
      if(isTRUE(merge))
      {
        extracted_values = merge(extracted_values,env.dataframe,by=env.id)
      }
    }
  }
  
  return(extracted_values)
}
