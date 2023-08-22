#' Create the settings for Boruta feature selection
#'
#' @details
#' From: https://doi.org/10.18637/jss.v036.i11
#'
#' @param nJobs       How many jobs to do in parallel
#' @param maxDepth    Max depth of each tree in the RandomForest used
#' @param nTrees      How many trees to use, default is `auto`
#'
#' @return
#' An object of class \code{featureEngineeringSettings}
#' @export
createBorutaFeatureSelection <- function(nJobs = 10L, 
                                         maxDepth = 5L,
                                         nTrees = "auto",
                                         verbosity = 2L,
                                         iterations = 100L,
                                         randomState = 42L
                                         ){
  # TODO check python env is correct here
  
  if (inherits(nJobs, "numeric")) {
    nJobs <- as.integer(nJobs)
  }
  if (inherits(maxDepth, "numeric")) {
    maxDepth <- as.integer(maxDepth)
  }
  if (inherits(nTrees, "numeric")) {
    nTrees <- as.integer(nTrees)
  }
  
  checkIsClass(nTrees, c('integer', "character"))
  checkIsClass(maxDepth, c('integer'))
  if (inherits(nTrees, c("integer"))) {
    checkHigher(nTrees, 0)  
  } else {
    if (nTrees != "auto") {
      stop("nTrees should be either an integer or 'auto'")
    }
  } 
  
  checkHigher(maxDepth, 0)
  
  featureEngineeringSettings <- list(
    nTrees = nTrees,
    maxDepth = maxDepth
  )
  
  attr(featureEngineeringSettings, "fun") <- "borutaFeatureSelection"
  class(featureEngineeringSettings) <- "featureEngineeringSettings"
  
  return(featureEngineeringSettings)
}

borutaFeatureSelection <- function(
    trainData, 
    featureEngineeringSettings,
    covariateIdsInclude = NULL
){

  if(is.null(covariateIdsInclude)){
    #convert data into matrix:
    mappedData <- toSparseM(trainData)
    
    matrixData <- mappedData$dataMatrix
    labels <- mappedData$labels
    covariateMap <- mappedData$covariateMap
    
    X <- reticulate::r_to_py(matrixData)$toarray()
    y <- reticulate::r_to_py(matrix(labels$outcomeCount, ncol=1))
    
    sklearn <- reticulate::import('sklearn')
    BorutaPy <- reticulate::import('boruta')$BorutaPy

          

    
    rf = sklearn$ensemble$RandomForestClassifier(
      max_depth = featureEngineeringSettings$maxDepth,
      n_jobs = featureEngineeringSettings$nJobs, 
    )
    
    featureSelector <- BorutaPy(rf, n_estimators=featureEngineeringSettings$nTrees,
                                verbose=2L, random_state=42L)
    
    featureSelector$fit(X, y$squeeze())
    
    includedFeatures <- featureSelector$support_
    
    
    covariateIdsInclude <- covariateMap %>% 
      dplyr::filter(.data$columnId %in% which(includedFeatures)) %>%
      dplyr::select("covariateId") %>% dplyr::arrange("covariateId") %>%
      dplyr::pull()
  } 
  
  trainData$covariateData$covariates <- trainData$covariateData$covariates %>% 
    dplyr::filter(.data$covariateId %in% covariateIdsInclude)
  
  trainData$covariateData$covariateRef <- trainData$covariateData$covariateRef %>% 
    dplyr::filter(.data$covariateId %in% covariateIdsInclude)
  
  
  featureEngineering <- list(
    funct = 'borutaFeatureSelection',
    settings = list(
      featureEngineeringSettings = featureEngineeringSettings,
      covariateIdsInclude = covariateIdsInclude
    )
  )
  
  attr(trainData, 'metaData')$featureEngineering = listAppend(
    attr(trainData, 'metaData')$featureEngineering,
    featureEngineering
  )
  
  return(trainData)
  
}