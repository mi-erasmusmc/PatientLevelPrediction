#' Create the settings for Boruta feature selection
#'
#' @details
#' From: https://doi.org/10.18637/jss.v036.i11
#'
#' @param nJobs       How many jobs to do in parallel
#' @param maxDepth    Max depth of each tree in the RandomForest used
#' @param nTrees      How many trees to use, default is `auto`
#' @param verbosity   0 for silent, 1 to display iteration number, 2 to display features selected as well
#' @param iterations  How many iterations to run `Boruta` for. Default: 100
#' @param randomState Either `NULL` or an integer. If integer it is the seed used by the random number generator
#'
#' @return
#' An object of class \code{featureEngineeringSettings}
#' @export
createBorutaFeatureSelection <- function(nJobs = -1L, 
                                         maxDepth = 5L,
                                         nTrees = "auto",
                                         verbosity = 2L,
                                         iterations = 100L,
                                         randomState = 42L
                                         ){
  # check python environment
  tryCatch(reticulate::import('numpy'), error = function(e) stop("Numpy must be available in python environment"))
  tryCatch(reticulate::import('boruta'), error= function(e) stop("Boruta must be installed in the python environment"))
  tryCatch(reticulate::import('sklearn'), error= function(e) stop("sklearn must be installed in the python environment"))

  # check input variables  
  if (inherits(nJobs, "numeric")) {
    nJobs <- as.integer(nJobs)
  }
  if (inherits(maxDepth, "numeric")) {
    maxDepth <- as.integer(maxDepth)
  }
  if (inherits(nTrees, "numeric")) {
    nTrees <- as.integer(nTrees)
  }
  if (inherits(verbosity, "numeric")) {
    verbosity <- as.integer(verbosity)
  }
  if (inherits(iterations, "numeric")) {
    iterations <- as.integer(iterations)
  }
  if (inherits(randomState, "numeric")) {
    randomState <- as.integer(randomState)
  }
  
  
  checkIsClass(nTrees, c('integer', "character"))
  checkIsClass(maxDepth, c('integer'))
  checkIsClass(randomState, c("integer", "NULL"))
  if (inherits(nTrees, c("integer"))) {
    checkHigher(nTrees, 0)  
  } else {
    if (nTrees != "auto") {
      stop("nTrees should be either an integer or 'auto'")
    }
  } 
  checkHigher(maxDepth, 0)
  checkHigher(iterations, 0)
  
  if (!verbosity %in% c(0L, 1L, 2L)) {
    stop(paste0("verbosity must be one of 0, 1, 2. You supplied: ", verbosity))
  }
  
  featureEngineeringSettings <- list(
    nTrees = nTrees,
    maxDepth = maxDepth,
    iterations = iterations,
    verbosity = verbosity,
    randomState = randomState,
    nJobs = nJobs
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
    
    X <- reticulate::r_to_py(matrixData)
    y <- reticulate::r_to_py(matrix(labels$outcomeCount, ncol=1))
    
    sklearn <- reticulate::import('sklearn')
    BorutaPy <- reticulate::import('boruta')$BorutaPy

    rf = sklearn$ensemble$RandomForestClassifier(
      max_depth = featureEngineeringSettings$maxDepth,
      n_jobs = featureEngineeringSettings$nJobs, 
    )
    
    featureSelector <- BorutaPy(rf, n_estimators=featureEngineeringSettings$nTrees,
                                verbose=featureEngineeringSettings$verbosity, 
                                random_state=featureEngineeringSettings$randomState,
                                max_iter=featureEngineeringSettings$iterations)
    
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