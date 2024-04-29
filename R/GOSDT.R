
#' Create setting for GOSDT with python
#'
#' @param ntrees    (list) The number of trees to build 
#' @param maxDepth   (list) The maximum depth of the tree. If NULL, then nodes are expanded until all leaves are pure or until all leaves contain less than minSamplesSplit samples.
#' @param seed  A seed when training the final model
#'
#' @examples
#' \dontrun{
#' model.gosdt <- setGOSDT(...)
#' }       
#' Adapted from SklearnClassifierSettings.R                    
#' @export
setGOSDT <- function(
    objective = list("auc"),
    regularization = list(0.005),
    maxDepth = list(0), # depth_budget
    warmLB = list(TRUE),
    pathToLabels = list("warm_label.tmp"),
    timeLimit = list(3600),
    similarSupport = list(FALSE),
    seed = sample(100000,1)
){
  
  checkIsClass(seed, c('numeric','integer'))
  checkIsClass(objective, c('list'))
  checkIsClass(regularization, c('list'))
  checkIsClass(maxDepth, c('list'))
  checkIsClass(warmLB, c('list'))
  checkIsClass(pathToLabels, c('list'))
  checkIsClass(timeLimit, c('list'))
  checkIsClass(similarSupport, c('list'))
  
  # convert to integer when needed
  for(i in 1:length(maxDepth)){
    if(inherits(x = maxDepth[[i]], what =  c("numeric", "integer"))){
      maxDepth[[i]] <- as.integer(maxDepth[[i]])
    }
  }
  
  # add value checks
  paramGrid = list(
    objective = objective,
    regularization =  regularization,
    maxDepth = maxDepth,
    warmLB = warmLB,
    pathToLabels = pathToLabels,
    timeLimit = timeLimit,
    similarSupport = similarSupport,
    seed = list(as.integer(seed[[1]]))
  )
  param <- listCartesian(paramGrid)
  
  attr(param, 'settings') <- list(
    modelType = 'GOSDT',
    seed = seed[[1]],
    paramNames = names(paramGrid), #use this for logging params
    requiresDenseMatrix = T,
    name = "General Optimal Sparse Decision Tree",
    pythonModule = 'gosdt.model.gosdt',
    pythonClass = 'GOSDT'
  ) 
  
  attr(param, 'saveToJson') <- F
  attr(param, 'saveType') <- 'file'
  
  result <- list(
    fitFunction = "fitSklearn",
    param = param
  )
  class(result) <- "modelSettings"
  
  return(result)
}

# Adapted from SklearnClassifierSettings.R
GOSDTInputs <- function(classifier, param){
  
  model <- classifier(
    configuration=list(
      objective = param[[which.max(names(param)=='objective')]],
      regularization = param[[which.max(names(param)=='regularization')]],
      depth_budget = param[[which.max(names(param)=='maxDepth')]],
      warm_LB = param[[which.max(names(param)=='warmLB')]],
      path_to_labels = param[[which.max(names(param)=='pathToLabels')]],
      time_limit = param[[which.max(names(param)=='timeLimit')]],
      similar_support = param[[which.max(names(param)=='similarSupport')]],
      seed = param[[which.max(names(param)=='seed')]],
      verbose = F
    )
  )
  
  return(model)
}
