
#' Create setting for GOSDT with python
#'
#' @param objective  Loss function (options include "acc", "bacc", "f1", "auc", "pauc")
#' @param regularization Note: We highly recommend setting the regularization to a value larger than 1/num_samples. A small regularization could lead to a longer training time.
#' @param maxDepth Used to set the maximum tree depth for solutions, counting a tree with just the root node as depth 1. 0 means unlimited.
#' @param warmLB 
#' @param pathToLabels 
#' @param timeLimit A time limit upon which the algorithm will terminate. If the time limit is reached, the algorithm will terminate with an error. When set to 0, no time limit is imposed.
#' @param similarSupport Enables the similar support bound imeplemented via the distance index. 
#' @param seed A seed when training the final model
#'
#' @return
#' @export
#'
#' @examples
setGOSDT <- function(
    objective = list("auc"),
    regularization = list(0.005), 
    maxDepth = list(0), 
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
