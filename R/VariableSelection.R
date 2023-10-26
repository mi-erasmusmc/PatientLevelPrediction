# @file VariableSelection.R
#
# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of PatientLevelPrediction
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


#' Create settings using normalized joint mutual information maximization (NJMIM) for feature selection
#' @param k number of variables to select, default `20`
#' @export
njmimSettings <- function(k = 20) {
  
  ensure_installed('praznik')
  
  checkIsClass(k, c('integer', 'numeric'))
  checkHigherEqual(k, 0)
  
  featureEngineeringSettings <- list(k = k) 
  
  attr(featureEngineeringSettings, "fun") <- "njmimFeatureSelection"
  class(featureEngineeringSettings) <- "featureEngineeringSettings"
  
  return(featureEngineeringSettings)
}

njmimFeatureSelection <- function(trainData,
                                  featureEngineeringSettings,
                                  covariateIdsInclude = NULL) {
  if (is.null(covariateIdsInclude)) {
    sparseData <- toSparseM(trainData, trainData$labels)
    denseMatrix <- as.matrix(sparseData$dataMatrix)
    dataFrame <- as.data.frame(denseMatrix)
    y <- sparseData$labels$outcomeCount
    
    selection <- praznik::NJMIM(X = dataFrame, 
                                Y = y,
                                k = featureEngineeringSettings$k)
    
    covariateIdsInclude <- sparseData$covariateMap %>%
      dplyr::filter(columnId %in% selection$selection) %>%
      dplyr::pull(covariateId)
  }
  
  trainData$covariateData$covariates <- trainData$covariateData$covariates %>%
    dplyr::filter(covariateId %in% covariateIdsInclude)
  trainData$covariateData$covariateRef <- trainData$covariateData$covariateRef %>%
    dplyr::filter(covariateId %in% covariateIdsInclude)
  
  featureEngineering <- list(
    funct = 'njmimFeatureSelection',
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

#' Create settings using univariate statistics for feature selection
#' @param corMethod which type of correlation to use, `pearson`, `kendall` or `spearman`. default `pearson`
#' @param k number of variables to select, default `20`
#'
#' @export
univariateSettings <- function(k = 20, # TODO: set dynamic number based on elbow
                               corMethod = "pearson"){
  
  if (inherits(k, 'numeric')) {
    k <- as.integer(k)
  }
  
  checkIsClass(k, 'integer')
  checkHigherEqual(k, 0)
  
  if (!corMethod %in% c("pearson", "kendall", "spearman")) {
    stop("corMethod needs to be either 'pearson', 'kendall' or 'spearman'")
  }
  
  featureEngineeringSettings <- list(corMethod = corMethod,
                                     k = k) 
  
  attr(featureEngineeringSettings, "fun") <- "univariateFeatureSelection"
  class(featureEngineeringSettings) <- "featureEngineeringSettings"
  
  return(featureEngineeringSettings)
}


univariateFeatureSelection <- function(trainData,
                                       featureEngineeringSettings,
                                       covariateIdsInclude = NULL) {
  if (is.null(covariateIdsInclude)) {
    sparseData <- toSparseM(trainData, trainData$labels)
    denseMatrix <- as.matrix(sparseData$dataMatrix)
    dataFrame <- as.data.frame(denseMatrix)
    y <- sparseData$labels$outcomeCount
    
    if (ncol(dataFrame) > featureEngineeringSettings$k) {
      # Select features based on univariate association with outcome
      correlation <- sapply(1:ncol(dataFrame), function(col) {
        stats::cor(dataFrame[,col], y, method = featureEngineeringSettings$corMethod)
      })
      names(correlation) <- 1:ncol(dataFrame)
      
      # Order from high to low absolute correlation
      correlation <-  correlation[order(abs(correlation), decreasing = TRUE)]
      
      # Select variables
      selected <- names(correlation)[1:min(featureEngineeringSettings$k, length(correlation))]
      
      covariateIdsInclude <- sparseData$covariateMap %>%
        dplyr::filter(columnId %in% selected) %>%
        dplyr::pull(covariateId)
      
    }
  }
  
  trainData$covariateData$covariates <- trainData$covariateData$covariates %>%
    dplyr::filter(covariateId %in% covariateIdsInclude)
  trainData$covariateData$covariateRef <- trainData$covariateData$covariateRef %>%
    dplyr::filter(covariateId %in% covariateIdsInclude)
  
  featureEngineering <- list(
    funct = 'univariateFeatureSelection',
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



#' Create settings using stepwise selection for feature selection
#' @param k number of variables to select, default `20`
#' @param selectMethod `backward` or `forward` selection
#' @param kStart number of variables to select initially
#' @param stepSize  How many variables to add/remove in each step
#' @param modelSettings settings of model to use in fit after selecting variables
#' @export
stepwiseSettings <- function(k = 20, # TODO: set dynamic number based on elbow
                             selectMethod = "backward",
                             kStart = 100, # TODO: set dynamic number based on elbow
                             stepSize = 1,
                             modelSettings = PatientLevelPrediction::setLassoLogisticRegression()) {
  
  checkIsClass(k, c('integer', 'numeric'))
  checkHigherEqual(k, 0)
  
  if (!selectMethod %in% c("forward", "backward")) {
    stop("selectMethod needs to be either 'forward' or 'backward")
  }
  
  checkIsClass(kStart, c('integer', 'numeric'))
  checkHigherEqual(kStart, 0)
  
  checkIsClass(stepSize, c('integer', 'numeric'))
  checkHigherEqual(stepSize, 0)
  
  featureEngineeringSettings <- list(k = k,
                                     selectMethod = selectMethod,
                                     kStart = kStart,
                                     stepSize = stepSize,
                                     modelSettings = modelSettings) 
  
  attr(featureEngineeringSettings, "fun") <- "stepwiseFeatureSelection"
  class(featureEngineeringSettings) <- "featureEngineeringSettings"
  
  return(featureEngineeringSettings)
}


stepwiseFeatureSelection <- function(trainData,
                                     featureEngineeringSettings,
                                     covariateIdsInclude = NULL) {
  if (is.null(covariateIdsInclude)) {
    search <- 999 # TODO: which number to use?
    analysisId <- 999 # TODO: which number to use?
    
    # Initial fit PLP model
    fun <- eval(parse(text = featureEngineeringSettings$modelSettings$fitFunction))
    args <- list(
      trainData = trainData,
      modelSettings = featureEngineeringSettings$modelSettings,
      search = search,
      analysisId = analysisId
    )
    plpModel <- do.call(fun, args)
    
    if (featureEngineeringSettings$selectMethod == "backward") {
      # Initial selection
      # Select non-zero coefficients for LASSO
      covIds <- plpModel$model$coefficients$covariateIds[plpModel$model$coefficients$betas != 0] 
      
      # TODO: make generic across algorithms based on var importance (e.g. randomForest)
      # covIds <- TODO
      
      fullCovariates <- NULL
      updateIteration <- "remove"
      
    } else if (featureEngineeringSettings$selectMethod == "forward") {
      # Initial selection (empty)
      covIds <- NULL
      
      fullCovariates <- as.data.frame(trainData$covariateData$covariates)
      updateIteration <- "add"
      
    } else {
      stop("Variable selection stopped: selectMethod not implemented.")
    }
    
    update <- "select" # Only first update of covariates
    update_trainData <- trainData
    
    while(!is.null(covIds) | update == "select") { # Stop when selected is NULL
      # Update covariate data
      update_trainData$covariateData <- updateCovariateData(update_trainData$covariateData, covIds=covIds, update=update, fullCovariates=fullCovariates)
      
      # Backward or forward select variables
      update <- updateIteration 
      covIds <- selectVariables(featureEngineeringSettings, update_trainData, fullCovariates, search, analysisId)
    }
    
    covariates <- as.data.frame(update_trainData$covariateData$covariates)
    covariateIdsInclude <- unique(covariates$covariateId)
  }
  
  trainData$covariateData$covariates <- trainData$covariateData$covariates %>%
    dplyr::filter(covariateId %in% covariateIdsInclude)
  trainData$covariateData$covariateRef <- trainData$covariateData$covariateRef %>%
    dplyr::filter(covariateId %in% covariateIdsInclude)
  
  featureEngineering <- list(
    funct = 'univariateFeatureSelection',
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

# Help function stepwiseFeatureSelection
selectVariables <- function(featureEngineeringSettings, trainData, fullCovariates, search, analysisId) {
  
  # Settings variable selection
  covariates <- as.data.frame(trainData$covariateData$covariates)
  
  if (featureEngineeringSettings$selectMethod == "backward") {
    update <- "remove"
    start <- (length(unique(covariates$covariateId)) > featureEngineeringSettings$k) # Too many variables
    covariateList <- unique(covariates$covariateId)
  } else if (featureEngineeringSettings$selectMethod == "forward") {
    update <- "add"
    start <- (length(unique(covariates$covariateId)) < featureEngineeringSettings$k) # Not enough variables
    covariateList <- unique(fullCovariates$covariateId)[!(unique(fullCovariates$covariateId) %in% unique(covariates$covariateId))]
  }
  # TODO: make correction for featureEngineeringSettings$stepSize to come to exactly the right number of variables
  
  if (start) {
    performance <- sapply(covariateList, function(covId) {
      # Temporarily update covarite data
      tempData <- trainData
      tempData$covariateData <- updateCovariateData(tempData$covariateData, covIds=covId, update=update, fullCovariates=fullCovariates)
      
      # Re-fit PLP model
      fun <- eval(parse(text = featureEngineeringSettings$modelSettings$fitFunction))
      args <- list(
        trainData = tempData,
        modelSettings = featureEngineeringSettings$modelSettings,
        search = search,
        analysisId = analysisId
      )
      tempModel <- do.call(fun, args)
      
      # Return performance for current covariate data
      # TODO: extend to other selection criteria? 
      metric <- tempModel$model$log_likelihood
      
      return(metric)
    })
    names(performance) <- covariateList
    
    if (featureEngineeringSettings$selectMethod == "backward") {
      # Order from low to high performance (minimum negative log likelihood is high)
      performance <- performance[order(performance, decreasing = TRUE)]
      
      # Select variables to remove
      covIds <- names(performance)[1:min(featureEngineeringSettings$stepSize, length(performance))]
      
    } else if (featureEngineeringSettings$selectMethod == "forward") {
      # Order from high to low performance (minimum negative log likelihood is high)
      performance <- performance[order(performance, decreasing = FALSE)]
      
      # Select variables to add
      covIds <- names(performance)[1:min(featureEngineeringSettings$stepSize, length(performance))]
    }
    
    return(covIds)
  }
  
  return(NULL) # Return NULL to initiate stop
}

# Help function stepwiseFeatureSelection
updateCovariateData <- function(covariateData, covIds, update="select", fullCovariates=NULL) {
  # FeatureExtraction -> excludedCovariateConceptIds: A list of concept IDs that should NOT be used to construct covariates.
  newCovariates <- as.data.frame(covariateData$covariates) # TODO: try without
  print(paste0("Rows covariates - before: ", nrow(newCovariates)))
  
  if (update == "select") {
    newCovariates <- newCovariates[newCovariates$covariateId %in% covIds,]
  } else if (update == "remove") {
    # print(paste0("Remove covIds: ", paste0(covIds, collapse = ", ")))
    newCovariates <- newCovariates[!(newCovariates$covariateId %in% covIds),]
  } else if (update == "add") {
    # print(paste0("Add covIds: ", paste0(covIds, collapse = ", ")))
    
    if (is.null(fullCovariates)) {
      stop("fullCovariates missing, needed for adding covariates.")
    }
    newCovariates <- fullCovariates[fullCovariates$covariateId %in% c(unique(newCovariates$covariateId), covIds),]
  }
  print(paste0("Rows covariates - after: ", nrow(newCovariates)))
  
  # metaData <- list(call = match.call())
  result <- Andromeda::andromeda(covariates = newCovariates,
                                 labels=covariateData$labels,
                                 covariateRef = covariateData$covariateRef,
                                 analysisRef = covariateData$analysisRef)
  # attr(result, "metaData") <- metaData
  class(result) <- "CovariateData"
  
  return(result)
}


