# @file EXPLORE.R
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
# Unless required by applicable v or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

#' Create setting for EXPLORE
#' @param variableSelection 
#' @param variableNumber 
#' @param startRulelength 
#' @param endRulelength 
#' @param operatorMethod 
#' @param cutoffMethod 
#' @param featureInclude 
#' @param maximize 
#' @param accuracy 
#' @param specificity 
#' @param printSettings 
#' @param printPerformance 
#' @param subsumption 
#' @param branchBound 
#' @param parallel 
#'
#' @export
setExplore <- function( # TODO: check default settings
  variableSelection = PatientLevelPrediction::setUnivariateSelection(),
  variableNumber = 10, # TODO: can be removed???
  startRulelength = 1,
  endRulelength = 3,
  operatorMethod = "EXHAUSTIVE",
  cutoffMethod = "RVAC",
  featureInclude = "",
  maximize = "BALANCEDACCURACY",
  accuracy = 0,
  specificity = 0,
  printSettings = TRUE,
  printPerformance = TRUE,
  subsumption = TRUE,
  branchBound = TRUE,
  parallel = FALSE,
  modelsCurve = FALSE,
  sort_by = "none",
  saveDirectory = getwd()){
  
  # TODO: check input
  
  param <- list(variableSelection = variableSelection,
                startRulelength = startRulelength,
                endRulelength = endRulelength,
                operatorMethod = operatorMethod,
                cutoffMethod = cutoffMethod,
                featureInclude = featureInclude,
                maximize = maximize,
                accuracy = accuracy,
                specificity = specificity,
                printSettings = printSettings,
                printPerformance = printPerformance,
                subsumption = subsumption,
                branchBound = branchBound,
                parallel = parallel,
                modelsCurve = modelsCurve,
                sort_by = sort_by,
                saveDirectory = saveDirectory)
  
  attr(param, 'settings') <- list(
    modelType = 'explore',
    modelName = 'EXPLORE'
  )
  
  attr(param, 'modelType') <- 'binary'
  attr(param, 'saveType') <- 'RtoJson'
  
  result <- list(
    fitFunction = "fitExplore",
    param = param
  )
  class(result) <- 'modelSettings' 
  
  return(result)
}

#' @export
fitExplore <- function(trainData,
                       modelSettings,
                       search = 'none',
                       analysisId,
                       ...) {
  
  param <- modelSettings$param
  
  # check plpData is coo format:
  if (!FeatureExtraction::isCovariateData(trainData$covariateData)){
    stop("Needs correct covariateData")
  }
  
  settings <- attr(param, 'settings')
  
  start <- Sys.time()
  
  exploreData <- convertToExploreData(trainData, param$variableSelection, search, analysisId, param$saveDirectory)
  
  exploreData <- sortData(exploreData, sort_by = param$sort_by)
  
  # train model
  fit <- tryCatch({
    ParallelLogger::logInfo('Running Explore')
    Explore::trainExplore(output_path = file.path(param$saveDirectory, "Explore"), train_data = exploreData,
                          ClassFeature = "'outcomeCount'", PositiveClass = "1",
                          StartRulelength = param$startRulelength, EndRulelength = param$endRulelength, 
                          OperatorMethod = param$operatorMethod, CutoffMethod = param$cutoffMethod,
                          FeatureInclude = param$featureInclude, Maximize = param$maximize,
                          Accuracy = param$accuracy, Specificity = param$specificity,
                          Subsumption = param$subsumption, BranchBound = param$branchBound,
                          Parallel = param$parallel)
  },
  finally = ParallelLogger::logInfo('Done.')
  )
  
  ParallelLogger::logTrace('Returned from fitting EXPLORE')
  comp <- Sys.time() - start
  
  ParallelLogger::logTrace('Getting variable importance')
  # Get the features selected using EXPLORE
  # TODO: add if fit NULL?
  vars <- unlist(stringr::str_match_all(fit, "'\\d*'"))
  vars <- stringr::str_remove_all(vars, "'")
  featureNames <- colnames(exploreData)[-1]
  varImp <- data.frame(
    covariateId = as.double(featureNames),
    value = sapply(featureNames, function(f) ifelse(f %in% vars, 1, 0))
  )
  
  variableImportance <- data.frame(
    covariateId = varImp$covariateId,
    covariateValue = varImp$value
  )
  
  # variableImportance <- trainData$covariateData$covariateRef %>% 
  #   dplyr::collect()  %>%
  #   dplyr::left_join(varImp, by = 'covariateId') %>%
  #   # dplyr::mutate(covariateValue =.data$value) %>%
  #   dplyr::mutate(covariateValue = ifelse(is.na(.data$value), 0, .data$value)) %>%
  #   dplyr::select(-.data$value) %>%
  #   dplyr::arrange(-abs(.data$covariateValue)) %>%
  #   dplyr::collect()
  
  modelTrained <-  list(fit = fit,
                        coefficients = vars)
  
  # Getting predictions on train set:
  tempModel <- list(model = modelTrained)
  attr(tempModel, "modelType") <- attr(param, 'modelType')
  prediction <- predictExplore(
    plpModel = tempModel,
    cohort = trainData$labels, 
    data = trainData
  )
  prediction$evaluationType <- 'Train'
  
  # Generate models for AUC curve:
  if (param$modelsCurve) {
    ParallelLogger::logInfo('Running Explore for different sensitivities/specificities')
    modelTrained[["modelsCurve"]] <- Explore::modelsCurveExplore(output_path = file.path(param$saveDirectory, "Explore"), train_data = exploreData,
                                                                 ClassFeature = "'outcomeCount'", PositiveClass = "1",
                                                                 StartRulelength = param$startRulelength, EndRulelength = param$endRulelength, 
                                                                 OperatorMethod = param$operatorMethod, CutoffMethod = param$cutoffMethod,
                                                                 FeatureInclude = param$featureInclude, Maximize = param$maximize,
                                                                 Accuracy = param$accuracy, Specificity = param$specificity,
                                                                 Subsumption = param$subsumption, BranchBound = param$branchBound,
                                                                 Parallel = param$parallel)
    # saveRDS(models, file = file.path(param$saveDirectory, "Explore", "modelsCurve"))
  }
  
  # TODO: save these models for later
  # How to evaluate this for test case? Use predict?
  
  # cross-validation
  # for(i in 1:max(population$indexes)) {
  #   hold_out <- which(population[population$indexes > 0,]$indexes==i)
  # 
  #   subset_fit <- Explore::trainExplore(output_path = param$output_path, train_data = var_sel[-hold_out,], ClassFeature = "'y'", PositiveClass = 1,
  #                                       StartRulelength = param$start_rule_length, EndRulelength = param$end_rule_length,
  #                                       FeatureInclude = param$feature_include, Specificity = param$specificity)
  # 
  #   # print model
  #   print(subset_fit)
  # 
  #   subset_predict <- as.numeric(Explore::predictExplore(model = subset_fit, test_data = var_sel[hold_out,]))
  #   pred$value[hold_out] <- subset_predict
  # 
  #   auc <- aucWithoutCi(subset_predict, pred$outcomeCount[pred$rowId %in% hold_out])
  #   writeLines(paste0('Model obtained CV AUC of ', auc, ' in fold ', i))
  # }
  # 
  # auc <- computeAuc(pred)
  # writeLines(paste0('Model obtained CV AUC of ', auc))
  
  result <- list(
    model = modelTrained,
    
    preprocessing = list(
      featureEngineering = attr(trainData, "metaData")$featureEngineering,#learned mapping
      tidyCovariates = attr(trainData$covariateData, "metaData")$tidyCovariateDataSettings,  #learned mapping
      requireDenseMatrix = F
    ),
    
    prediction = prediction,
    
    modelDesign = PatientLevelPrediction::createModelDesign(
      targetId = attr(trainData, "metaData")$targetId, # added
      outcomeId = attr(trainData, "metaData")$outcomeId, # added
      restrictPlpDataSettings = attr(trainData, "metaData")$restrictPlpDataSettings, # made this restrictPlpDataSettings
      covariateSettings = attr(trainData, "metaData")$covariateSettings,
      populationSettings = attr(trainData, "metaData")$populationSettings, 
      featureEngineeringSettings = attr(trainData, "metaData")$featureEngineeringSettings,
      preprocessSettings = attr(trainData$covariateData, "metaData")$preprocessSettings,
      modelSettings = modelSettings, # modified
      splitSettings = attr(trainData, "metaData")$splitSettings,
      sampleSettings = attr(trainData, "metaData")$sampleSettings
    ),
    
    trainDetails = list(
      analysisId = analysisId, 
      analysisSource = '', #TODO: add from model
      developmentDatabase = attr(trainData, "metaData")$cdmDatabaseSchema,
      attrition = attr(trainData, "metaData")$attrition, 
      trainingTime =  paste(as.character(abs(comp)), attr(comp,'units')),
      trainingDate = Sys.Date(),
      modelName = settings$modelType,
      finalModelParameters = list() #TODO: add parameters
      # hyperParamSearch = cvPerFold
    ),
    
    covariateImportance = variableImportance
  )
  
  class(result) <- 'plpModel'
  attr(result, 'predictionFunction') <- 'predictExplore'
  attr(result, 'modelType') <- attr(param, 'modelType')
  attr(result, 'saveType') <- attr(param, 'saveType')
  return(result)
}

# TODO: change to convertToDenseData?
convertToExploreData <- function(trainData, modelSettings, search, analysisId, saveDirectory) {
  
  # Apply pre-variable selection
  if (!is.null(modelSettings)) {
    # Adjust settings to return data
    modelSettings$param$param$returnData <- TRUE
    modelSettings$param$param$saveDirectory <- stringr::str_remove(saveDirectory, analysisId)
    
    # TODO: always if # covariate > ?
    fun <- eval(parse(text = modelSettings$fitFunction))
    args <- list(
      trainData = trainData,
      modelSettings,
      search = search,
      analysisId = analysisId
    )
    trainData <- do.call(fun, args)
  }
  
  # Convert to dense covariates
  covariates <- as.data.frame(trainData$covariateData$covariates)
  denseData <- reshape2::dcast(covariates, rowId ~ covariateId, value.var = 'covariateValue', fill = 0)
  
  exploreData <- merge(trainData$labels[c("rowId", "outcomeCount")], denseData, by = 'rowId', all.x = TRUE)
  exploreData[is.na(exploreData)] <- 0
  exploreData$rowId <- NULL
  
  return(exploreData)
}


sortData <- function(exploreData, sort_by = "none", corMethod = "pearson") {
  
  if (sort_by == "none") {
    return(exploreData)
  } else if (sort_by == "random") {
    
    # Random shuffle of covariates (except outcomeCount)
    shuffle <- sample(2:ncol(exploreData), replace = FALSE)
    
    # Order randomly
    exploreData <-  exploreData[c(1, shuffle)]
    
    return(exploreData)
  } else if (sort_by == "correlation") {
    
    # Compute univariate correlation
    correlation <- sapply(2:ncol(exploreData), function(covId) {
      cor(exploreData[covId], exploreData$outcomeCount, method = corMethod)
    })
    
    # Order from high to low absolute correlation
    exploreData <-  exploreData[c(1, order(abs(correlation), decreasing = TRUE) + 1)]
    
    return(exploreData)
  }
  
}

predictExplore <- function(plpModel, data, cohort) {
  
  # Convert to dense covariates
  covariates <- as.data.frame(data$covariateData$covariates)
  covariates <- covariates[covariates$covariateId %in% plpModel$model$coefficients,] # Select only covariates included in model
  denseData <- reshape2::dcast(covariates, rowId ~ covariateId, value.var = 'covariateValue', fill = 0)
  
  exploreData <- merge(cohort[c("rowId", "outcomeCount")], denseData, by = 'rowId', all.x = TRUE)
  exploreData[is.na(exploreData)] <- 0
  exploreData[c("rowId", "outcomeCount")] <- NULL
  
  prediction <- data.frame(rowId=cohort$rowId, value=as.numeric(Explore::predictExplore(model = plpModel$model$fit, test_data = exploreData)))
  
  # return the cohorts as a data frame with the prediction added as 
  # a new column with the column name 'value'
  prediction <- merge(cohort, prediction, by='rowId', all.x=T)
  attr(prediction, "metaData")$modelType <-  plpModel$model$modelType
  
  return(prediction)
}
