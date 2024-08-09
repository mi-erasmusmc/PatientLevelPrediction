# @file RIPPER.R
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

#' Create setting for RIPPER
#' @param variableSelection 
#' @param variableNumber 
#'
#' @export
setRIPPER <- function( # TODO: check default settings
  variableSelection = PatientLevelPrediction::setUnivariateSelection(),
  variableNumber = 10,
  saveDirectory = getwd()){
  
  # TODO: check input
  
  param <- list(variableSelection = variableSelection,
                variableNumber = variableNumber,
                saveDirectory = saveDirectory)
  
  attr(param, 'settings') <- list(
    modelType = 'ripper',
    modelName = 'RIPPER'
  )
  
  attr(param, 'modelType') <- 'binary'
  attr(param, 'saveType') <- 'RtoJson'
  
  result <- list(
    fitFunction = "fitRIPPER",
    param = param
  )
  class(result) <- 'modelSettings' 
  
  return(result)
}

#' @export
fitRIPPER <- function(trainData,
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
  
  denseData <- convertToDenseData(trainData, param$variableSelection, search, analysisId, param$saveDirectory)
  
  # convert age to groups
  if ("1002" %in% colnames(denseData)) {
    denseData['1002'] <- cut(denseData[['1002']],breaks=c(0,25,50,75,100),labels=c('0-25','25-50','50-75','75-100'))
  }
  
  # convert to factors (JRip cannot handle numeric features)
  # binary_cols <- sapply(1:ncol(denseData), function(c) all(denseData[[c]] %in% 0:1))
  # denseData[binary_cols] <- as.data.frame(sapply(denseData[binary_cols], function(col) factor(col, levels = c(0,1))), stringsAsFactors = TRUE)
  denseData <- as.data.frame(sapply(denseData, function(col) factor(col, levels = unique(col))), stringsAsFactors = TRUE)
  
  # train model
  fit <- tryCatch({
    ParallelLogger::logInfo('Running RIPPER')
    JRip(outcomeCount ~ . , data = denseData)
  },
  finally = ParallelLogger::logInfo('Done.')
  )
  
  ParallelLogger::logTrace('Returned from fitting RIPPER')
  comp <- Sys.time() - start
  
  ParallelLogger::logInfo(paste0("RIPPER rule: ", fit$classifier$toString()))
  
  ParallelLogger::logTrace('Getting variable importance')
  # Get the features selected using RIPPER
  featureNames <- colnames(denseData)[-1]
  varImp <- data.frame(
    covariateId = as.double(featureNames),
    value = sapply(featureNames, function(f) ifelse(grepl(f, fit$classifier$toString()), 1, 0))
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
                        coefficients = varImp$covariateId[varImp$value==1])
  
  # Getting predictions on train set:
  tempModel <- list(model = modelTrained)
  attr(tempModel, "modelType") <- attr(param, 'modelType')
  prediction <- predictRIPPER(
    plpModel = tempModel,
    cohort = trainData$labels, 
    data = trainData
  )
  prediction$evaluationType <- 'Train'
  
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
  attr(result, 'predictionFunction') <- 'predictRIPPER'
  attr(result, 'modelType') <- attr(param, 'modelType')
  attr(result, 'saveType') <- attr(param, 'saveType')
  return(result)
}

convertToDenseData <- function(trainData, modelSettings, search, analysisId, saveDirectory) {
  
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
  
  denseData <- merge(trainData$labels[c("rowId", "outcomeCount")], denseData, by = 'rowId', all.x = TRUE)
  denseData[is.na(denseData)] <- 0
  denseData$rowId <- NULL
  
  return(denseData)
}

predictRIPPER <- function(plpModel, data, cohort) {
  
  varSelection <- names(attr(plpModel$model$fit$terms, "dataClasses"))
  
  # Convert to dense covariates
  covariates <- as.data.frame(data$covariateData$covariates)
  covariates <- covariates[covariates$covariateId %in% varSelection,] # Select only covariates included in model
  denseData <- reshape2::dcast(covariates, rowId ~ covariateId, value.var = 'covariateValue', fill = 0)
  
  # convert age to groups
  if ("1002" %in% colnames(denseData)) {
    denseData['1002'] <- cut(denseData[['1002']],breaks=c(0,25,50,75,100),labels=c('0-25','25-50','50-75','75-100'))
  }
  
  # convert to factors (JRip cannot handle numeric features)
  # binary_cols <- sapply(1:ncol(denseData), function(c) all(denseData[[c]] %in% 0:1))
  # denseData[binary_cols] <- as.data.frame(sapply(denseData[binary_cols], function(col) factor(col, levels = c(0,1))), stringsAsFactors = TRUE)
  denseData <- as.data.frame(sapply(denseData, function(col) factor(col, levels = unique(col))), stringsAsFactors = TRUE)
  
  # Check if all covariates in data (in case no observations in test set with record)
  addCols <- varSelection[!(varSelection %in% c(colnames(denseData), "outcomeCount"))]
  denseData[addCols] <- 0
  
  denseData <- merge(cohort[c("rowId", "outcomeCount")], denseData, by = 'rowId', all.x = TRUE)
  denseData[is.na(denseData)] <- 0
  denseData[c("rowId", "outcomeCount")] <- NULL
  
  prediction <- data.frame(rowId=cohort$rowId, value=as.numeric(stats::predict(plpModel$model$fit, denseData)==1))
  
  # return the cohorts as a data frame with the prediction added as 
  # a new column with the column name 'value'
  prediction <- merge(cohort, prediction, by='rowId', all.x=T)
  attr(prediction, "metaData")$modelType <-  plpModel$model$modelType
  
  return(prediction)
}
