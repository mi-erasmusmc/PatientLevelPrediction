# @file EXPLORE.R
#
# Copyright 2020 Observational Health Data Sciences and Informatics
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

#' Create setting for EXPLORE (https://github.com/mi-erasmusmc/explore/)
#'
#' @param variableNumer       An option to add a seed when training the final model
#'
#' @examples
#' model.explore <- setEXPLORE(variableNumber=2000)
#'
#' @export
setEXPLORE <- function(output_path, variableNumber=2000){
  
  # check input
  if(length(variableNumber)!=1)
    stop('Can only currently enter a single value for variableNumber')
  
  if(!class(variableNumber) %in% c("numeric", "integer"))
    stop('Can incorrect class for variableNumber - must be numeric')
  
  result <- list(model='fitEXPLORE', param= list('output_path'=output_path, 'variableNumber'=variableNumber), name='EXPLORE')
  class(result) <- 'modelSettings' 
  
  return(result)
}


# EXPLORE
fitEXPLORE <- function(population, plpData, param, quiet=F,
                      outcomeId, cohortId, ...){
  
  # check plpData is libsvm format or convert if needed
  if (!FeatureExtraction::isCovariateData(plpData$covariateData)){
    stop("Needs correct covariateData")
  }
  
  if(colnames(population)[ncol(population)]!='indexes'){
    warning('indexes column not present as last column - setting all index to 1')
    population$indexes <- rep(1, nrow(population))
  }
  
  # check logger
  if(length(ParallelLogger::getLoggers())==0){
    logger <- ParallelLogger::createLogger(name = "SIMPLE",
                                           threshold = "INFO",
                                           appenders = list(ParallelLogger::createConsoleAppender(layout = ParallelLogger::layoutTimestamp)))
    ParallelLogger::registerLogger(logger)
  }
  
  if(!quiet)
    ParallelLogger::logTrace('Training EXPLORE model')
  
  start <- Sys.time()
  
  # make sure population is ordered?
  prediction <- population
  population$rowIdPython <- population$rowId - 1  # -1 to account for python/r index difference
  pPopulation <- as.matrix(population[,c('rowIdPython','outcomeCount','indexes')])
  
  covariateRef <- as.data.frame(plpData$covariateData$covariateRef)
  
  # convert plpData in coo to python:
  x <- toSparseM(plpData, population, map = NULL)
  
  # save the model to outLoc TODO: make this an input or temp location?
  outLoc <- createTempModelLoc()
  # clear the existing model pickles
  for(file in dir(outLoc))
    file.remove(file.path(outLoc,file))
  
  # then run standard python code
  e <- environment()
  reticulate::source_python(system.file(package='PatientLevelPrediction','python','featureSelection.py'), envir = e)
  pdata <- reticulate::r_to_py(x$data)
  
  # initial variable selection
  selection <- univariate_feature_selection(population=pPopulation, 
                                            plpData=pdata, 
                                            variableNumber = as.integer(param$variableNumber),
                                            quiet = quiet)
  
  pred <- as.data.frame(selection[[3]])
  colnames(pred) <- c('rowId','outcomeCount','indexes') # todo: add last column value
  attr(pred, "metaData") <- list(predictionType="binary")
  
  # add 1 to rowId from python:
  pred$rowId <- pred$rowId+1
  pred$value <- NA
  
  var_sel <- as.data.frame(cbind(selection[[2]],selection[[1]]))
  # var_sel <- as.data.frame(sapply(var_sel, function(col) factor(col, levels = c(0,1))), stringsAsFactors = TRUE)
  colnames(var_sel)[1] <- "y"
  
  # train model
  library('Explore')
  modelTrained <- Explore::trainExplore(output_path = param$output_path, train_data = var_sel, ClassFeature = "'y'", PositiveClass = 1)
  
  # cross validation
  for(i in 1:max(population$indexes)) {
    hold_out <- which(population[population$indexes > 0,]$indexes==i)
    
    subset_fit <- Explore::trainExplore(output_path = param$output_path, train_data = var_sel[-hold_out,], ClassFeature = "'y'", PositiveClass = 1)
    
    # print model
    print(subset_fit) 
    
    subset_predict <- as.numeric(Explore::predictExplore(model = subset_fit, test_data = var_sel[hold_out,]))
    pred$value[hold_out] <- subset_predict
    
    auc <- aucWithoutCi(subset_predict, pred$outcomeCount[pred$rowId %in% hold_out])
    writeLines(paste0('Model obtained CV AUC of ', auc, ' in fold ', i))
  }
  
  auc <- computeAuc(pred)
  writeLines(paste0('Model obtained CV AUC of ', auc))
  
  # get prediction on test set:
  ParallelLogger::logInfo('Getting predictions on train set')
  prediction <- merge(prediction, pred[,c('rowId', 'value')], by='rowId')
  
  # get the univeriate selected features (explore requires dense so need feat sel)
  varImp <- selection[[4]]
  varImp[is.na(varImp)] <- 0
  if(mean(varImp)==0)
    stop('No important variables - seems to be an issue with the data')
  
  topN <- varImp[order(-varImp)][param$variableNumber]
  inc <- which(varImp>=topN, arr.ind=T)
  
  incs <- rep(0, nrow(covariateRef))
  incs[inc] <- 1
  covariateRef$included <- incs
  covariateRef$covariateValue <- varImp
  
  comp <- start-Sys.time()
  
  result <- list(model = modelTrained,
                 modelSettings = list(model='EXPLORE', modelParameters=param),
                 trainCVAuc = auc,
                 hyperParamSearch = NULL,
                 metaData = plpData$metaData,
                 populationSettings = attr(population, 'metaData'),
                 outcomeId=outcomeId,
                 cohortId=cohortId,
                 varImp = covariateRef,
                 trainingTime=comp,
                 dense=1,
                 covariateMap=x$map,
                 predictionTrain = prediction
  )
  
  class(result) <- 'plpModel'
  attr(result, 'type') <- 'explore'
  attr(result, 'predictionType') <- 'binary'
  return(result)
}
