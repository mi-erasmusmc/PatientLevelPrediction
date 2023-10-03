# Copyright 2023 Observational Health Data Sciences and Informatics
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
#' @param k number of features to select 
#' @export
njmimSettings <- function(k=20) {
  
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