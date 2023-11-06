# @file SaveLoadData.R
#
# Copyright 2021 Observational Health Data Sciences and Informatics
#
# This file is part of CohortMethod
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


#' Save the cohort data to folder
#'
#' @description
#' \code{saveTrainTestData} saves an data object before runModelDevelopment to folder.
#'
#' @param data               Data object before runModelDevelopment.
#' @param file               The name of the folder where the data will be written. The folder should
#'                           not yet exist.
#' @param envir              The environment for to evaluate variables when saving.
#' @param overwrite          Whether to force overwrite an existing file
#' @details
#' The data will be written to a set of files in the folder specified by the user.
#'
#' @examples
#' # todo
#'
saveTrainTestData <- function(data, file, test=T, envir=NULL, overwrite=F) {
  if (missing(data)){
    stop("Must specify data")
  }
  if (missing(file)){
    stop("Must specify file")
  }
  if(dir.exists(file.path(file, "train-covariates"))){
    stop('Folder to save train covariates already exists...')
  }
  if(dir.exists(file.path(file, "test-covariates"))){
    stop('Folder to save test covariates already exists...')
  }
  if(!dir.exists(file)){
    dir.create(file, recursive = T)
  }

  if (test) { # Save train and test
    Andromeda::saveAndromeda(data$Train$covariateData, file = file.path(file, "train-covariates"), maintainConnection = T)
    saveRDS(data$Train$labels, file = file.path(file, "train-labels.rds"))
    saveRDS(data$Train$folds, file = file.path(file, "train-folds.rds"))
    
    Andromeda::saveAndromeda(data$Test$covariateData, file = file.path(file, "test-covariates"), maintainConnection = T)
    saveRDS(data$Test$labels, file = file.path(file, "test-labels.rds"))
    
    saveRDS(attr(data$Train, "metaData"), file = file.path(file, "metadata.rds"))
    saveRDS(attr(data$Train$covariateData, "metaData"), file = file.path(file, "covariatedata_metadata.rds"))
    
  } else {  # Save only train
    Andromeda::saveAndromeda(data$covariateData, file = file.path(file, "train-covariates"), maintainConnection = T)
    saveRDS(data$labels, file = file.path(file, "train-labels.rds"))
    saveRDS(data$folds, file = file.path(file, "train-folds.rds"))
    
    saveRDS(attr(data, "metaData"), file = file.path(file, "metadata.rds"))
    saveRDS(attr(data$covariateData, "metaData"), file = file.path(file, "covariatedata_metadata.rds"))
    
  }
}

#' Load the cohort data from a folder
#'
#' @description
#' \code{loadTrainTestData} loads an data object before runModelDevelopment from a folder in the file.
#' system.
#'
#' @param file       The name of the folder containing the data.
#' @param readOnly   If true, the data is opened read only.
#'
#' @details
#' The data will be written to a set of files in the folder specified by the user.
#'
#' @return
#' Data object before runModelDevelopment.
#'
#' @examples
#' # todo
#'
loadTrainTestData <- function(file, readOnly = TRUE) {
  if (!file.exists(file))
    stop(paste("Cannot find folder", file))
  if (!file.info(file)$isdir)
    stop(paste("Not a folder", file))
  
  result <- list(Train = list(covariateData = FeatureExtraction::loadCovariateData(file = file.path(file, "train-covariates")),
                              labels = readRDS(file.path(file, "train-labels.rds")),
                              folds = readRDS(file.path(file, "train-folds.rds"))),
                 Test = list(covariateData = FeatureExtraction::loadCovariateData(file = file.path(file, "test-covariates")),
                             labels = readRDS(file.path(file, "test-labels.rds"))))
  attr(result$Train, "metaData") <- readRDS(file.path(file, "metadata.rds"))
  attr(result$Train$covariateData, "metaData") <- readRDS(file.path(file, "covariatedata_metadata.rds"))
  
  class(result$Train) <- "plpData"
  class(result$Test) <- "plpData"
  
  return(result)
}