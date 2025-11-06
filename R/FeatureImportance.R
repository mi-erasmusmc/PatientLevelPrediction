# @file FeatureEngineering.R
# Copyright 2025 Observational Health Data Sciences and Informatics
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

# reticulate::virtualenv_install(packages="shap")

#' Model Agnostic Feature Importance
#' @export
computeFI <- function(plpResult, population, plpData, method="PFI", params = NULL,
                      convergence = TRUE, covariates = NULL, cores = NULL, log = NULL,
                      logthreshold = "INFO") {
  if (method == "PFI_R") {
    fi <- pfi(plpResultEunomia, populationEunomia, plpDataEunomia)
  } else if (method %in% c("PFI", "SHAP")) {
    fi <- fiSklearn(plpResultEunomia, populationEunomia, plpDataEunomia, method, params, convergence)
  }
  
  return(fi)
}

#' Model Agnostic Feature Importance
#' @export
fiSklearn <- function(plpResult, population, plpData, method="PFI", params = NULL,
                      convergence = TRUE, covariates = NULL, cores = NULL, log = NULL,
                      logthreshold = "INFO") {
  if (!is.null(log)) {
    appender <- ParallelLogger::createFileAppender(
      layout = ParallelLogger::layoutParallel,
      fileName = log
    )
    
    logger <- ParallelLogger::createLogger(
      name = "PAR",
      threshold = logthreshold,
      appenders = list(appender)
    )
    ParallelLogger::registerLogger(logger)
  }
  
  if (is.null(covariates)) { # TODO: why is this restricted to positive covariateImportance? 
    covariates <- plpResult$model$covariateImportance %>%
      # dplyr::filter(.data$covariateValue != 0) %>%
      dplyr::select("covariateId") %>%
      dplyr::pull()
  }
  
  # TODO: check if input is sklearn model (or "pythonModule" exists) ?
  # TODO: check correct input params per FI method
  
  params <- listCartesian(params)
  
  # if convergence = FALSE use first params
  if (!convergence) {
    params <- list(params[[1]])
  }
  # TODO: how to sort for convergence?
  
  # convert
  matrixObjects <- toSparseM(
    plpData = plpData,
    cohort = population,
    map = plpResult$model$covariateImportance %>%
      dplyr::select("columnId", "covariateId")
  )
  matrixData <- matrixObjects$dataMatrix
  labels <- matrixObjects$labels
  
  # load model
  if (attr(plpResult$model, "saveToJson")) {
    modelLocation <-
      reticulate::r_to_py(file.path(plpResult$model$model, "model.json"))
    model <- sklearnFromJson(path = modelLocation)
  } else {
    os <- reticulate::import("os")
    joblib <- reticulate::import("joblib", convert = FALSE)
    modelLocation <-
      reticulate::r_to_py(file.path(plpResult$model$model, "model.pkl"))
    model <- joblib$load(os$path$join(modelLocation))
  }
  
  X <- reticulate::r_to_py(matrixData)
  y <- reticulate::r_to_py(matrix(labels$outcomeCount, ncol = 1))
  
  names <- matrixObjects$covariateMap
  names <- names[order(names$columnId),]
  
  featureImportanceList <- list()
  converged <- FALSE
  p <- 0
  metrics <- c()
  
  if (method == "PFI") {
    reticulate::py_require("scikit-learn")
    sklearn <- reticulate::import("sklearn")
    
    while(!converged && p < length(params)) {
      p <- p + 1 # next set of params
      
      res <- sklearn$inspection$permutation_importance(estimator = model,
                                                       X = X$toarray(),
                                                       y = y$ravel(),
                                                       scoring = params[[p]]$scoring,
                                                       n_repeats = as.integer(params[[p]]$n_repeats))
      values <- res$importances_mean
      
      featureImportanceList[[p]] <- values[names$covariateId %in% covariates]
      converged <- checkConvergence(featureImportanceList)
    }
    
  } else if (method == "SHAP") {
    reticulate::py_require("shap")
    reticulate::py_require("numpy")
    
    shap <- reticulate::import("shap")
    np <- reticulate::import("numpy")
    
    while(!converged && p < length(params)) {
      p <- p + 1 # next set of params
      
      N <- nrow(matrixData)
      background_size <- as.integer(ceiling(N * params[[p]]$background_sample))
      X_background <- X[1L:background_size, ] # TODO: extend with random sample 
      
      explainer = shap$Explainer(model = model,
                                 masker = X_background$toarray(),
                                 link=shap$links$identity)
      
      X_np <- np$array(as.matrix(matrixData))
      shap_values = explainer(X_np)
      
      feature_importance <- np$mean(np$abs(shap_values$values), axis=as.integer(0)) # compute average
      
      featureImportanceList[[p]] <- feature_importance[,1] # TODO: check columns
      converged <- checkConvergence(featureImportanceList)
    }
  }
  
  # create plot when converged or last params done
  if (convergence) {
    convergencePlot(featureImportanceList, params)
  }
  
  varImp <- data.frame(
    covariateId = covariates,
    fi = featureImportanceList[[length(featureImportanceList)]]
  )
  
  return(varImp)
}


checkConvergence <- function(featureImportanceList, metric = "distance", stop = 0.005) {
  
  if (length(featureImportanceList) > 1) {
    # compute pairwise similarity of most recent FIs
    fi_prev <- featureImportanceList[[length(featureImportanceList)-1]]
    fi_curr <- featureImportanceList[[length(featureImportanceList)]]
    
    value <- calculateConvergence(fi_prev, fi_curr, metric)
    
    # check stop criterion
    if ((metric == "correlation" && value >= stop) ||
        (metric %in% c("rmse", "max_diff", "distance") && value <= stop)) {
      ParallelLogger::logInfo("Converged at params ", length(featureImportanceList), " (value = ", value, ")")
      return(TRUE)
    }
  } 
  
  return(FALSE)
}

calculateConvergence <- function(fi_prev, fi_curr, metric = "distance") {
  if (metric == "correlation") {
    value <- cor(fi_prev, fi_curr)
  } else if (metric == "rmse") {
    value <- sqrt(mean((fi_prev - fi_curr)^2))
  } else if (metric == "max_diff") {
    value <- max(abs(fi_prev - fi_curr))
  } else if (metric == "distance") {
    value <- sqrt(sum((fi_prev-fi_curr)^2))
  } else {
    stop("Metric unknown")
  }
  return(as.numeric(value))
}


convergencePlot <- function(featureImportanceList, params, metric = "distance", stop = NULL) {
  
  if (length(featureImportanceList) > 1) {
    
    # Compute for all consecutive pairs
    values <- numeric(length(featureImportanceList)-1)
    for (i in 2:length(featureImportanceList)) {
      v <- calculateConvergence(featureImportanceList[[i-1]], featureImportanceList[[i]], metric)
      # print(v, digits=10)   # Print the raw value with higher precision
      values[i-1] <- v
    }
    
    # convert to df for plot
    param_matrix <- do.call(rbind, lapply(params, function(x) unlist(x)))
    df_plot <- as.data.frame(param_matrix)
    df_plot <- df_plot[2:(length(values)+1),]
    df_plot$value <- values
    
    # TODO: remove packages here
    library(ggplot2)
    library(dplyr)
    
    # 1 plot per setting
    settings <- colnames(df_plot)[-length(colnames(df_plot))]
    for (s in settings) {
      df_plot_s <- df_plot %>%
        group_by(!!sym(s)) %>%
        summarise(avg_value = mean(value, na.rm = TRUE))
      
      plot <- ggplot2::ggplot(df_plot_s, aes(x = !!sym(s), y = avg_value)) +
        geom_line(group=1) +
        geom_point() +
        labs(
          x = paste0("Parameter ", s),
          y = paste("FI ", metric),
          title = "Convergence Analysis"
        ) +
        theme_minimal() + geom_hline(yintercept = stop, col = "red", linetype = 2)
      
      ggsave(paste0("plot_", s, "_", metric, ".png"), plot = plot)
    }
  } else {
    ParallelLogger::logInfo("Need at least two sets of calculated FI.")
  }
}

#' Permutation Feature Importance
#'
#' @description
#' Calculate the permutation feature importance (pfi) for a PLP model.
#' @details
#' The function permutes the each covariate/features \code{repeats} times and
#' calculates the mean AUC change caused by the permutation.
#' @param plpResult                         An object of type \code{runPlp}
#' @param population                       The population created using createStudyPopulation() who will have their risks predicted
#' @param plpData                          An object of type \code{plpData} - the patient level prediction
#'                                         data extracted from the CDM.
#' @param repeats                          The number of times to permute each covariate
#' @param covariates                       A vector of covariates to calculate the pfi for.  If NULL it uses all covariates included in the model.
#' @param cores                            Number of cores to use when running this (it runs in parallel)
#' @param log                              A location to save the log for running pfi
#' @param logthreshold                     The log threshold (e.g., INFO, TRACE, ...)
#'
#' @return
#' A dataframe with the covariateIds and the pfi (change in AUC caused by permuting the covariate) value
#' @examples
#' \donttest{ \dontshow{ # takes too long }
#' library(dplyr)
#' # simulate some data
#' data("simulationProfile")
#' plpData <- simulatePlpData(simulationProfile, n = 1000, seed = 42)
#' # now fit a model
#' saveLoc <- file.path(tempdir(), "pfi")
#' plpResult <- runPlp(plpData, outcomeId = 3, saveDirectory = saveLoc)
#' population <- createStudyPopulation(plpData, outcomeId = 3)
#  # calculate permutation feature importance for nonzero coefficients
#' pfi(plpResult, population, plpData, repeats = 1, cores = 1)
#' # compare to model coefficients
#' plpResult$model$covariateImportance %>% dplyr::filter(.data$covariateValue != 0)
#' # clean up
#' unlink(saveLoc, recursive = TRUE)
#' }
#' @export
pfi <- function(plpResult, population, plpData, repeats = 1,
                covariates = NULL, cores = NULL, log = NULL,
                logthreshold = "INFO") {
  if (!is.null(log)) {
    appender <- ParallelLogger::createFileAppender(
      layout = ParallelLogger::layoutParallel,
      fileName = log
    )
    
    logger <- ParallelLogger::createLogger(
      name = "PAR",
      threshold = logthreshold,
      appenders = list(appender)
    )
    ParallelLogger::registerLogger(logger)
  }
  
  
  if (is.null(covariates)) { # TODO: why is this restricted to positive covariateImportance? 
    covariates <- plpResult$model$covariateImportance %>%
      # dplyr::filter(.data$covariateValue != 0) %>%
      dplyr::select("covariateId") %>%
      dplyr::pull()
  }
  
  # add code to format covariateData based on plpModel
  
  if (!is.null(plpResult$model$preprocessing$featureEngineering)) {
    # do feature engineering/selection
    ParallelLogger::logInfo("Running FE in model")
    plpData <- do.call(
      applyFeatureEngineering,
      list(
        plpData = plpData,
        settings = plpResult$model$preprocessing$featureEngineering
      )
    )
  } else {
    ParallelLogger::logInfo("No FE in model")
  }
  
  if (!is.null(plpResult$model$preprocessing$tidyCovariates)) {
    # do preprocessing
    ParallelLogger::logInfo("Applying data tidying done in model")
    plpData$covariateData <- do.call(
      applyTidyCovariateData,
      list(
        covariateData = plpData$covariateData,
        preprocessSettings = plpResult$model$preprocessing$tidyCovariates
      )
    )
  } else {
    ParallelLogger::logInfo("No data tidying done in model")
  }
  
  # apply prediction function
  pred <- do.call(
    attr(plpResult$model, "predictionFunction"),
    list(
      plpModel = plpResult$model,
      data = plpData,
      cohort = population
    )
  )
  
  auc <- computeAuc(pred)
  
  # do permulation and savePlpData to temp loc
  plpDataLocation <- file.path(tempdir(), paste0("data"))
  savePlpData(plpData, file = plpDataLocation)
  
  if (is.null(cores)) {
    rlang::check_installed("parallel")
    ParallelLogger::logInfo(paste0("Number of cores not specified"))
    cores <- parallel::detectCores()
    ParallelLogger::logInfo(paste0("Using all ", cores))
    ParallelLogger::logInfo(paste0("Set cores input to use fewer..."))
  }
  getVpiSettings <- function(i) {
    result <- list(
      plpModel = plpResult$model,
      population = population,
      plpDataLocation = plpDataLocation,
      covariateId = covariates[i],
      repeats = repeats
    )
    return(result)
  }
  if (cores > 1) {
    cluster <- ParallelLogger::makeCluster(numberOfThreads = cores)
    ParallelLogger::clusterRequire(cluster, c("PatientLevelPrediction", "Andromeda"))
    
    vpiSettings <- lapply(1:length(covariates), getVpiSettings)
    
    aucP <- ParallelLogger::clusterApply(
      cluster = cluster,
      x = vpiSettings,
      fun = permutePerf,
      stopOnError = FALSE,
      progressBar = TRUE
    )
    ParallelLogger::stopCluster(cluster)
  } else {
    ParallelLogger::logInfo("Running in serial")
    aucP <- lapply(1:length(covariates), function(i) {
      permutePerf(getVpiSettings(i))
    })
  }
  aucP <- do.call(c, aucP)
  varImp <- data.frame(
    covariateId = covariates,
    pfi = auc - aucP
  )
  return(varImp)
}


# function to take plpModel,plpData location, load data and permute then calcualte auc
permutePerf <- function(settings) {
  auc <- c()
  for (i in 1:settings$repeats) {
    ParallelLogger::logInfo(paste0("Starting to permute data for covariate: ", settings$covariateId))
    plpData <- tryCatch(
      {
        suppressWarnings(permute(
          settings$plpDataLocation, settings$covariateId,
          settings$population
        ))
      },
      warning = function(war) {
        ParallelLogger::logInfo(paste0("a warning: ", war))
      },
      error = function(err) {
        ParallelLogger::logError(paste0("an error: ", err))
        return(NULL)
      }
    )
    
    if (is.null(plpData)) {
      ParallelLogger::logInfo(paste0("plpData NULL for covariate: ", settings$covariateId))
      return(0)
    }
    ParallelLogger::logInfo(paste0("Calculating prediction for permutation of covariate: ", settings$covariateId))
    
    # need to stop preprocessing and do it once...
    pred <- do.call(
      attr(settings$plpModel, "predictionFunction"),
      list(
        plpModel = settings$plpModel,
        data = plpData,
        cohort = settings$population
      )
    )
    
    auct <- computeAuc(pred)
    auc <- c(auc, auct)
  }
  return(mean(auc))
}


permute <- function(plpDataLocation, cId, population) {
  plpData <- suppressWarnings(PatientLevelPrediction::loadPlpData(plpDataLocation))
  
  # get analysisId
  aId <- plpData$covariateData$covariateRef %>%
    dplyr::filter(.data$covariateId == !!cId) %>%
    dplyr::select("analysisId") %>%
    dplyr::collect()
  
  # if analysis id is not 3 (age group), 4 (race) or 5 (ethnicity)
  if (!aId$analysisId %in% c(3, 4, 5)) {
    # select covariateId data
    coi <- plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId == !!cId) %>%
      dplyr::collect()
    nSamp <- length(coi$rowId)
    
    # find a new random selection of people and give them the covariate and value
    newPlp <- sample(population$rowId, nSamp)
    newData <- dplyr::as_tibble(cbind(rowId = newPlp, coi[, -1]))
    
    # swap old covariate data with new
    plpData$covariateData$covariates <- plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId != !!cId) %>%
      dplyr::collect()
    Andromeda::appendToTable(plpData$covariateData$covariates, newData)
  } else {
    # do some code for the connected variables... with more than 2 options
    
    # get the ppl with covariates
    haveCidData <- plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId == !!cId) %>%
      dplyr::collect()
    nSamp <- length(haveCidData$rowId)
    
    # sample the pop to replace
    swapPlp <- sample(population$rowId, nSamp)
    haveCidDataSwapped <- dplyr::as_tibble(cbind(rowId = swapPlp, haveCidData[, -1]))
    
    # find the swapped people to switch
    connectedCovs <- plpData$covariateData$covariateRef %>%
      dplyr::filter(.data$analysisId == !!aId$analysisId) %>%
      dplyr::group_by(.data$covariateId) %>%
      dplyr::select("covariateId") %>%
      dplyr::collect()
    plpToSwap <- plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId %in% !!connectedCovs$covariateId) %>%
      dplyr::filter(.data$rowId %in% swapPlp) %>%
      dplyr::collect()
    
    swappedForCid <- dplyr::as_tibble(cbind(rowId = haveCidData$rowId[1:nrow(plpToSwap)], plpToSwap[, -1]))
    
    
    # swap old covariate data with new
    plpData$covariateData$covariates <- plpData$covariateData$covariates %>%
      dplyr::filter(.data$covariateId != !!cId) %>%
      dplyr::filter(!.data$rowId %in% !!swapPlp) %>%
      dplyr::collect()
    Andromeda::appendToTable(
      plpData$covariateData$covariates,
      rbind(haveCidDataSwapped, swappedForCid)
    )
  }
  
  return(plpData)
}
