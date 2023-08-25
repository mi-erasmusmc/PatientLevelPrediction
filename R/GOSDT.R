
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
    regularization = list(0.001),
    maxDepth = list(2), # depth_budget
    warmLB = list(TRUE),
    pathToLabels = "/tmp/warm_lb_labels/warm_label.tmp", #
    timeLimit = list(60),
    similarSupport = list(FALSE),
    seed = sample(100000,1)
){
  
  checkIsClass(seed, c('numeric','integer'))
  checkIsClass(regularisation, c('list'))
  checkIsClass(maxDepth, c('list'))
  checkIsClass(warmLB, c('list'))
  checkIsClass(pathToLabels, c('string'))
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
    requiresDenseMatrix = F,
    name = "General Optimal Sparse Decision Tree",
    pythonImport = 'gosdt',
    pythonImportSecond = NULL,
    pythonClassifier = 'model.gosdt.GOSDT'
  ) 
  
  attr(param, 'saveToJson') <- T
  attr(param, 'saveType') <- 'file'
  
  result <- list(
    fitFunction = "fitGOSDT",
    param = param
  )
  class(result) <- "modelSettings"
  
  return(result)
}

# Adapted from SklearnClassifierSettings.R
GOSDTClassifierInputs <- function(classifier, param){
  
  model <- classifier(
    regularization =  param[[which.max(names(param)=='regularization')]],
    maxDepth = param[[which.max(names(param)=='maxDepth')]],
    warmLB = warmLB, # TODO: what to add here?
    pathToLabels = pathToLabels, # TODO: what to add here?
    timeLimit = param[[which.max(names(param)=='timeLimit')]],
    similarSupport = similarSupport, # TODO: what to add here?
    random_state = param[[which.max(names(param)=='seed')]],
    verbose = 0L
  )
  
  return(model)
}


# Adapted from SklearnClassifier.R
fitGOSDT <- function(
    trainData,
    modelSettings,
    search = "grid",
    analysisId,
    ...) { 
    
  param <- modelSettings$param
  
  # check covariate data
  if(!FeatureExtraction::isCovariateData(trainData$covariateData)){stop("Needs correct covariateData")}
  
  # get the settings from the param
  pySettings <- attr(param, 'settings')
  
  # make sure the inputs are valid
  checkPySettings(pySettings)
  
  start <- Sys.time()
  
  if(!is.null(trainData$folds)){
    trainData$labels <- merge(trainData$labels, trainData$fold, by = 'rowId')
  }
  
  # convert the data to a sparse R Matrix and then use reticulate to convert to python sparse
  # need to save the covariateMap so we know the covariateId to columnId when applying model
  mappedData <- toSparseM(trainData)
  
  matrixData <- mappedData$dataMatrix
  labels <- mappedData$labels
  covariateRef <- mappedData$covariateRef
  
  # save the model to outLoc
  outLoc <- createTempModelLoc()
  
  # functions does CV and fits final models
  # returns: prediction (Train/CV),
  #          finalParam (optimal hyper-parameters)
  #          variableImportance (final model) 
  #          paramGridSearch list with performance and params for complete grid search
  # at the moment it uses AUC as performance but this could be modified to let user
  # specify the performance metric to optimise
  cvResult <- do.call( 
    what = gridCvPython,
    args = list(
      matrixData = matrixData,
      labels = labels,
      seed = pySettings$seed,
      requiresDenseMatrix = pySettings$requiresDenseMatrix,
      modelName = pySettings$name,
      pythonImport = pySettings$pythonImport,
      pythonImportSecond = pySettings$pythonImportSecond,
      pythonClassifier = pySettings$pythonClassifier,
      modelLocation = outLoc,
      paramSearch = param,
      saveToJson = attr(param, 'saveToJson')
    )
  )
  
  hyperSummary <- do.call(rbind, lapply(cvResult$paramGridSearch, function(x) x$hyperSummary))
  
  prediction <- cvResult$prediction
  
  variableImportance <- cvResult$variableImportance
  variableImportance[is.na(variableImportance)] <- 0
  
  incs <- rep(1, nrow(covariateRef))
  covariateRef$included <- incs
  covariateRef$covariateValue <- unlist(variableImportance) # check this is correct order
  
  comp <- start - Sys.time()
  
  result <- list(
    model = file.path(outLoc),
    
    preprocessing = list(
      featureEngineering = attr(trainData, "metaData")$featureEngineering,
      tidyCovariates = attr(trainData$covariateData, "metaData")$tidyCovariateDataSettings, 
      requireDenseMatrix = attr(param, 'settings')$requiresDenseMatrix
    ),
    
    prediction = prediction,
    
    modelDesign = PatientLevelPrediction::createModelDesign(
      targetId = attr(trainData, "metaData")$targetId,
      outcomeId = attr(trainData, "metaData")$outcomeId,
      restrictPlpDataSettings = attr(trainData, "metaData")$restrictPlpDataSettings,
      covariateSettings = attr(trainData, "metaData")$covariateSettings,
      populationSettings = attr(trainData, "metaData")$populationSettings,
      featureEngineeringSettings = attr(trainData$covariateData, "metaData")$featureEngineeringSettings,
      preprocessSettings = attr(trainData$covariateData, "metaData")$preprocessSettings,
      modelSettings = modelSettings,
      splitSettings = attr(trainData, "metaData")$splitSettings,
      sampleSettings = attr(trainData, "metaData")$sampleSettings
    ),
    
    trainDetails = list(
      analysisId = analysisId,
      analysisSource = '', #TODO add from model
      developmentDatabase = attr(trainData, "metaData")$cdmDatabaseSchema,
      attrition = attr(trainData, "metaData")$attrition, 
      trainingTime = paste(as.character(abs(comp)), attr(comp,'units')),
      trainingDate = Sys.Date(),
      modelName = pySettings$name, 
      finalModelParameters = cvResult$finalParam,
      hyperParamSearch = hyperSummary
    ),
    
    covariateImportance = covariateRef
  )
  
  class(result) <- "plpModel"
  attr(result, "predictionFunction") <- "predictPythonSklearn"
  attr(result, "modelType") <- "binary"
  attr(result, "saveType") <- attr(param, 'saveType') # in save/load plp
  attr(result, "saveToJson") <- attr(param, 'saveToJson') # when saving in reticulate
  
  return(result)
}


# Adapted from SklearnClassifier.R
predictPythonGOSDT <- function(
    plpModel, 
    data, 
    cohort
){
  
  if(inherits(data, 'plpData')){
    # convert
    matrixObjects <- toSparseM(
      plpData = data, 
      cohort = cohort,
      map = plpModel$covariateImportance %>% 
        dplyr::select("columnId", "covariateId")
    )
    
    newData <- matrixObjects$dataMatrix
    cohort <- matrixObjects$labels
    
  }else{
    newData <- data
  }
  
  # load model
  if(attr(plpModel,'saveToJson')){
    modelLocation <- reticulate::r_to_py(file.path(plpModel$model,"model.json"))
    model <- sklearnFromJson(path=modelLocation)
  } else{
    os <- reticulate::import('os')
    joblib <- reticulate::import('joblib', convert=FALSE)
    modelLocation <- reticulate::r_to_py(file.path(plpModel$model,"model.pkl"))
    model <- joblib$load(os$path$join(modelLocation)) 
  }
  included <- plpModel$covariateImportance$columnId[plpModel$covariateImportance$included>0] # does this include map?
  pythonData <- reticulate::r_to_py(newData[,included, drop = F])
  
  # make dense if needed
  if(plpModel$preprocessing$requireDenseMatrix){
    pythonData <- pythonData$toarray()
  }
  
  cohort <- predictValues(
    model = model, 
    data = pythonData, 
    cohort = cohort, 
    type = attr(plpModel, 'modelType')
  )
  
  return(cohort)
}
