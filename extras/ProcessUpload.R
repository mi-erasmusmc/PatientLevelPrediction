ParallelLogger::logInfo(paste0('Starting R Script...'))

library(PatientLevelPrediction)
library(DatabaseConnector)
library(jsonlite)
library(devtools)
library(jsonlite)
library(zip)

Sys.setenv("DATABASECONNECTOR_JAR_FOLDER" = "C:/[ DEV ]/[ DRIVER ]")
Sys.setenv("PROJECT_FOLDER" = "C:/[ DEV ]/PredictionLibraryApp/PredictionLibraryR")

myArgs = commandArgs(trailingOnly=TRUE)
myDirectory = gsub("\\\\", "/", myArgs)
hasDirectory <- !rlang::is_empty(myDirectory)
resetDatabase = F
saveLoc <- NULL

if (hasDirectory) {
  ParallelLogger::logInfo(paste0("Importing data from: ", gsub("\\\\", "/", myDirectory), sep="", collapse=NULL))
  setwd(myDirectory)
} else {
  saveLoc <- paste0(Sys.getenv("PROJECT_FOLDER"), "/upload/", gsub('[.]', '', format(Sys.time(),format="%Y%m%d%H%M%S%OS")))
  dir.create(saveLoc)
  setwd(saveLoc)
}

randVar <- rawToChar(as.raw(sample(c(65:90,97:122), 5, replace=T)))
appendRandom <- function(x, rand = randVar){
  return(paste(rand, x, sep=''))
}

randString <- function(characters=0, numbers=0, symbols=0, lowerCase=0, upperCase=0) {
  ASCII <- NULL
  if(symbols>0)    ASCII <- c(ASCII, sample(c(33:47, 58:34, 91:96, 123:126), symbols))
  if(numbers>0)    ASCII <- c(ASCII, sample(48:57, numbers))
  if(upperCase>0)  ASCII <- c(ASCII, sample(65:90, upperCase))
  if(lowerCase>0)  ASCII <- c(ASCII, sample(97:122, lowerCase))
  if(characters>0) ASCII <- c(ASCII, sample(c(65:90, 97:122), characters))
  
  return( rawToChar(as.raw(sample(ASCII, length(ASCII)))) )
}

is_installed <- function(pkg, version = 0) {
  installed_version <- tryCatch(utils::packageVersion(pkg),
                                error = function(e) NA
  )
  !is.na(installed_version) && installed_version >= version
}

ensure_installed <- function(pkg) {
  if (!is_installed(pkg)) {
    msg <- paste0(sQuote(pkg), " must be installed for this functionality.")
    if (interactive()) {
      inform(paste(msg, "Would you like to install it?", sep = "\n"))
      if (menu(c("Yes", "No")) == 1) {
        install.packages(pkg)
      } else {
        stop(msg, call. = FALSE)
      }
    } else {
      stop(msg, call. = FALSE)
    }
  }
}

# Download the PostreSQL driver ---------------------------
# If DATABASECONNECTOR_JAR_FOLDER exists, assume driver has been downloaded
jarFolder <- Sys.getenv("DATABASECONNECTOR_JAR_FOLDER", unset = "")
if (jarFolder == "") {
  tempJarFolder <- tempfile("jdbcDrivers")
  dir.create(tempJarFolder)
  Sys.setenv("DATABASECONNECTOR_JAR_FOLDER" = tempJarFolder)
  downloadJdbcDrivers("postgresql")
  
  withr::defer({
    unlink(tempJarFolder, recursive = TRUE, force = TRUE)
    Sys.unsetenv("DATABASECONNECTOR_JAR_FOLDER")
  }, testthat::teardown_env())
}

travis <- T

if(ifelse(is.null(Sys.info()), T, Sys.info()['sysname'] != 'Windows')){
  # configure and activate python
  PatientLevelPrediction::configurePython(envname = 'r-reticulate', envtype = "conda")
  PatientLevelPrediction::setPythonEnvironment(envname = 'r-reticulate', envtype = "conda")
  
  # if mac install nomkl -- trying to fix github actions
  if(ifelse(is.null(Sys.info()), F, Sys.info()['sysname'] == 'Darwin')){
    reticulate::conda_install(envname = 'r-reticulate', packages = c('nomkl'),
                              forge = TRUE, pip = FALSE, pip_ignore_installed = TRUE,
                              conda = "auto")
  }
}

if (!hasDirectory) {
  plpResult <- NULL
  plpDataSimulationProfile <- NULL
  data(plpDataSimulationProfile, envir = environment())
  
  sampleSize <- 2500+sample(1000,1)
  plpData <- PatientLevelPrediction:::simulatePlpData(plpDataSimulationProfile, n = sampleSize)
  plpData$metaData$cohortId <- plpData$metaData$databaseDetails$cohortId
  
  populationSettings <- PatientLevelPrediction::createStudyPopulationSettings(
    firstExposureOnly = FALSE,
    washoutPeriod = 0,
    removeSubjectsWithPriorOutcome = FALSE,
    priorOutcomeLookback = 99999,
    requireTimeAtRisk = T,
    minTimeAtRisk=10,
    riskWindowStart = 0,
    startAnchor = 'cohort start',
    riskWindowEnd = as.integer(randString(0, 3, 0, 0, 0)),
    endAnchor = 'cohort start'
  )
  
  lrSet <- setLassoLogisticRegression()
  plpResult <- runPlp(
    plpData = plpData,
    outcomeId = plpData$outcomes$outcomeId[1],
    analysisId = paste(Sys.Date(), plpData$outcomes$outcomeId[1], sep = "-"),
    analysisName = 'LASSO LR Testing analysis',
    populationSettings = populationSettings,
    splitSettings = createDefaultSplitSetting(),
    sampleSettings = createSampleSettings(),
    featureEngineeringSettings = createFeatureEngineeringSettings(),
    preprocessSettings = createPreprocessSettings(),
    modelSettings = lrSet,
    logSettings = createLogSettings(verbosity = 'TRACE'),
    executeSettings = createDefaultExecuteSettings(),
    saveDirectory = file.path(paste0(unlist(strsplit(saveLoc, '/'))[-length(unlist(strsplit(saveLoc, '/')))], collapse = '/'), 'generated')
  )
}

if (!hasDirectory) {
  saveDirectory = file.path(paste0(unlist(strsplit(saveLoc, '/'))[-length(unlist(strsplit(saveLoc, '/')))], collapse = '/'), 'export')
  savePlpShareable(result = plpResult, saveDirectory = file.path(saveDirectory))
  saveZip <- paste0(saveLoc, "/", "model-",gsub("^.*/", "", saveLoc), ".zip")
  setwd(saveDirectory)

  zipFiles <- list.files(path = saveDirectory, pattern = ".")
  zip(zipfile = saveZip, files = zipFiles)
  ParallelLogger::logInfo(paste0('Zip saved to: ', saveZip))
} else {
  ParallelLogger::logInfo(paste0('Reading model into runPlp object...'))
  
  connectionDetails  <- DatabaseConnector::createConnectionDetails(
    dbms = "postgresql",
    user = "postgres",
    password = "12345",
    pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER"),
    server = "127.0.0.1/postgres"
  )
  con <- DatabaseConnector::connect(connectionDetails = connectionDetails)
  
  createPlpResultTables(
    conn = con,
    resultSchema = "covid_vaccination_plp",
    targetDialect = 'postgresql',
    deleteTables = resetDatabase,
    createTables = resetDatabase
  )
  
  runPlp <- loadPlpShareable(loadDirectory = myDirectory)
  addRunPlpToDatabase(
    databaseList = NULL,
    conn = con,
    runPlp = runPlp,
    modelSaveLocation = myDirectory,
    databaseSchemaSettings = createDatabaseSchemaSettings(
      resultSchema = "covid_vaccination_plp",
      targetDialect = "postgresql",
      tablePrefix = "",
    ),
    cohortDefinitions = data.frame(
      cohortName = c(strtoi(format(Sys.time(),format="%M%S%OS"))+1, strtoi(format(Sys.time(),format="%M%S%OS"))+2, strtoi(format(Sys.time(),format="%M%S%OS"))+3), 
      cohortId = c(1,2,3), 
      json = rep('bla',3)
    )
  )
}