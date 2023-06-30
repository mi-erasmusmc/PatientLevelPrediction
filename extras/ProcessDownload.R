ParallelLogger::logInfo(paste0('Starting R Script...'))

library(PatientLevelPrediction)
library(DatabaseConnector)
library(jsonlite)
library(devtools)
library(jsonlite)
library(zip)

Sys.setenv("DATABASECONNECTOR_JAR_FOLDER" = "C:/[ DEV ]/[ DRIVER ]")
Sys.setenv("PROJECT_FOLDER" = "C:/[ DEV ]/PredictionLibraryApp/PredictionLibraryR")

modelId = 1
myArgs = commandArgs(trailingOnly=TRUE)
myDirectory = gsub("\\\\", "/", myArgs)
hasDirectory <- !rlang::is_empty(myDirectory)
resetDatabase = F
saveLoc <- NULL

if (hasDirectory) {
  ParallelLogger::logInfo(paste0("Exporting data to: ", gsub("\\\\", "/", myDirectory), sep="", collapse=NULL))
  setwd(myDirectory)
} else {
  saveLoc <- paste0(Sys.getenv("PROJECT_FOLDER"), "/download/", modelId)
  dir.create(saveLoc)
  setwd(saveLoc)
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

connectionDetails  <- DatabaseConnector::createConnectionDetails(
  dbms = "postgresql",
  user = "postgres",
  password = "12345",
  pathToDriver = Sys.getenv("DATABASECONNECTOR_JAR_FOLDER"),
  server = "127.0.0.1/postgres"
)
con <- DatabaseConnector::connect(connectionDetails = connectionDetails)

getPlpDataSettings <- list(
  databaseDetails = con
)

plpDataSettings <- list(
  databaseDetails = connectionDetails,
  covariateSettings = ParallelLogger::convertJsonToSettings(dataSettings$covariateSettings[i]),
  restrictPlpDataSettings = ParallelLogger::convertJsonToSettings(dataSettings$restrictPlpDataSettings[i])
)


plpData <- tryCatch({
  do.call(getPlpData, databaseDetails = list(targetId = 1, outcomeIds = NULL))
},
error = function(e){ParallelLogger::logError(e); return(NULL)}
)

# 
# 
# savePlpShareable(result = plpResult, saveDirectory = file.path(saveLoc))
# saveZip <- paste0(saveLoc, "/", "model-",gsub("^.*/", "", saveLoc), ".zip")
# setwd(saveDirectory)
# 
# zipFiles <- list.files(path = saveDirectory, pattern = ".")
# zip(zipfile = saveZip, files = zipFiles)
# ParallelLogger::logInfo(paste0('Zip saved to: ', saveZip))

