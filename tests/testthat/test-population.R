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


library("testthat")
context("Population")

### Help functions
getPlpData <- function(){
  return(readRDS(paste0(getwd(), "/tests/testthat/plpData_sample.rds")))
}

getDefaultPopulationSettings <- function(outcomeId = 2){
  default_settings <- formals(PatientLevelPrediction::createStudyPopulation)
  default_settings$plpData <- getPlpData()
  default_settings$outcomeId <- outcomeId
  return(as.list(default_settings))
}

default_settings <- getDefaultPopulationSettings()

isvalid_studyPopulation <- function(studyPopulation) {
  expect_is(studyPopulation, "data.frame")
  
  # check if expected columns are present
  expect_equal(sum((c("rowId","subjectId","cohortId","cohortStartDate","daysFromObsStart","daysToCohortEnd","daysToObsEnd","ageYear","gender","outcomeCount","timeAtRisk","daysToEvent","survivalTime") %in% colnames(studyPopulation))==FALSE),0)
  
  # check if row count > 0
  # TODO: this is not strictly required but I think we should stop already at the plp object if that is empty
  expect_gt(nrow(studyPopulation), 0)
  
}

# TODO: add more informative error messages (return actual numbers?)
# includedRowIds = rowIds you know should BE IN the study population (could be more)
# excludedRowIds = rowIds you know should NOT BE IN the study population anymore (could be more)
# rowIdsWithOutcome = rowIds you know should HAVE the outcome (could be more)
# rowIdsWithoutOutcome = rowIds you know should NOT HAVE the outcome (could be more)
iscorrect_studyPopulation <- function(studyPopulation, includedRowIds = c(), excludedRowIds = c(), rowIdsWithOutcome = c(), rowIdsWithoutOutcome = c()) {
  expect_equal(sum((includedRowIds %in% studyPopulation$rowId)==FALSE), 0, info = "includedRowIds NOT in study population")
  expect_equal(sum((excludedRowIds %in% studyPopulation$rowId)==TRUE), 0, info = "excludedRowIds STILL in study population")
  expect_equal(sum(studyPopulation$outcomeCount[studyPopulation$rowId %in% rowIdsWithOutcome] == 0), 0, info = "rowIdsWithOutcome do not have outcome")
  expect_equal(sum(studyPopulation$outcomeCount[studyPopulation$rowId %in% rowIdsWithoutOutcome] == 1), 0, info = "rowIdsWithoutOutcome do have outcome")
}

###  Test functions

# TODO: create weird value test
# test the abnormal input cases eg.) -365 / cohort stard / -25 / cohort ent
test_that("Check the riskWindow parameters: 365/start/365/start", {
  settings <- default_settings
  
  # Set test parameters
  settings$riskWindowStart <- -365
  settings$removeSubjectsWithPriorOutcome <- F
  settings$startAnchor <- 'cohort start'
  settings$riskWindowEnd <- 0
  settings$endAnchor <- 'cohort start'
  
  # Test
  expect_error(studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings))
})

# TODO: not sure if this is necessary in every test, maybe test once in beginning?
# isvalid_studyPopulation(studyPopulation)   

test_that("Patients with outcome but without full time-at-risk are included", {
  settings <- default_settings
  
  # Set test parameters
  settings$includeAllOutcomes <- TRUE
  settings$requireTimeAtRisk <- TRUE
  settings$minTimeAtRisk <- 100
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(1, 3, 4),
                            excludedRowIds = c(2),
                            rowIdsWithOutcome = c(1, 3),
                            rowIdsWithoutOutcome = c(4))
})

test_that("Patients with outcome but without full time-at-risk are excluded", {
  settings <- default_settings
  
  # Set test parameters
  settings$includeAllOutcomes <- FALSE
  settings$requireTimeAtRisk <- TRUE
  settings$minTimeAtRisk <- 100
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(3, 4),
                            excludedRowIds = c(1, 2),
                            rowIdsWithOutcome = c(3),
                            rowIdsWithoutOutcome = c(4))
})

test_that("Patients without outcome but with full time-at-risk are included", {
  settings <- default_settings
  
  # Set test parameters
  settings$includeAllOutcomes <- FALSE
  settings$requireTimeAtRisk <- TRUE
  settings$minTimeAtRisk <- 80
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(1, 2, 3, 4),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(1, 3),
                            rowIdsWithoutOutcome = c(2, 4))
})

test_that("Patients without outcome and without full time-at-risk are excluded", {
  settings <- default_settings
  
  # Set test parameters
  settings$includeAllOutcomes <- FALSE
  settings$requireTimeAtRisk <- TRUE
  settings$minTimeAtRisk <- 100
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(3, 4),
                            excludedRowIds = c(1, 2),
                            rowIdsWithOutcome = c(3),
                            rowIdsWithoutOutcome = c(4))
})


test_that("Check the riskWindow parameters: 0/start/365/start", {
  settings <- default_settings
  
  # Set test parameters
  settings$riskWindowStart <- 0
  settings$startAnchor <- 'cohort start'
  settings$riskWindowEnd <- 365
  settings$endAnchor <- 'cohort start'
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(1,3,6,8,10,12,15,16,17,19,20),
                            rowIdsWithoutOutcome = c())
})

test_that("Check the riskWindow parameters: 0/start/10/start", {
  settings <- default_settings
  
  # Set test parameters
  settings$riskWindowStart <- 0
  settings$startAnchor <- 'cohort start'
  settings$riskWindowEnd <- 10
  settings$endAnchor <- 'cohort start'
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(17),
                            rowIdsWithoutOutcome = c())
})

test_that("Check the riskWindow parameters: 0/start/10/end", {
  settings <- default_settings
  
  # Set test parameters
  settings$riskWindowStart <- 0
  settings$startAnchor <- 'cohort start'
  settings$riskWindowEnd <- 10
  settings$endAnchor <- 'cohort end'
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(18),
                            rowIdsWithoutOutcome = c())
})

test_that("Check the riskWindow parameters: 0/start/0/start", {
  settings <- default_settings
  
  # Set test parameters
  settings$riskWindowStart <- 0
  settings$startAnchor <- 'cohort start'
  settings$riskWindowEnd <- 0
  settings$endAnchor <- 'cohort start'
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(19),
                            rowIdsWithoutOutcome = c())
})

test_that("Check the riskWindow parameters: 365/start/365/start", {
  settings <- default_settings
  
  # Set test parameters
  settings$riskWindowStart <- 365
  settings$startAnchor <- 'cohort start'
  settings$riskWindowEnd <- 365
  settings$endAnchor <- 'cohort start'
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(20),
                            rowIdsWithoutOutcome = c())
})

## first exposure only

testthat::test_that("Testing firstExposureOnly - include patients with more than one exposure", {
  settings <- default_settings
  
  # Set test parameters
  settings$firstExposureOnly <- FALSE # Default
  
  # Test
  studyPopulation <- do.call(createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(5, 6),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(6),
                            rowIdsWithoutOutcome = c(5))
}) 

testthat::test_that("Testing firstExposureOnly - exclude patients with more than one exposure", {
  settings <- default_settings
  
  # Set test parameters
  settings$firstExposureOnly <- TRUE 
  
  # Test
  studyPopulation <- do.call(createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(5),
                            excludedRowIds = c(6),
                            rowIdsWithOutcome = c(),
                            rowIdsWithoutOutcome = c(5))
}) 

testthat::test_that("Testing washout period - patients with prior observation period of at least 364 days", {
  settings <- default_settings
  
  # Set test parameters
  settings$washoutPeriod <- 364
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(7, 8, 9, 10, 11, 12),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(8, 10, 12),
                            rowIdsWithoutOutcome = c(7, 9, 11))
}) 

testthat::test_that("Testing washout period - patients with prior observation period of at least 365 days", {
  settings <- default_settings
  
  # Set test parameters
  settings$washoutPeriod <- 365
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(7, 8, 11, 12),
                            excludedRowIds = c(9, 10),
                            rowIdsWithOutcome = c(8, 12),
                            rowIdsWithoutOutcome = c(7, 11))
})

testthat::test_that("Testing washout period - patients with prior observation period of at least 366 days", {
  settings <- default_settings
  
  # Set test parameters
  settings$washoutPeriod <- 366
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(11, 12),
                            excludedRowIds = c(7, 8, 9, 10),
                            rowIdsWithOutcome = c(12),
                            rowIdsWithoutOutcome = c(11))
}) 

testthat::test_that("Testing washout period - patients with prior observation period of at least 367 days", {
  settings <- default_settings
  
  # Set test parameters
  settings$washoutPeriod <- 367
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(),
                            excludedRowIds = c(7, 8, 9, 10, 11, 12),
                            rowIdsWithOutcome = c(),
                            rowIdsWithoutOutcome = c())
}) 


test_that("Patients with prior outcome should not be excluded", {
  settings <- default_settings
  
  # Set test parameters
  settings$removeSubjectsWithPriorOutcome <- FALSE
  settings$priorOutcomeLookback <- 99999
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(13, 14, 15, 16),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(15,16),
                            rowIdsWithoutOutcome = c(13,14))
})


test_that("Patients with outcome 10 days before index should be excluded", {
  settings <- default_settings
  
  # Set test parameters
  settings$removeSubjectsWithPriorOutcome <- TRUE
  settings$priorOutcomeLookback <- 10
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(13, 14, 15, 16),
                            excludedRowIds = c(),
                            rowIdsWithOutcome = c(15,16),
                            rowIdsWithoutOutcome = c(13,14))
})


test_that("Patients with outcome 13 days before index should be excluded", {
  settings <- default_settings
  
  # Set test parameters
  settings$removeSubjectsWithPriorOutcome <- TRUE
  settings$priorOutcomeLookback <- 13
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(13, 14, 15),
                            excludedRowIds = c(16),
                            rowIdsWithOutcome = c(15),
                            rowIdsWithoutOutcome = c(13,14))
})



test_that("Patients with outcome 15 days before index should be excluded", {
  settings <- default_settings
  
  # Set test parameters
  settings$removeSubjectsWithPriorOutcome <- TRUE
  settings$priorOutcomeLookback <- 15
  
  # Test
  studyPopulation <- do.call(PatientLevelPrediction::createStudyPopulation, settings)
  iscorrect_studyPopulation(studyPopulation,
                            includedRowIds = c(13, 15),
                            excludedRowIds = c(14, 16),
                            rowIdsWithOutcome = c(15),
                            rowIdsWithoutOutcome = c(13))
})
