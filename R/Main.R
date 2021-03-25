# Copyright 2020 Observational Health Data Sciences and Informatics
#
# This file is part of CovidVaccinePrediction
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

#' Execute the Study
#'
#' @details
#' This function executes the CovidVaccinePrediction Study.
#' 
#' @param connectionDetails    An object of type \code{connectionDetails} as created using the
#'                             \code{\link[DatabaseConnector]{createConnectionDetails}} function in the
#'                             DatabaseConnector package.
#' @param cdmDatabaseSchema    Schema name where your patient-level data in OMOP CDM format resides.
#'                             Note that for SQL Server, this should include both the database and
#'                             schema name, for example 'cdm_data.dbo'.
#' @param cdmDatabaseName      Shareable name of the database 
#' @param cohortDatabaseSchema Schema name where intermediate data can be stored. You will need to have
#'                             write priviliges in this schema. Note that for SQL Server, this should
#'                             include both the database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTable          The name of the table that will be created in the work database schema.
#'                             This table will hold the target population cohorts used in this
#'                             study.
#' @param oracleTempSchema     Should be used in Oracle to specify a schema where the user has write
#'                             priviliges for storing temporary tables.
#' @param outputFolder         Name of local folder to place results; make sure to use forward slashes
#'                             (/). Do not use a folder on a network drive since this greatly impacts
#'                             performance.
#' @param createProtocol       Creates a protocol based on the analyses specification                             
#' @param createCohorts        Create the cohortTable table with the target population and outcome cohorts?
#' @param fetchData        Only fetch data for the analyses without fitting models. Setting this flag will overwrite your input provided to the runAnalyses and createCohorts parameters.
#' @param runDevelopment          Run the model development
#' @param sampleSize           The number of patients in the target cohort to sample (if NULL uses all patients)
#' @param verbosity            Sets the level of the verbosity. If the log level is at or higher in priority than the logger threshold, a message will print. The levels are:
#'                                         \itemize{
#'                                         \item{DEBUG}{Highest verbosity showing all debug statements}
#'                                         \item{TRACE}{Showing information about start and end of steps}
#'                                         \item{INFO}{Show informative information (Default)}
#'                                         \item{WARN}{Show warning messages}
#'                                         \item{ERROR}{Show error messages}
#'                                         \item{FATAL}{Be silent except for fatal errors}
#'                                         }                              
#' @param cdmVersion           The version of the common data model                       
#'
#' @examples
#' \dontrun{
#' connectionDetails <- createConnectionDetails(dbms = "postgresql",
#'                                              user = "joe",
#'                                              password = "secret",
#'                                              server = "myserver")
#'
#' execute(connectionDetails,
#'         cdmDatabaseSchema = "cdm_data",
#'         cdmDatabaseName = 'shareable name of the database'
#'         cohortDatabaseSchema = "study_results",
#'         cohortTable = "cohort",
#'         oracleTempSchema = NULL,
#'         outputFolder = "c:/temp/study_results", 
#'         createCohorts = T,
#'         fetchData = T,
#'         runDevelopment = T,
#'         sampleSize = 10000,
#'         verbosity = "INFO",
#'         cdmVersion = 5)
#' }
#'
#' @export
execute <- function(connectionDetails,
                    cdmDatabaseSchema,
                    cdmDatabaseName = 'friendly database name',
                    cohortDatabaseSchema = cdmDatabaseSchema,
                    cohortTable = "cohort",
                    oracleTempSchema = cohortDatabaseSchema,
                    outputFolder,
                    createCohorts = F,
                    fetchData = F,
                    runDevelopment = F,
                    sampleSize = NULL,
                    verbosity = "INFO",
                    cdmVersion = 5) {
  
  if (!file.exists(outputFolder))
    dir.create(outputFolder, recursive = TRUE)
  
  ParallelLogger::addDefaultFileLogger(file.path(outputFolder, "log.txt"))
  
  
  if (createCohorts) {
    ParallelLogger::logInfo("Creating cohorts")
    createCohorts(connectionDetails = connectionDetails,
                  cdmDatabaseSchema = cdmDatabaseSchema,
                  cohortDatabaseSchema = cohortDatabaseSchema,
                  cohortTable = cohortTable,
                  oracleTempSchema = oracleTempSchema,
                  outputFolder = outputFolder)
  }
  
  if(fetchData){
    ParallelLogger::logInfo("Fetching data")

    covSet <- FeatureExtraction::createCovariateSettings(useDemographicsGender = T, 
                                                         useDemographicsAge = T, 
                                                         useConditionGroupEraLongTerm = T, 
                                                         useDrugGroupEraLongTerm = T, 
                                                         longTermStartDays = -1, 
                                                         endDays = -365)
    for(targetId in c(1001,2001,3001)){
      fileName <- file.path(outputFolder,'data', paste0('T_',targetId))
      if(!dir.exists(fileName )){
        dir.create(fileName, recursive = T)
      }
      
      plpData <- PatientLevelPrediction::getPlpData(connectionDetails = connectionDetails, 
                                         cdmDatabaseSchema = cdmDatabaseSchema, 
                                         oracleTempSchema = oracleTempSchema, 
                                         cohortId = targetId, 
                                         outcomeIds = c(1002,2002,3002,1003,2003,3003,1004,2004,3004,1005,2005,3005,1006,2006,3006,1007,2007,3007,4007,1008,2008,3008,1009,2009,4009,1010,2010,3010), 
                                         cohortDatabaseSchema = cohortDatabaseSchema, 
                                         cohortTable = cohortTable, 
                                         outcomeDatabaseSchema = cohortDatabaseSchema, 
                                         outcomeTable = cohortTable, 
                                         cdmVersion = cdmVersion, 
                                         covariateSettings = covSet, 
                                         washoutPeriod = 365, 
                                         firstExposureOnly = T, 
                                         sampleSize = 1000000)
      
      PatientLevelPrediction::savePlpData(plpData, fileName)
    }
    
  }
  
  if(runDevelopment){
    ParallelLogger::logInfo("Running predictions")
    
    targets <- c(1001,2001,3001)
    outcomes <- list(c(1002,2002,3002),
                     c(1003,2003,3003),
                     c(1004,2004,3004),
                     c(1005,2005,3005),
                     c(1006,2006,3006),
                     c(1007,2007,3007, 4007),
                     c(1008,2008,3008),
                     c(1009,2009,4009),
                     c(1010,2010,3010))
    
    developInPar(targets = targets,
                 outcomes = outcomes,
                 outputFolder = outputFolder, 
                 seed = 1022 )
    
  }
  
  invisible(NULL)
}




