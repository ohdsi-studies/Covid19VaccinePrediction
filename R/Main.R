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
#' @param cdmDatabaseSchemas    Schema name where your patient-level data in OMOP CDM format resides.
#'                             Note that for SQL Server, this should include both the database and
#'                             schema name, for example 'cdm_data.dbo'.
#' @param cdmDatabaseNames      Shareable name of the database 
#' @param cohortDatabaseSchemas Schema name where intermediate data can be stored. You will need to have
#'                             write priviliges in this schema. Note that for SQL Server, this should
#'                             include both the database and schema name, for example 'cdm_data.dbo'.
#' @param cohortTables          The name of the table that will be created in the work database schema.
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
#' @param runValidate          Run the model validation
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
#'         cdmDatabaseSchemas = "cdm_data",
#'         cdmDatabaseNames = 'shareable name of the database'
#'         cohortDatabaseSchemas = "study_results",
#'         cohortTables = "cohort",
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
                    cdmDatabaseSchemas,
                    cdmDatabaseNames = 'friendly database name',
                    cohortDatabaseSchemas = cdmDatabaseSchema,
                    cohortTables = "cohort",
                    oracleTempSchema = cohortDatabaseSchema,
                    outputFolder,
                    createCohorts = F,
                    fetchData = F,
                    runDevelopment = F,
                    runValidate = F,
                    sampleSize = NULL,
                    verbosity = "INFO",
                    cdmVersion = 5) {
  
  if (!file.exists(outputFolder))
    dir.create(outputFolder, recursive = TRUE)
  
  ParallelLogger::addDefaultFileLogger(file.path(outputFolder, "log.txt"))
  
  
  if (createCohorts) {
    for(i in 1:length(cdmDatabaseNames)){
    ParallelLogger::logInfo(paste0("Creating cohorts for ",cdmDatabaseNames[i]))
      createConnectionDetails <- DatabaseConnector::createConnectionDetails
      connectionDetailsT <- do.call('createConnectionDetails',connectionDetails[[i]])
    createCohorts(connectionDetails = connectionDetailsT,
                  cdmDatabaseSchema = cdmDatabaseSchemas[i],
                  cohortDatabaseSchema = cohortDatabaseSchemas[i],
                  cohortTable = cohortTables[i],
                  oracleTempSchema = oracleTempSchema,
                  outputFolder = file.path(outputFolder,cdmDatabaseNames[i])
                  )
    }
  }
  
  if(fetchData){
    
    for(i in 1:length(cdmDatabaseNames)){
      databaseName <- cdmDatabaseNames[i]
      ParallelLogger::logInfo(paste0("Fetching data for ",databaseName))
      
      covSet <- FeatureExtraction::createCovariateSettings(useDemographicsGender = T, 
                                                           useDemographicsAgeGroup  = T, 
                                                           useConditionGroupEraLongTerm = T, 
                                                           useDrugGroupEraLongTerm = T, 
                                                           longTermStartDays = -365, 
                                                           endDays = -1)
      for(targetId in c(1001,2001,3001)){
        fileName <- file.path(outputFolder,databaseName,'data', paste0('T_',targetId))
        fileName2 <- file.path(outputFolder,databaseName,'data', paste0('T_',targetId,'_ag'))
        if(!dir.exists(fileName )){
          dir.create(fileName, recursive = T)
        }
        
        if(!file.exists(file.path(fileName,'covariates'))){
          ParallelLogger::logInfo(paste0('Extracting for ',targetId))
          createConnectionDetails <- DatabaseConnector::createConnectionDetails
          connectionDetailsT <- do.call('createConnectionDetails',connectionDetails[[i]])
          plpData <- tryCatch({PatientLevelPrediction::getPlpData(connectionDetails = connectionDetailsT, 
                                                        cdmDatabaseSchema = cdmDatabaseSchemas[i], 
                                                        oracleTempSchema = oracleTempSchema, 
                                                        cohortId = targetId, 
                                                        outcomeIds = c(1002,2002,3002,1003,2003,3003,1004,2004,3004,1005,2005,3005,1006,2006,3006,1007,2007,3007,4007,1008,2008,3008,1009,2009,4009,1010,2010,3010), 
                                                        cohortDatabaseSchema = cohortDatabaseSchemas[i], 
                                                        cohortTable = cohortTables[i], 
                                                        outcomeDatabaseSchema = cohortDatabaseSchemas[i], 
                                                        outcomeTable = cohortTables[i], 
                                                        cdmVersion = cdmVersion, 
                                                        covariateSettings = covSet, 
                                                        washoutPeriod = 365, 
                                                        firstExposureOnly = T, 
                                                        sampleSize = 2000000)},
                              error = function(e){ParallelLogger::logInfo(e);return(NULL)})
          
          if(!is.null(plpData)){
            PatientLevelPrediction::savePlpData(plpData, fileName)
          }
        } else{
          ParallelLogger::logInfo(paste0('Data exists for ',targetId))
        }
        
        
        # restrict to age/gender only data
        if(!file.exists(file.path(fileName2,'covariates')) & file.exists(file.path(fileName,'covariates')) ){
          ParallelLogger::logInfo(paste0('Creating age/gender only data for ',targetId))
          plpData <- PatientLevelPrediction::loadPlpData(fileName)
          
          plpData$covariateData$covariateIO <- data.frame(covariateId = c(8532001, 8507001, (1:24)*1000+3))
          plpData$covariateData$covariates <- plpData$covariateData$covariates %>% 
            dplyr::inner_join(plpData$covariateData$covariateIO, by = 'covariateId')
          
          plpData$covariateData$covariateRef <- plpData$covariateData$covariateRef %>% 
            dplyr::inner_join(plpData$covariateData$covariateIO, by = 'covariateId')
          
          PatientLevelPrediction::savePlpData(plpData, fileName2)
          
        } else{
          ParallelLogger::logInfo(paste0('Age/gender data exists for ',targetId))
        }
        
      }
    }
    
  }
  
  if(runDevelopment){
    
    for(i in 1:length(cdmDatabaseNames)){
      ParallelLogger::logInfo(paste0("Running predictions for ",cdmDatabaseNames[i]))
      
      targets <- c(1001,2001,3001)
      #outcomes <- list(c(1002,2002,3002),
      #                 c(1003,2003,3003),
      #                 c(1004,2004,3004),
      #                 c(1005,2005,3005),
      #                 c(1006,2006,3006),
      #                 c(1007,2007,3007, 4007),
      #                 c(1008,2008,3008),
      #                 c(1009,2009,4009),
      #                 c(1010,2010,3010))
      
      outcomes <- list(c(1002,2002,3002,1003,2003,3003),
                       c(1004,2004,3004,1005,2005,3005),
                       c(1006,2006,3006,1007,2007,3007),
                       c(1008,2008,3008,1009,2009,4009),
                       c(1010,2010,3010, 4007))
      
      ParallelLogger::logInfo("Full Models")
      developInPar(targets = targets,
                   outcomes = outcomes,
                   ageGenderOnly = F,
                   outputFolder = outputFolder, 
                   databaseName = cdmDatabaseNames[i],
                   seed = 1022 )
      
      ParallelLogger::logInfo("Benchmark Models")
      developInPar(targets = targets,
                   outcomes = outcomes,
                   ageGenderOnly = T,
                   outputFolder = outputFolder, 
                   databaseName = cdmDatabaseNames[i],
                   seed = 1022 )
    }
    
  }
  
  if(runValidate){
 
    for(i in 1:length(cdmDatabaseNames)){
      ParallelLogger::logInfo(paste0("Running validations for ",cdmDatabaseNames[i]))
    # internal other Ts and other Os (create settings per models and valdiate in para)
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
    internalValidateWrapper(resultLocation = outputFolder, 
                            devDatabaseName = cdmDatabaseNames[i],
                            phenotypeGroup = outcomes,
                            targets = targets)
    
    # external: all Ts and all phenos (need data location )
    externalValidateWrapper(resultLocation = outputFolder, 
                            devDatabaseName = cdmDatabaseNames[i],
                            valDatabaseNames = cdmDatabaseNames[-i],
                            phenotypeGroup = outcomes,
                            targets = targets)
      }
    }
  
  invisible(NULL)
}




