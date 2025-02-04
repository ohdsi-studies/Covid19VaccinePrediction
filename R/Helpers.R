
developInPar <- function(targets,outcomes, ageGenderOnly = F, outputFolder, databaseName, seed , cores = NULL){
  model <- PatientLevelPrediction::setLassoLogisticRegression()
  
  for(targetId in targets ){
      
      ParallelLogger::logInfo(paste0('Running models for: target ', targetId, 
                                     ' and outcomes ', paste0(outcomes, collapse = '', sep=',')))
      
      # create the cluster
      if(is.null(cores)){
        ParallelLogger::logInfo(paste0('Number of cores not specified'))
        cores <- length(outcomes) 
      }
      ParallelLogger::logInfo(paste0('Using this many cores ', cores))
      ParallelLogger::logInfo(paste0('Set cores input to use fewer...'))
      
      
      cluster <- ParallelLogger::makeCluster(numberOfThreads = cores)
      ParallelLogger::clusterRequire(cluster, c("PatientLevelPrediction", "Andromeda"))
      
      logger <- ParallelLogger::createLogger(name = "PAR",
                                             threshold = 'INFO', 
                                             appenders = list(ParallelLogger::createFileAppender(layout = ParallelLogger::layoutParallel, 
                                                                                                 fileName = file.path(outputFolder,'parlog.txt'))))
      ParallelLogger::registerLogger(logger)
      
      # settings:
      getRunSettings <- function(outcomeId){
        result <- list(plpDataLoc = file.path(outputFolder,databaseName,'data',paste0('T_', targetId, ifelse(ageGenderOnly,'_ag','') )),
                       outcomeId = outcomeId,
                       modelSettings = model,
                       splitSettings = PatientLevelPrediction::createDefaultSplitSetting(
                         splitSeed = seed
                         ),
                       populationSettings = PatientLevelPrediction::createStudyPopulationSettings(
                         firstExposureOnly = FALSE,
                         washoutPeriod = 0,
                         removeSubjectsWithPriorOutcome = T, # remove if they have outcome in prior year?
                         priorOutcomeLookback = 99999,
                         requireTimeAtRisk = T,
                         minTimeAtRisk= 1,
                         riskWindowStart = 0,
                         startAnchor = 'cohort start',
                         riskWindowEnd = 365,
                         endAnchor = 'cohort start'
                       ),
                       executeSettings = PatientLevelPrediction::createExecuteSettings(
                         runSplitData = T, 
                         runPreprocessData = T, 
                         runModelDevelopment = T, 
                         runCovariateSummary = F
                         ),
                       preprocessSettings = PatientLevelPrediction::createPreprocessSettings(),
                       saveDirectory= file.path(outputFolder,databaseName,'models'), 
                       analysisId = paste0('T_',targetId,'_O_',outcomeId,ifelse(ageGenderOnly,'_ag','')))
        return(result)
      }
      runSettings <- lapply(outcomes, getRunSettings)
      
      allResults <- ParallelLogger::clusterApply(cluster = cluster, 
                                                 x = runSettings, 
                                                 fun = runPlpI, 
                                                 stopOnError = FALSE,
                                                 progressBar = TRUE)
      ParallelLogger::stopCluster(cluster)
      
      # garbage collection?
      gc(verbose=TRUE, full=TRUE)
      
    }
    
}

runPlpI <- function(settings){
  
  if(dir.exists(file.path(settings$saveDirectory,settings$analysisId, 'plpResult'))){
    ParallelLogger::logInfo(paste0('Results exist in ',file.path(settings$saveDirectory,settings$analysisId), ' so skipping'))
    return(NULL)
  }
  if(!dir.exists(file.path(settings$plpDataLoc))){
    ParallelLogger::logInfo(paste0('No Data in ',file.path(settings$plpDataLoc), ' so skipping'))
    return(NULL)
  }
  
  settings$plpData <- tryCatch({PatientLevelPrediction::loadPlpData(settings$plpDataLoc)},
                               error = function(e){ParallelLogger::logInfo(e); return(NULL)})
  settings$plpDataLoc <- NULL
  
  result <- tryCatch({do.call(PatientLevelPrediction::runPlp, settings)},
                     error = function(e){ParallelLogger::logInfo(e); return(NULL)})
  return(result)
}



internalValidateWrapper <- function(resultLocation, devDatabaseName,phenotypeGroup,targets){
  
  
  for(phenotypes in phenotypeGroup){
    for(phenotype in phenotypes){
      for(target in targets){
        
        internalValidateInPar(resultLocation = resultLocation, 
                              devDatabaseName = devDatabaseName, 
                              phenotypes = phenotypes, 
                              targets = targets, 
                              phenotype = phenotype, 
                              target = target, 
                              benchmark = F)
      }
    }
  }
  
  
  
}

internalValidateInPar <- function(resultLocation, devDatabaseName, phenotypes, targets, phenotype, target, benchmark){
  # internal other Ts and other Os (create settings per models and valdiate in para)
  modelLoc <- file.path(resultLocation,devDatabaseName,'models', paste0('T_',target,'_O_',phenotype,ifelse(benchmark,'_ag', '')), 'plpResult')
  
  plpResult <- tryCatch({PatientLevelPrediction::loadPlpResult(modelLoc)}, 
                        error = function(e){ParallelLogger::logInfo(e);return(NULL)})
  
  if(!is.null(plpResult)){
    
    if(!dir.exists(file.path(resultLocation,devDatabaseName,'validation'))){
      dir.create(file.path(resultLocation,devDatabaseName,'validation'), recursive = T)
    }
    
    dataToValidate <- expand.grid(devDatabaseName,phenotypes,targets)
    colnames(dataToValidate) <- c('databaseName','outcomeId','T')
    
    # remove model data:
    ind <- dataToValidate$outcomeId == phenotype & dataToValidate$T == target
    dataToValidate <- dataToValidate[!ind,]
    
    dataToValidate$plpDataLoc <- file.path(resultLocation,dataToValidate$databaseName,'data', paste0('T_',dataToValidate$T,ifelse(benchmark,'_ag', '')))
    dataToValidate$modelLoc <- modelLoc
    dataToValidate$saveDirectory <- file.path(resultLocation,devDatabaseName,'validation',dataToValidate$databaseName)
    dataToValidate$analysisId <- paste0('Td_',target,'_Od_',phenotype,ifelse(benchmark,'_ag', ''),'_Tv_',dataToValidate$T,'_Ov_',dataToValidate$outcomeId,ifelse(benchmark,'_ag', ''))
    
    ParallelLogger::logInfo(paste0('Internally validating model across phenotypes ', modelLoc ,' ', nrow(dataToValidate), ' times'))
    
    # create the cluster
    ParallelLogger::logInfo(paste0('Number of cores not specified'))
    cores <- min(nrow(dataToValidate), parallel::detectCores() - 3)## set to max 
    ParallelLogger::logInfo(paste0('Using this many cores ', cores))
    ParallelLogger::logInfo(paste0('Set cores input to use fewer...'))
    
    
    cluster <- ParallelLogger::makeCluster(numberOfThreads = cores)
    ParallelLogger::clusterRequire(cluster, c("PatientLevelPrediction", "Andromeda"))
    
    logger <- ParallelLogger::createLogger(name = "PAR",
                                           threshold = 'INFO', 
                                           appenders = list(ParallelLogger::createFileAppender(layout = ParallelLogger::layoutParallel, 
                                                                                               fileName = file.path(resultLocation,devDatabaseName,'validation','parlog.txt'))))
    ParallelLogger::registerLogger(logger)
    
    # settings: convert dataToValidate to list
    dataToValidateList <- split(dataToValidate, 1:nrow(dataToValidate)) #!!!
    
    allResults <- ParallelLogger::clusterApply(cluster = cluster, 
                                               x = dataToValidateList, 
                                               fun = valPlpI, 
                                               stopOnError = FALSE,
                                               progressBar = TRUE)
    ParallelLogger::stopCluster(cluster)
    
  }
  
  
}


externalValidateWrapper <- function(resultLocation, 
                                    devDatabaseName, 
                                    valDatabaseNames,
                                    phenotypeGroup, 
                                    targets){
  
  
  for(phenotypes in phenotypeGroup){
    for(phenotype in phenotypes){
      for(target in targets){
        
        ParallelLogger::logInfo('Eternally validating full model')
        externalValidateInPar(resultLocation = resultLocation, 
                              devDatabaseName = devDatabaseName,
                              valDatabaseNames = valDatabaseNames,
                              phenotypes = phenotypes, 
                              targets = targets, 
                              phenotype = phenotype, 
                              target = target, 
                              benchmark = F)
        
      }
    }
  }
  
}

externalValidateInPar <- function(resultLocation, devDatabaseName, 
                                  valDatabaseNames, phenotypes, 
                                  targets, phenotype, 
                                  target, benchmark ){
  # internal other Ts and other Os (create settings per models and valdiate in para)
  
  modelLoc <- file.path(resultLocation,devDatabaseName,'models', paste0('T_',target,'_O_',phenotype,ifelse(benchmark,'_ag', '')), 'plpResult')
  
  plpResult <- tryCatch({PatientLevelPrediction::loadPlpResult(modelLoc)}, 
                        error = function(e){ParallelLogger::logInfo(e);return(NULL)})
  

  if(!is.null(plpResult)){
    
    if(!dir.exists(file.path(resultLocation,devDatabaseName,'validation'))){
      dir.create(file.path(resultLocation,devDatabaseName,'validation'), recursive = T)
    }
    
    dataToValidate <- expand.grid(valDatabaseNames,phenotypes,targets)
    colnames(dataToValidate) <- c('databaseName','outcomeId','T')
    dataToValidate$plpDataLoc <- file.path(resultLocation,dataToValidate$databaseName,'data', paste0('T_',dataToValidate$T,ifelse(benchmark,'_ag', '')))
    dataToValidate$modelLoc <- modelLoc
    dataToValidate$saveDirectory <- file.path(resultLocation,devDatabaseName,'validation',dataToValidate$databaseName)
    dataToValidate$analysisId <- paste0('Td_',target,'_Od_',phenotype,ifelse(benchmark,'_ag', ''),'Tv_',dataToValidate$T,'_Ov_',dataToValidate$outcomeId,ifelse(benchmark,'_ag', ''))
    
    ParallelLogger::logInfo(paste0('Validating model ', modelLoc ,' ', nrow(dataToValidate), ' times'))
    
    # create the cluster
    ParallelLogger::logInfo(paste0('Number of cores not specified'))
    cores <- min(nrow(dataToValidate),parallel::detectCores() -3) ## set to max 
    ParallelLogger::logInfo(paste0('Using this many cores ', cores))
    ParallelLogger::logInfo(paste0('Set cores input to use fewer...'))
    
    
    cluster <- ParallelLogger::makeCluster(numberOfThreads = cores)
    ParallelLogger::clusterRequire(cluster, c("PatientLevelPrediction", "Andromeda"))
    
    logger <- ParallelLogger::createLogger(name = "PAR",
                                           threshold = 'INFO', 
                                           appenders = list(ParallelLogger::createFileAppender(layout = ParallelLogger::layoutParallel, 
                                                                                               fileName = file.path(resultLocation,devDatabaseName,'validation','parlog.txt'))))
    ParallelLogger::registerLogger(logger)
    
    # settings: convert dataToValidate to list
    dataToValidateList <- split(dataToValidate, 1:nrow(dataToValidate)) #!!!

    allResults <- ParallelLogger::clusterApply(cluster = cluster, 
                                               x = dataToValidateList, 
                                               fun = valPlpI, 
                                               stopOnError = FALSE,
                                               progressBar = TRUE)
    ParallelLogger::stopCluster(cluster)
    
  }
  
}


valPlpI <- function(settings){
  
  # settings needs: plpDataLoc, modelLoc, outcomeId saveDirectory and analysisId
  
  if(file.exists(file.path(settings$saveDirectory,settings$analysisId, 'valResult.rds'))){
    ParallelLogger::logInfo(paste0('Results exist in ',file.path(settings$saveDirectory,settings$analysisId), ' so skipping'))
    return(NULL)
  }
  
  ParallelLogger::logInfo('Loading data and model')
  plpData <- tryCatch({PatientLevelPrediction::loadPlpData(settings$plpDataLoc)}, 
                      error = function(e){ParallelLogger::logInfo(e); return(NULL)})
  if(is.null(plpData)){return(NULL)}
  plpResult <- tryCatch({PatientLevelPrediction::loadPlpResult(settings$modelLoc)}, 
                        error = function(e){ParallelLogger::logInfo(e); return(NULL)})
  if(is.null(plpResult)){return(NULL)}
  
  ParallelLogger::logInfo("Restricting data to model")
  plpData$covariateData$model <- plpResult$covariateSummary[plpResult$covariateSummary$covariateValue !=0, 'covariateId']
  plpData$covariateData$covariates <- plpData$covariateData$covariates %>% dplyr::inner_join(plpData$covariateData$model)
  plpData$covariateData$model <- NULL
  
  ParallelLogger::logInfo('Running validation')
  popSet <- plpResult$inputSettings$populationSettings
  popSet$plpData = plpData
  popSet$outcomeId <- settings$outcomeId
  population <- do.call(PatientLevelPrediction::createStudyPopulation,popSet)
  
  if(sum(population$outcomeCount>0)<10){
    ParallelLogger::logInfo(paste0('Insufficient outcome count'))
    return(NULL)
  }
  
  result <- tryCatch({PatientLevelPrediction::applyModel(population = population,
                                                         plpData = plpData,
                                                         plpModel = plpResult$model,
                                                         calculatePerformance = T,
                                                         databaseOutput = NULL,
                                                         silent = F)},
                     error = function(e){ParallelLogger::logInfo(e); return(NULL)})
  
  if(is.null(result)){
    return(NULL)
  }
  
  if(!dir.exists(file.path(settings$saveDirectory, settings$analysisId))){
    dir.create(file.path(settings$saveDirectory, settings$analysisId), recursive = T)
  }
  saveRDS(result, file.path(settings$saveDirectory, settings$analysisId, 'valResult.rds'))
  return(result)
}
