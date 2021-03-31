
developInPar <- function(targets,outcomes, ageGenderOnly = F, outputFolder, databaseName, seed ){
  model <- PatientLevelPrediction::setLassoLogisticRegression()
  
  for(targetId in targets ){
    for(outcomeSets in outcomes){
      
      ParallelLogger::logInfo(paste0('Running models for: target ', targetId, 
                                     ' and outcomes ', paste0(outcomeSets, collapse = '', sep=',')))
      
      # create the cluster
      ParallelLogger::logInfo(paste0('Number of cores not specified'))
      cores <- length(outcomeSets) 
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
                       minCovariateFraction = 0.0001,
                       saveDirectory= file.path(outputFolder,databaseName,'models'), 
                       savePlpData=F, 
                       savePlpResult=T, 
                       savePlpPlots = F, 
                       saveEvaluation = F, splitSeed = seed,
                       analysisId = paste0('T_',targetId,'_O_',outcomeId,ifelse(ageGenderOnly,'_ag','')))
        return(result)
      }
      runSettings <- lapply(outcomeSets, getRunSettings)
      
      allResults <- ParallelLogger::clusterApply(cluster = cluster, 
                                                 x = runSettings, 
                                                 fun = runPlpI, 
                                                 stopOnError = FALSE,
                                                 progressBar = TRUE)
      ParallelLogger::stopCluster(cluster)
      
    }
    
  }
}

runPlpI <- function(settings){
  
  if(dir.exists(file.path(settings$saveDirectory,settings$analysisId, 'plpResult'))){
    ParallelLogger::logInfo(paste0('Results exist in ',file.path(settings$saveDirectory,settings$analysisId), ' so skipping'))
    return(NULL)
  }
  
  settings$plpData <- PatientLevelPrediction::loadPlpData(settings$plpDataLoc)
  settings$plpDataLoc <- NULL
  
  settings$population <- PatientLevelPrediction::createStudyPopulation(settings$plpData,
                                               outcomeId = settings$outcomeId,
                                               firstExposureOnly = FALSE,
                                               washoutPeriod = 0,
                                               removeSubjectsWithPriorOutcome = FALSE,
                                               priorOutcomeLookback = 99999,
                                               requireTimeAtRisk = T,
                                               minTimeAtRisk=1,
                                               riskWindowStart = 0,
                                               startAnchor = 'cohort start',
                                               riskWindowEnd = 365,
                                               endAnchor = 'cohort start')
  
  if(sum(settings$population$outcomeCount>0)<25){
    ParallelLogger::logInfo(paste0('Insufficient outcome count'))
    return(NULL)
  }
  
  settings$outcomeId<- NULL
  
  getMinCov <- function(x){
    if(x>1000){
      return(0.0001)
    }
    if(x>500){
      return(0.001)
    }
    if(x>200){
      return(0.01)
    }
    return(0.05)
  }
  
  settings$minCovariateFraction <- getMinCov(sum(settings$population$outcomeCount>0))
  
  result <- tryCatch({do.call(PatientLevelPrediction::runPlp, settings)},
                     error = function(e){ParallelLogger::logInfo(e); return(NULL)})
  return(result)
}



internalValidateWrapper <- function(resultLocation, devDatabaseName,phenotypeGroup,targets){
  
  
  for(phenotypes in phenotypeGroups){
    for(phenotype in phenotypes){
      for(target in targets){
        
        internalValidateInPar(resultLocation = resultLocation, 
                              devDatabaseName = devDatabaseName, 
                              phenotypes = phenotypes, 
                              targets = targets, 
                              phenotype = phenotype, 
                              target = target, 
                              benchmark = F)
        
        internalValidateInPar(resultLocation = resultLocation, 
                              devDatabaseName = devDatabaseName, 
                              phenotypes = phenotypes, 
                              targets = targets, 
                              phenotype = phenotype, 
                              target = target,
                              benchmark = T)
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
    ind <- dataToValidate$outcomeId == phenotype && dataToValidate$T == target
    dataToValidate <- dataToValidate[!ind,]
    
    dataToValidate$plpDataLoc <- file.path(resultLocation,dataToValidate$databaseName,'data', paste0('T_',dataToValidate$T,ifelse(benchmark,'_ag', '')))
    dataToValidate$modelLoc <- modelLoc
    dataToValidate$saveDirectory <- file.path(resultLocation,devDatabaseName,'validation',dataToValidate$databaseName)
    dataToValidate$analysisId <- paste0('Td_',target,'_Od_',phenotype,ifelse(benchmark,'_ag', ''),'_Tv_',dataToValidate$T,'_Ov_',dataToValidate$outcomeId,ifelse(benchmark,'_ag', ''))
    
    ParallelLogger::logInfo(paste0('Internally validating model across phenotypes ', modelLoc ,' ', nrow(dataToValidate), ' times'))
    
    # create the cluster
    ParallelLogger::logInfo(paste0('Number of cores not specified'))
    cores <- parallel::detectCores() - 1## set to max 
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
  
  
  for(phenotypes in phenotypeGroups){
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
        
        ParallelLogger::logInfo('Eternally validating benchmark model')
        externalValidateInPar(resultLocation = resultLocation, 
                              devDatabaseName = devDatabaseName,
                              valDatabaseNames = valDatabaseNames,
                              phenotypes = phenotypes, 
                              targets = targets, 
                              phenotype = phenotype, 
                              target = target,
                              benchmark = T)
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
    cores <- parallel::detectCores() -1 ## set to max 
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
  
  plpData <- PatientLevelPrediction::loadPlpData(settings$plpDataLoc)
  plpResult <- PatientLevelPrediction::loadPlpData(settings$modelLoc)
  
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
  
  if(!dir.exists(file.path(settings$saveDirectory, settings$analysisId))){
    dir.create(file.path(settings$saveDirectory, settings$analysisId), recursive = T)
  }
  saveRDS(result, file.path(settings$saveDirectory, settings$analysisId, 'valResult.rds'))
  return(result)
}
