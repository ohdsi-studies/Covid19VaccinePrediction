
developInPar <- function(targets,outcomes,outputFolder, seed ){
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
        result <- list(plpDataLoc = file.path(outputFolder,'data',paste0('T_', targetId)),
                       outcomeId = outcomeId,
                       modelSettings = model,
                       saveDirectory= file.path(outputFolder,'models'), 
                       savePlpData=F, 
                       savePlpResult=T, 
                       savePlpPlots = F, 
                       saveEvaluation = F, splitSeed = seed,
                       analysisId = paste0('T_',targetId,'_O_',outcomeId))
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
  settings$outcomeId<- NULL
  
  result <- do.call(PatientLevelPrediction::runPlp, settings)
  return(result)
}
