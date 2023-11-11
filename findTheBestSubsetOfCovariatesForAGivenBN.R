rm(list=ls())
server = T
source("header.R")
targetPath = paste("results", "/", sep = "")
load(file=paste(targetPath, 'environment.RData', sep = ""))
################################################# Choose the best structure
# The one with the highest AUC on test
numberOfLevelsResultingInMaxAUC = 
  names(AUCOverLevelsTest[,1])[match(max(AUCOverLevelsTest[,1]), AUCOverLevelsTest[,1])]
# Print it
cat(numberOfLevelsResultingInMaxAUC, 'levels of quantization resulted in the maximum AUC of ', 
    round(max(AUCOverLevelsTest[,1]),3), '\n')
eval(parse(text = paste0('net = net', numberOfLevelsResultingInMaxAUC)))
################################################# Parallel: Find the subset of best covariates
if (server) {
  # The cores:
  cores=detectCores()
  cl <- makeCluster(min(cores[1]-1, 10), outfile = "") #not to overload your computer
  registerDoParallel(cl)
    clusterExport(cl=cl, varlist=c("compile", 'as.grain','querygrain', 'setEvidence',
                                   'roc.curve',
                                   'test', 't', 'net', 'fastAucScoreFunctionOnTest',
                                   'fastAucScoreFunctionOnTest'))
  ## The matrix to store the AUC of each subset of covariates
  # n: number of test covariates
  n = dim(test)[2]-1
  AUCScoreOfTheCovariateSubset <- matrix(nrow = 2^n-1)
  colnames(AUCScoreOfTheCovariateSubset) = c("AUC")
  # On test or train?
  datasetToCheck = 'test'
  # loop through the covariates
  AUCScoreOfTheCovariateSubset = 
    foreach (subsetSize = 1:n, .combine=rbind, .packages = c('utils')) %dopar% {  
      # Test or train?
      temp = t(combn(1:n, subsetSize, fastAucScoreFunctionOnTest))
      #temp = t(combn(1:n, subsetSize, fastAucScoreFunctionOnTrain))
      # report
      text = paste0(datasetToCheck, 
                    ': Done with subsetsize ', subsetSize, '. Max AUC = ', max(temp), 
                    ', resulting in the covariates ', 
                    paste(names(test)[t(combn(1:n, subsetSize))[which(temp %in% max(temp)),]], collapse = ","),
                    '\n')
      cat(text)
      cat(text, file = paste0(targetPath, 'whichCovariateSubsetResultsInTheHighestAUCOn',
                              datasetToCheck, '.tex'),
          append = T)
      # return it
      temp
  }
  stopCluster(cl)
} 
################################################# 