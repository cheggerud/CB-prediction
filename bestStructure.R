rm(list=ls())
server = T
if (!server)
  setwd("/Users/pouriaramazi/Desktop/CyB")
source("header.R")
################################################# User's input 
discritizationLevelsRange = 3
## Options for discretization methods: 'hartemink', 'frequency' or 'interval'. 
discretizationMethod = 'interval'
## Partitioning the dataset into train and test:
# For each lake, take the full last year as a test insance--up to 3 years--and only allow up to
# the following number of lakes that only show up in the test dataset:
maxNumberOfLakesThatCanBeExclusiveToTheTestDataset = 5
# What portion of the original dataset does the test dataset takes?
ratioOfTestToOriginalData = .15
# number of digits to be saved in the tables: (e.g., 2 implies 0.98)
digitAccuracy = 2
################################################# User's input for debugging, speed, etc.
# Read the train and test from previous calculations, or do the partitioning to obtain them?
readTrainAndTest = T
# Find the best Bayesian network structure, or find good predictors?
findTheBestStructure = T
findBestPredictors = T
if (findBestPredictors) 
  methods = c('GLM','NB', 'GBM', 'KNN', 'NN', 'MM')
# which variables are not covariates?
nonCovariates = c("lake", "Year")
################################################# Setup 
# Where to save the data?
targetPath = paste("results", "/", sep = "")
dir.create(targetPath, recursive = TRUE, showWarnings = FALSE)
# Where to save the cpt figures of teh best structure?
bestStructureFigurestargetPath = paste(targetPath, "/figures/", sep = "")
dir.create(bestStructureFigurestargetPath, recursive = TRUE, showWarnings = FALSE)
# Emtpy this write file
write("", file = paste(targetPath, "MarkovBlanket.tex", sep = ""))
## AUC scores of differnt number of levels
AUCOverLevelsTrain = matrix(ncol = 3, nrow = length(discritizationLevelsRange))
colnames(AUCOverLevelsTrain) = c("AUROC", "AUPR", "AUPRforARandomClassifier")
rownames(AUCOverLevelsTrain) = discritizationLevelsRange
AUCOverLevelsTest = AUCOverLevelsTrain
################################################# Partition the data
# this should be dataForPrediction.csv.csv
data <- read.csv(file = "data.csv",
                 header = TRUE,
                 sep = ",")
dataBeforeFactorizing = data
################################################# Loop through different discretization levels
for (discritizationLevels in discritizationLevelsRange) {
  ################################################# Divide the data into training (t) and validation (v) parts
  partitionedDadta = trainAndTestPartitionForLakes(
    data,
    ratioOfTestToOriginalData,
    discretizationMethod,
    discritizationLevels,
    readTrainAndTest
  )
  t = partitionedDadta[[1]]
  test = partitionedDadta[[2]]
  str(t)
  str(test)
  ################################################# Find the best-structure 
  if (findTheBestStructure) {
    # bnstruct package requires a special format of the data
    numeric_t = t
    eval(parse(text = paste('numeric_t$',  names(t), ' = as.numeric(numeric_t$', names(t), ")", sep = "")))
    nodeSize = 0
    for (i in 1:dim(t)[2]) {
      nodeSize[i] = nlevels(t[,i])
    }
    datasetForUsingSM <- BNDataset(data = numeric_t,
                                   discreteness = rep('d',dim(t)[2]),
                                   variables = names(t),
                                   node.sizes = nodeSize)
    # Now it can be used
    
    bestStructureOriginal <- learn.network(datasetForUsingSM,
                                   algo = "sm",
                                   scoring.func = "BIC")#,
                                   #mandatory.edges = adj)
    save.image(file=paste(targetPath, 'environment.RData', sep = ""))
    # Convert it to bnlearn format
    temp = empty.graph(names(t))
    amat(temp) = bnstruct::dag(bestStructureOriginal)
    bestStructure = temp
    # Visualize it
    nameOfTheTargetNode = "Target"
    visualizeAndSaveAStructure(bestStructure, 'dot', 
                               filePath = bestStructureFigurestargetPath, 
                               targetNode = nameOfTheTargetNode, 
                               main = paste0("bestStructure", discritizationLevels), 
                               save = TRUE)
    visualizeAndSaveAStructure(bestStructure, 'fdp', 
                               filePath = bestStructureFigurestargetPath, 
                               targetNode = nameOfTheTargetNode, 
                               main = paste0("bestStructure", discritizationLevels), 
                               save = TRUE)
    # What's its AUC score
    net = bn.fit(bestStructure, t, method = 'bayes')
    eval(parse(text = paste0("net", discritizationLevels," = net")))
    net = bn.fit(bestStructure, t, method = 'bayes')
    # make sure Target is numeric 
    test$Target = as.numeric(test$Target)
    if (2 %in% test$Target)
      test$Target = test$Target-1
    # Store the scores
    AUCOverLevelsTrain[as.character(discritizationLevels),] = 
      aucScoreFunction(net, 'Target', t)
    AUCOverLevelsTest[as.character(discritizationLevels),] = 
      aucScoreFunction(net, 'Target', test)
    # Report and save
    temp1 = paste("AUROC and AUPR scores of the best structure with", 
                  discritizationLevels, 
                  "levels on test dataset: ", 
                  paste(round(AUCOverLevelsTest[as.character(discritizationLevels),],3), 
                        collapse = " "))
    temp2 = paste(" on train dataset: ", 
                  paste(round(AUCOverLevelsTrain[as.character(discritizationLevels),],3), 
                        collapse = " "))  
    temp = paste(temp1, temp2, collapse = "\n ")
    write(temp, file = paste(targetPath, "/selectedCovariates.tex", sep = ""), append = TRUE)
    print(temp)
}
    # Save the whole data
    save.image(file=paste(targetPath, 'environment.RData', sep = ""))
}
################################################# Plot the AUC scores
bNNames = c("AUROC", "AUPR", "AUPR of a random classifier")
pdf(paste(targetPath, "AUCScoreOverLevels.pdf", sep = ""))
par(family = "serif")
matplot(AUCOverLevelsTest, type = "n", xlab = "Number of discretization levels" , ylab = "AUC scores", 
        xaxt = "none", cex.lab = 1.3)
axis(side = 1, at = 1:6, labels = 2:7)
grid( )
abline(v = 1:19,  lty = 3, col = "lightgray") 
par(xpd=TRUE, family = "serif")
matplot(AUCOverLevelsTest, type = "b", xaxt = "none", add = T,  lty = 1:3, col = c(1,2,4), pch = 15:17)
legend("topright", legend=bNNames,  lty =1:3, col = c(1,2,4), pch = 15:17, inset=c(0,.3),
       bg= ("white"), horiz=F, pt.cex = 1, cex = 1) 
par(xpd = FALSE) 
dev.off() 
################################################# Choose the best structure
# The one with the highest AUC on test
numberOfLevelsResultingInMaxAUC = 
  names(AUCOverLevelsTest[,1])[match(max(AUCOverLevelsTest[,1]), AUCOverLevelsTest[,1])]
# Print it
cat(numberOfLevelsResultingInMaxAUC, 'levels of quantization resulted in the maximum AUC of ', 
    round(max(AUCOverLevelsTest[,1]),3), '\n')
eval(parse(text = paste0('net = net', numberOfLevelsResultingInMaxAUC)))
################################################# Identify the d-separations of your chosen structure
## List all nodes Y that can be d-separated from 'Target' given some set of nodes Z
singleDseparableNodes(
  net, 'Target', filePath = paste(targetPath, 
    "singleNodeDSeparations.tex", sep = ""))
## Find the maximum number of nodes that can be d-separated from 'Target'
maximumDSeparatedNodesWithMinimumKnownNodes(
  net, 'Target', filePath = paste(targetPath, 
    "maximumNumberOfDSeparatedNodesWithMinimumNumberOfKnownNodes.tex", sep = ""))
################################################# How does the AUC value change as we miss some variables?
# Which nodes we want to delete and then calculate the AUC of the resulting network?
nodesToBeDeleted = mb(net, nameOfTheTargetNode)
## This is the AUROC and AUPRC scores of the best Bayesian Netowork structure 
#as a single node is deleted 
AUCScoreWhenANodeIsDeleted <- matrix(ncol = 3, nrow = length(nodesToBeDeleted))
colnames(AUCScoreWhenANodeIsDeleted) = c("AUROC", "AUPR", "AUPRforARandomClassifier")
rownames(AUCScoreWhenANodeIsDeleted) = nodesToBeDeleted
# save v
originalV = test
for (node in nodesToBeDeleted) {
  v = originalV
  eval(parse(text = paste('v$', node, ' = ""', sep = "")))
  AUCScoreWhenANodeIsDeleted[match(node, nodesToBeDeleted), ] = aucScoreFunction(net, 'Target', v)
}  
# Make sure v has not changed.
v = originalV
# Save the AUC table in latex format
print(xtable(AUCScoreWhenANodeIsDeleted, digits = digitAccuracy), 
      file = paste(targetPath, "AUCScoreWhenANodeIsDeleted.tex", sep =""))
################################################# plot all conditional probabilities for 'Target'
plotAndPrintAllSingleNodeConditionedProbabilities(net, 'Target', t,
                                                  plotFilePath = bestStructureFigurestargetPath,
                                                  tableFilePath = bestStructureFigurestargetPath, 
                                                  plotSave = TRUE, tableSave = TRUE)
################################################# plot enquiry
junction = compile(as.grain(net))
# see the nodes
nodes(net)[1:(length(nodes(net))-1)]
# Which ones should be y (joint or conditioned on)?
yIndices = c(4)
# Which ones should be known variables?
zIndices = c(3)
levels(t[,zIndices[1]])
# What levels?
zlevelIndices = c(3)
plotAndPrintAConditionalProbability(net, 'Target', yIndices, zIndices, zlevelIndices, t,
                                    printTable = TRUE)
################################################# AUC upon missing a group of covariates
# Which nodes we want to delete and then calculate the AUC of the resulting network?
missingGroupOfCovariates_lakeMonitoring_monthIncluded = names(t)[c(1,2,3,4,5,6)]
missingGroupOfCovariates_weather_monthIncluded= names(t)[c(1,7,8,9,10)]
missingGroupOfCovariates_lakeFeatures_monthIncluded = names(t)[c(1,11,12)]
missingGroupOfCovariates_watershed_monthIncluded = names(t)[c(1,13:17)]
#
missingGroupOfCovariates_lakeMonitoring = names(t)[c(2,3,4,5,6)]
missingGroupOfCovariates_weather= names(t)[c(7,8,9,10)]
missingGroupOfCovariates_lakeFeatures = names(t)[c(11,12)]
missingGroupOfCovariates_watershed = names(t)[c(13:17)]
missingGroupOfCovariates_nonDynamicLakeCharacteristics = names(t)[c(6,11,17,5,12)] 
missingGroupOfCovariates_dynamicLakeCharacteristics = names(t)[c(3,4,2)] 
missingGroupOfCovariates_CyBGrowth = names(t)[c(3,4,7,10)] 
missingGroupOfCovariates_agricultureLand = names(t)[c(13,14)]
missingGroupOfCovariates_eutrophication = names(t)[c(3,4)]
missingGroupOfCovariates_localResidentIncluded = names(t)[-c(2,5,7,8,9,10,11,12,13,14,15,16,17,18)]
missingGroupOfCovariates_distantTouristIncluded = names(t)[-c(7,8,9,10,11,12,13,14,15,16,17,18)]
missingGroupOfCovariates_highScores = names(t)[c(2,3,4,5,9,11,13,14,15,16)]
secondaryCovariateIndices = c(5,6,7,8,9,10,11,12,14,15,16,17)
missingGroupOfCovariates_monthAndSecondaryVariables = names(t)[c(1,secondaryCovariateIndices)]
missingGroupOfCovariates_v0_currntAndSecondaryVariables = names(t)[c(2,secondaryCovariateIndices)]
missingGroupOfCovariates_v01_PAndSecondaryVariables = names(t)[c(3,secondaryCovariateIndices)]
missingGroupOfCovariates_v02_NAndSecondaryVariables = names(t)[c(4,secondaryCovariateIndices)]
missingGroupOfCovariates_v11_pastuAndSecondaryVariables = names(t)[c(13,secondaryCovariateIndices)]
#missingGroupOfCovariates_test1 = names(t)[c(14,15,16)]
# Included means AUC upon inclusion (instead of omission)
missingGroupOfCovariates_lakeMonitoringIncluded_monthIncluded = names(t)[-c(1,2,3,4,5,6,18)]
missingGroupOfCovariates_weatherIncluded_monthIncluded = names(t)[-c(1,7,8,9,10,18)]
missingGroupOfCovariates_lakeFeaturesIncluded_monthIncluded = names(t)[-c(1,11,12,18)]
missingGroupOfCovariates_watershedIncluded_monthIncluded = names(t)[-c(1,13:17,18)]
missingGroupOfCovariates_lakeMonitoringIncluded = names(t)[-c(2,3,4,5,6,18)]
missingGroupOfCovariates_weatherIncluded = names(t)[-c(7,8,9,10,18)]
missingGroupOfCovariates_lakeFeaturesIncluded = names(t)[-c(11,12,18)]
missingGroupOfCovariates_watershedIncluded = names(t)[-c(13:17,18)]
missingGroupsfOfCovaraites = 
  c('lakeMonitoring_monthIncluded',
    'weather_monthIncluded',
    'lakeFeatures_monthIncluded',
    'watershed_monthIncluded',
    'lakeMonitoring',
    'weather',
    'lakeFeatures',
    'watershed',
    'nonDynamicLakeCharacteristics',
    'dynamicLakeCharacteristics',
    'CyBGrowth',
    'agricultureLand',
    'eutrophication', 
    'highScores',
    'lakeMonitoringIncluded', 'weatherIncluded', 'lakeFeaturesIncluded', 'watershedIncluded',
    'lakeMonitoringIncluded_monthIncluded', 'weatherIncluded_monthIncluded', 'lakeFeaturesIncluded_monthIncluded', 'watershedIncluded_monthIncluded',
    'localResidentIncluded', 'distantTouristIncluded',
    'monthAndSecondaryVariables', 
    'v0_currntAndSecondaryVariables',
    'v01_PAndSecondaryVariables',
    'v02_NAndSecondaryVariables',
    'v11_pastuAndSecondaryVariables')
## This is the AUROC and AUPRC scores of the best Bayesian Netowork when the values 
#of a group is missing. 
AUCScoreWhenAGroupIsDeleted <- matrix(ncol = 3, nrow = length(missingGroupsfOfCovaraites))
colnames(AUCScoreWhenAGroupIsDeleted) = c("AUC", "AUPR", "AUPRforARandomClassifier")
rownames(AUCScoreWhenAGroupIsDeleted) = missingGroupsfOfCovaraites
## Go through them, each time delete a group and test the performance
for (group in missingGroupsfOfCovaraites) {
  # Construct the test dataset with missing values
  v = test
  eval(parse(text = paste0('nodesToBeDeleted = missingGroupOfCovariates_', group)))
  for (node in nodesToBeDeleted) 
    eval(parse(text = paste0('v$', node, ' = NULL')))
  # Now test it
  AUCScoreWhenAGroupIsDeleted[group,] = aucScoreFunction(net, 'Target', v)
}
# Make sure v has not changed
v = test
# Save the AUC table in latex format
print(xtable(AUCScoreWhenAGroupIsDeleted[,1, drop = F], digits = digitAccuracy), 
      file = paste(targetPath, "AUCScoreWhenAGroupIsDeleted.tex", sep =""))
## Test when only one covariate is present.
AUCScoreWhenOnlyOneNodeIsAvailabele <- matrix(ncol = 3, nrow = dim(t)[2]-1)
colnames(AUCScoreWhenOnlyOneNodeIsAvailabele) = c("AUC", "AUPR", "AUPRforARandomClassifier")
rownames(AUCScoreWhenOnlyOneNodeIsAvailabele) = names(t)[-dim(t)[2]]
# Construct the test dataset with missing values
for (availableNode in names(t)[1:(dim(t)[2]-1)]) {
  v = test
  for (node in exclude(names(t)[1:(dim(t)[2]-1)], availableNode)) 
    eval(parse(text = paste0('v$', node, ' = NULL')))
  # Now test it
  AUCScoreWhenOnlyOneNodeIsAvailabele[availableNode,] = 
      aucScoreFunction(net, 'Target', v)[1]
}
# Make sure v has not changed
v = test
# Save the AUC table in latex format
print(xtable(AUCScoreWhenOnlyOneNodeIsAvailabele[,1, drop = F], digits = digitAccuracy), 
      file = paste(targetPath, "AUCScoreWhenOnlyOneNodeIsAvailabele.tex", sep =""))
################################################# compare with null models
## Store the AUC results in matrix
AUCScoreOfNullModels <- matrix(ncol = 3, nrow = dim(t)[2]-1)
colnames(AUCScoreOfNullModels) = c("AUC", "AUPR", "AUPRforARandomClassifier")
rownames(AUCScoreOfNullModels) = names(t)[-dim(t)[2]]
# Go through all covariates:
for (covariateIndex in 1:(dim(t)[2]-1)) {
  covariate = names(t)[covariateIndex]
  # Construct the model and estimate the parameters
  nullModel = model2network(paste0("[", covariate, "][Target|", covariate,"]"))
  nullNet = bn.fit(nullModel, t[,c(covariateIndex, dim(t)[2]), drop = F], method = "bayes")
  # Construct the test dataset for just one covariate
  testDataForASingleCovariate = test
  for (node in exclude(names(t)[1:(dim(t)[2]-1)], covariate)) 
    eval(parse(text = paste0('testDataForASingleCovariate$', node, ' = NULL')))
  # Now test it
  AUCScoreOfNullModels[covariate,] = 
    aucScoreFunction(nullNet, 'Target', testDataForASingleCovariate)[1]
}
# Save the whole image
save.image(file=paste(targetPath, 'environment.RData', sep = ""))
# Save the AUC table in latex format
print(xtable(AUCScoreOfNullModels[,1, drop = F], digits = digitAccuracy), 
      file = paste(targetPath, "AUCScoreOfNullModels.tex", sep =""))
################################################# Parallel: Find the subset of best covariates
if (server) {
  # The cores:
  cores=detectCores()
  cl <- makeCluster(min(cores[1]-1, 10), outfile = "") #not to overload your computer
  registerDoParallel(cl)
  clusterExport(cl=cl, varlist=c("compile", 'as.grain','querygrain', 'setEvidence',
                                 'roc.curve',
                                 'test', 'net', 'fastAucScoreFunction'))
  ## The matrix to store the AUC of each subset of covariates
  # n: number of test covariates
  n = dim(test)[2]-1
  AUCScoreOfTheCovariateSubset <- matrix(nrow = 2^n)
  colnames(AUCScoreOfTheCovariateSubset) = c("AUC")
  # loop through the covariates
  AUCScoreOfTheCovariateSubset = 
    foreach (subsetSize = 1:n, .combine=rbind, .packages = c('utils')) %dopar% {  
      temp = combn(1:n, subsetSize, fastAucScoreFunction)
      # report
      text = paste0('Done with subsetsize ', subsetSize, '. Max AUC = ', max(temp), 
                    ' for the subsetSize = ',  
                    length(t(combn(1:n, subsetSize))[which(temp %in% max(temp)),]), 
                    ', resulting in the covariates ', 
                    paste(names(test)[t(combn(1:n, subsetSize))[which(temp %in% max(temp)),]], collapse = ","),
                    '\n')
      cat(text)
      cat(text, file = paste0(targetPath, 'whichCovariateSubsetResultsInTheHighestAUC.tex'), 
          append = T)
      # return
      temp
    }
  stopCluster(cl)
} else {
  ## The matrix to store the AUC of each subset of covariates
  # n: number of test covariates
  n = dim(test)[2]-1
  AUCScoreOfTheCovariateSubset <- matrix(nrow = 2^n)
  colnames(AUCScoreOfTheCovariateSubset) = c("AUC")
  # loop through the covariates
  counter = 0
  for (subsetSize in 1:2) {
    temp = combn(1:n, subsetSize, fastAucScoreFunctionOnTrain)
    AUCScoreOfTheCovariateSubset[(counter+1) : (counter + length(combn(1:n, subsetSize)))] =
      temp
    counter = counter + length(combn(1:n, subsetSize)) 
    # report
    cat(max(AUCScoreOfTheCovariateSubset, na.rm = T), '\n')
    text = paste0('Done with subsetsize ', subsetSize, '. Max AUC = ', max(temp), 
                  ' for the subsetSize = ',  
                  length(t(combn(1:n, subsetSize))[which(temp %in% max(temp)),]), 
                  ', resulting in the covariates ', 
                  paste(names(test)[t(combn(1:n, subsetSize))[which(temp %in% max(temp)),]], collapse = ","),
                  '\n')
    
    cat(text)
    cat(text, file = paste0(targetPath, 'whichCovariateSubsetResultsInTheHighestAUC.tex'), 
        append = T)
  }
}
################################################# Find best predictors
if (findBestPredictors) {
  ############################## manual input
  # lake, year, Target
  numberOfNonFeatureVariables = 3
  # Name of target node
  target = tail(names(tOriginal),1)
  ############################## Setup
  # How many features do we have?
  totalNumberOfFeatures = dim(tOriginal)[2] - numberOfNonFeatureVariables
  AUCScores = matrix(nrow = totalNumberOfFeatures, ncol = length(methods),
                     dimnames = list(NULL, methods))
  # Results on the test
  testAUCScores = 
    matrix(nrow = length(methods), ncol = 1, dimnames = list(methods,NULL))
  ############################## data preprocessing
  # Target node must be binary
  if (2 %in% t$Target) {
    tOriginal$Target = as.integer(t$Target) - 1
    test$Target = as.integer(test$Target) - 1
  }
  target = 'Target'
  ############################## Do feature selection
  # library mRMRe
  ## Do the mrmr on the train. So first, Read from the train 
  d = t
  n = dim(d)[2]
  # Turn everything into numeric
  for (i in 1:n) d[, i] = as.numeric(d[, i])
  # Apply MRMR (Pixel and Year are excluded)
  mRMRResults = mRMR.classic(data = mRMR.data(d[,1:n]), target_indices = n, feature_count = n-1)
  featureIndices = rev(unlist(mRMRResults@filters, use.names = F))
  # Print
  names(d[,1:n])[featureIndices]
  # Save the mrmr result
  cat(names(d[,1:n])[featureIndices], file = paste(targetPath, "mRMR_results.tex", sep =""))
  ############################## Go through the methods
  for (method in methods) {
    ########################  obtaining validation
    for (numberOfFeatures in 1:totalNumberOfFeatures) {
      data = tOriginal[,c(1:2,2+featureIndices[(1:numberOfFeatures)],dim(tOriginal)[2])]
      numberOfValidationSamples = round(ratioOfTestToOriginalData*dim(data)[1])
      # v: Test dataset
      indicesOfV = NULL
      numberOfLakesThatAreExclusiveToTheTestDataset = 0
      for (lake in unique(data$lake)) {
        numberOfYearsAssociatedToTheLakeInTheDataSet = length(unique(data[data$lake==lake,2]))
        if (numberOfLakesThatAreExclusiveToTheTestDataset >
            maxNumberOfLakesThatCanBeExclusiveToTheTestDataset) {
          if (numberOfYearsAssociatedToTheLakeInTheDataSet==1) 
            numberOfYearsToIncludeForTheLake = 0
        }
        if (numberOfYearsAssociatedToTheLakeInTheDataSet>2) 
          numberOfYearsToIncludeForTheLake = 2
        if (numberOfYearsAssociatedToTheLakeInTheDataSet==2) 
          numberOfYearsToIncludeForTheLake = 1
        if ( (numberOfYearsAssociatedToTheLakeInTheDataSet==1) & (numberOfYearsToIncludeForTheLake!=0) ){
          numberOfYearsToIncludeForTheLake = 1
          numberOfLakesThatAreExclusiveToTheTestDataset = numberOfLakesThatAreExclusiveToTheTestDataset + 1
        }
        # Now for each year
        for (year in tail(unique(data[data$lake==lake,2]), n = numberOfYearsToIncludeForTheLake)) {
          if (length(indicesOfV) > numberOfValidationSamples)
            break
          # Assumption: the lakes are sorted in ascending order of first, sampled year, and second, sampled month.
          indicesOfV = append(indicesOfV, 
                              intersect(which(data$lake %in% lake), which(data$Year %in% year)))
        }
      }
      v = data[indicesOfV,]
      t = data[setdiff(1:dim(data)[1], indicesOfV),]
      # Delete lake and year
      t[,c(1:2)] <- NULL
      v[,c(1:2)] <- NULL
      # Print tout the dimensions
      cat('The current t and v have ', dim(t)[2], ' and ', dim(v)[2], ' variables respectively.\n')
      ######################## GLM
      if (method == 'GLM') {
        glm.model = glm(as.formula(paste(target, " ~.")), t, family = binomial())
        outputs = predict(glm.model, newdata = v, type = 'response')
        # for MM
        eval(parse(text = paste("tMM", numberOfFeatures, 
                                " <<- data.frame(", target, " = t[,dim(t)[2]])", sep = "")))
        eval(parse(text = paste("vMM", numberOfFeatures, 
                                " <<- data.frame(", target, " = v[,dim(v)[2]])", sep = "")))
        
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) 
          outputsMM <<- predict(glm.model, newdata = t, type = 'response')
      }
      ######################## GBM
      if (method == 'GBM') {
        eval(parse(text = paste('gbm.model <- gbm(', target, '~ ., data = t, 
                            distribution = "bernoulli", n.trees = 100)', sep = ""))) 
        best.iter = gbm.perf(gbm.model, method = "OOB", plot.it = F) #, type = 'response'
        outputs = predict(gbm.model, newdata = v[,1:(dim(v)[2]-1), drop = F], 
                          n.trees = best.iter, type = 'response')
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) 
          outputsMM <<-predict(gbm.model, newdata = t[,1:(dim(v)[2]-1), drop = F], 
                               n.trees = best.iter, type = 'response')
      }
      ################################################# KNN
      if ('KNN' %in% method) {
        if ((numberOfFeatures == 1) | (numberOfFeatures == 2)) {
          outputs = rep(.5, dim(v)[1])
          outputsMM <<-rep(.5, dim(t)[1])
        } else {
          # get around the ``too many ties'' issue
          knn.values = knn(t[,-dim(t)[2], drop = F], v[, -dim(v)[2], drop = F], t[,dim(t)[2]], k = 5, l = 1, prob = TRUE, use.all = F)
          probabilities = attributes(knn.values)$prob
          outputs = rep(NaN, length(probabilities))
          for (i in 1:length(probabilities))
            outputs[i] = if (knn.values[i] == 1) probabilities[i] else (1-probabilities[i]) 
          # Construct the outputs for the training dataset of the 'MM' method
          if ('MM' %in% methods) {
            knn.values = knn(t[,-dim(t)[2], drop = F], t[, -dim(t)[2], drop = F], t[,dim(t)[2]], k = 20, l = 1, prob = TRUE, use.all = F)
            probabilities = attributes(knn.values)$prob
            outputsMM <<-rep(NaN, length(probabilities))
            for (i in 1:length(probabilities))
              outputsMM[i] = if (knn.values[i] == 1) probabilities[i] else (1-probabilities[i])  
          }
        }
      }
      ################################################# NN
      if ('NN' %in% method) {
        # scale the data
        procValues <- preProcess(t, method = c("center", "scale"))
        scaledT <-  predict(procValues, t)
        scaledV <-  predict(procValues, v)
        # train
        nn.model <- nn.train(as.matrix.data.frame(scaledT[,1:(dim(scaledT)[2]-1), drop = F]), 
                             as.matrix.data.frame(scaledT[,dim(scaledT)[2], drop = F]), 
                             hidden = c(round(dim(t)[2]/2)), activationfun = 'sigm')
        # plot(nn.model)
        cat('number of hidden nodes = ', round(dim(t)[2]/2),'\n')
        outputs = nn.predict(nn.model, as.matrix.data.frame(scaledV[,1:(dim(v)[2]-1), drop = F]))
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) 
          outputsMM <<-nn.predict(nn.model, as.matrix.data.frame(scaledT[,1:(dim(t)[2]-1), drop = F]))
      }
      ######################## NB
      if ('NB' %in% method) {
        tDiscrete = t
        vDiscrete = v
        # Factorize them 
        for (i in 1:dim(t)[2]) {
          tDiscrete[,i] = factor(tDiscrete[,i])
          vDiscrete[,i] = factor(vDiscrete[,i], levels = levels(tDiscrete[,i]))
        }
        # Train the Naive Bayesian
        eval(parse(text = paste("nbn.model <- naive.bayes(tDiscrete,'", target,"')", sep = "")))
        ## Outputs
        nbn.model = bn.fit(nbn.model, tDiscrete, method = "bayes")
        junction = compile(as.grain(nbn.model))
        eval(parse(text = paste('targets = v$', target, sep = "")))
        outputs = rep(0, length(targets))
        var = head(names(v), -1)
        values = matrix(nrow = 1, ncol = dim(v)[2]-1)
        for (i in 1:length(targets)){
          for (j in 1:(dim(v)[2]-1))
            values[j] = as.character(vDiscrete[i, j])
          outputs[i] = querygrain(setEvidence(junction, 
                                              nodes = var, 
                                              states = values), 
                                  nodes = target, 
                                  type = "joint")[2]
        }
        #outputs = predict(nbn.model, vDiscrete)
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) {
          outputsMM <<-rep(0, dim(t)[1])
          var = head(names(t), -1)
          values = matrix(nrow = 1, ncol = dim(t)[2]-1)
          for (i in 1:dim(t)[1]){
            for (j in 1:(dim(t)[2]-1))
              values[j] = as.character(tDiscrete[i, j])
            outputsMM[i] = querygrain(setEvidence(junction, 
                                                  nodes = var, 
                                                  states = values), 
                                      nodes = target, 
                                      type = "joint")[2]
          }
        }
      }
      ################################################# MM
      if ('MM' %in% method) {
        # data
        eval(parse(text = paste0("tMM = tMM", numberOfFeatures)))
        eval(parse(text = paste0("vMM = vMM", numberOfFeatures)))
        glm.model = glm(as.formula(paste(target, " ~.")), tMM, family = binomial())
        outputs = predict(glm.model, newdata = vMM, type = 'response')
        # just to have the returning `list` complete
        outputsMM = outputs
      }
      # Construct the dataset for MM
      if (('MM' %in% methods) && (method != 'MM')) {
        eval(parse(text = paste("tMM", numberOfFeatures,
                                "$", method, " <- outputsMM", sep = "")))
        eval(parse(text = paste("vMM", numberOfFeatures, 
                                "$", method, " <- outputs", sep = "")))
      }
      ################################################# AUC
      eval(parse(text = paste('targets = v$', target, sep = "")))
      ROC = roc.curve(scores.class0 = outputs[(targets == 1)], 
                      scores.class1 = outputs[(targets == 0)])
      # Store the AUC scores for the fold
      AUCScores[numberOfFeatures, method] = ROC$auc[[1]]
      # print
      cat('AUC score of ', method, ': ', ROC$auc[[1]], '.\n')
    }
  }
  ################################################# test the model(s) with the best AUROC scores
  ## Identify the highest AUC earning  history and feature numbers of each method
  # index of the highest scores of each method
  featureNumberOfHighestEarningAUCs = 
    matrix(nrow = length(methods), ncol = 1, dimnames = list(methods,c('numberOfFeatures')))
  for (method in methods)
    featureNumberOfHighestEarningAUCs[method,] = 
      which(AUCScores[, method] == max(AUCScores[, method], na.rm = T), arr.ind = F)[1]
  ## now train and test
  for (method in methods) {
    ########################  obtaining validation
      data = dataBeforeFactorizing
      # Target node must be binary
      if (2 %in% data$Target) 
        data$Target = as.integer(data$Target) - 1
      data = data[,c(1:2,2+featureIndices[(1:featureNumberOfHighestEarningAUCs[method,])],dim(data)[2])]
      numberOfValidationSamples = round(ratioOfTestToOriginalData*dim(data)[1])
      # v: Test dataset
      indicesOfV = NULL
      numberOfLakesThatAreExclusiveToTheTestDataset = 0
      for (lake in unique(data$lake)) {
        numberOfYearsAssociatedToTheLakeInTheDataSet = length(unique(data[data$lake==lake,2]))
        if (numberOfLakesThatAreExclusiveToTheTestDataset >
            maxNumberOfLakesThatCanBeExclusiveToTheTestDataset) {
          if (numberOfYearsAssociatedToTheLakeInTheDataSet==1) 
            numberOfYearsToIncludeForTheLake = 0
        }
        if (numberOfYearsAssociatedToTheLakeInTheDataSet>2) 
          numberOfYearsToIncludeForTheLake = 2
        if (numberOfYearsAssociatedToTheLakeInTheDataSet==2) 
          numberOfYearsToIncludeForTheLake = 1
        if ( (numberOfYearsAssociatedToTheLakeInTheDataSet==1) & (numberOfYearsToIncludeForTheLake!=0) ){
          numberOfYearsToIncludeForTheLake = 1
          numberOfLakesThatAreExclusiveToTheTestDataset = numberOfLakesThatAreExclusiveToTheTestDataset + 1
        }
        # Now for each year
        for (year in tail(unique(data[data$lake==lake,2]), n = numberOfYearsToIncludeForTheLake)) {
          if (length(indicesOfV) > numberOfValidationSamples)
            break
          # Assumption: the lakes are sorted in ascending order of first, sampled year, and second, sampled month.
          indicesOfV = append(indicesOfV, 
                              intersect(which(data$lake %in% lake), which(data$Year %in% year)))
        }
      }
      v = data[indicesOfV,]
      t = data[setdiff(1:dim(data)[1], indicesOfV),]
      # Delete lake and year
      t[,c(1:2)] <- NULL
      v[,c(1:2)] <- NULL
      # Print tout the dimensions
      cat('The current t and v have ', dim(t)[2], ' and ', dim(v)[2], ' variables respectively.\n')
      ######################## GLM
      if (method == 'GLM') {
        glm.model = glm(as.formula(paste(target, " ~.")), t, family = binomial())
        outputs = predict(glm.model, newdata = v, type = 'response')
        outputsGLM = outputs
        # for MM
        eval(parse(text = paste("tMM", numberOfFeatures, 
                                " <<- data.frame(", target, " = t[,dim(t)[2]])", sep = "")))
        eval(parse(text = paste("vMM", numberOfFeatures, 
                                " <<- data.frame(", target, " = v[,dim(v)[2]])", sep = "")))
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) 
          outputsMM <<- predict(glm.model, newdata = t, type = 'response')
      }
      ######################## GBM
      if (method == 'GBM') {
        eval(parse(text = paste('gbm.model <- gbm(', target, '~ ., data = t, 
                                distribution = "bernoulli", n.trees = 10000)', sep = ""))) 
        best.iter = gbm.perf(gbm.model, method = "OOB", plot.it = F) #, type = 'response'
        outputs = predict(gbm.model, newdata = v[,1:(dim(v)[2]-1), drop = F], 
                          n.trees = best.iter, type = 'response')
        outputsGBM = outputs
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) 
          outputsMM <<-predict(gbm.model, newdata = t[,1:(dim(v)[2]-1), drop = F], 
                               n.trees = best.iter, type = 'response')
      }
      ################################################# KNN
      if ('KNN' %in% method) {
        if ((numberOfFeatures == 1) | (numberOfFeatures == 2)) {
          outputs = rep(.5, dim(v)[1])
          outputsMM <<-rep(.5, dim(t)[1])
        } else {
          # get around the ``too many ties'' issue
          knn.values = knn(t[,-dim(t)[2], drop = F], v[, -dim(v)[2], drop = F], t[,dim(t)[2]], k = 5, l = 1, prob = TRUE, use.all = F)
          probabilities = attributes(knn.values)$prob
          outputs = rep(NaN, length(probabilities))
          for (i in 1:length(probabilities))
            outputs[i] = if (knn.values[i] == 1) probabilities[i] else (1-probabilities[i]) 
          outputsKNN = outputs
          # Construct the outputs for the training dataset of the 'MM' method
          if ('MM' %in% methods) {
            knn.values = knn(t[,-dim(t)[2], drop = F], t[, -dim(t)[2], drop = F], t[,dim(t)[2]], k = 20, l = 1, prob = TRUE, use.all = F)
            probabilities = attributes(knn.values)$prob
            outputsMM <<-rep(NaN, length(probabilities))
            for (i in 1:length(probabilities))
              outputsMM[i] = if (knn.values[i] == 1) probabilities[i] else (1-probabilities[i])  
          }
        }
      }
      ################################################# NN
      if ('NN' %in% method) {
        # scale the data
        procValues <- preProcess(t, method = c("center", "scale"))
        scaledT <-  predict(procValues, t)
        scaledV <-  predict(procValues, v)
        # train
        nn.model <- nn.train(as.matrix.data.frame(scaledT[,1:(dim(scaledT)[2]-1), drop = F]), 
                             as.matrix.data.frame(scaledT[,dim(scaledT)[2], drop = F]), 
                             hidden = c(round(dim(t)[2]/2)), activationfun = 'sigm')
        # plot(nn.model)
        cat('number of hidden nodes = ', round(dim(t)[2]/2),'\n')
        outputs = nn.predict(nn.model, as.matrix.data.frame(scaledV[,1:(dim(v)[2]-1), drop = F]))
        outputsNN = as.numeric(outputs)
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) 
          outputsMM <<-nn.predict(nn.model, as.matrix.data.frame(scaledT[,1:(dim(t)[2]-1), drop = F]))
      }
      ######################## NB
      if ('NB' %in% method) {
        tDiscrete = t
        vDiscrete = v
        # Factorize them 
        for (i in 1:dim(t)[2]) {
          tDiscrete[,i] = factor(tDiscrete[,i])
          vDiscrete[,i] = factor(vDiscrete[,i], levels = levels(tDiscrete[,i]))
        }
        # Train the Naive Bayesian
        eval(parse(text = paste("nbn.model <- naive.bayes(tDiscrete,'", target,"')", sep = "")))
        ## Outputs
        nbn.model = bn.fit(nbn.model, tDiscrete, method = "bayes")
        junction = compile(as.grain(nbn.model))
        eval(parse(text = paste('targets = v$', target, sep = "")))
        outputs = rep(0, length(targets))
        var = head(names(v), -1)
        values = matrix(nrow = 1, ncol = dim(v)[2]-1)
        for (i in 1:length(targets)){
          for (j in 1:(dim(v)[2]-1))
            values[j] = as.character(vDiscrete[i, j])
          outputs[i] = querygrain(setEvidence(junction, 
                                              nodes = var, 
                                              states = values), 
                                  nodes = target, 
                                  type = "joint")[2]
        }
        outputsNB = outputs
        #outputs = predict(nbn.model, vDiscrete)
        # Construct the outputs for the training dataset of the 'MM' method
        if ('MM' %in% methods) {
          outputsMM <<-rep(0, dim(t)[1])
          var = head(names(t), -1)
          values = matrix(nrow = 1, ncol = dim(t)[2]-1)
          for (i in 1:dim(t)[1]){
            for (j in 1:(dim(t)[2]-1))
              values[j] = as.character(tDiscrete[i, j])
            outputsMM[i] = querygrain(setEvidence(junction, 
                                                  nodes = var, 
                                                  states = values), 
                                      nodes = target, 
                                      type = "joint")[2]
          }
        }
      }
      ################################################# MM
      if ('MM' %in% method) {
        # data
        eval(parse(text = paste0("tMM = tMM", numberOfFeatures)))
        eval(parse(text = paste0("vMM = vMM", numberOfFeatures)))
        glm.model = glm(as.formula(paste(target, " ~.")), tMM, family = binomial())
        outputs = predict(glm.model, newdata = vMM, type = 'response')
        # just to have the returning `list` complete
        outputsMM = outputs
      }
      # Construct the dataset for MM
      if (('MM' %in% methods) && (method != 'MM')) {
        eval(parse(text = paste("tMM", numberOfFeatures,
                                "$", method, " <- outputsMM", sep = "")))
        eval(parse(text = paste("vMM", numberOfFeatures, 
                                "$", method, " <- outputs", sep = "")))
      }
      ################################################# AUC
      eval(parse(text = paste('targets = v$', target, sep = "")))
      ROC = roc.curve(scores.class0 = outputs[(targets == 1)], 
                      scores.class1 = outputs[(targets == 0)])
      # Store the AUC scores for the fold
      testAUCScores[method,1] = ROC$auc[[1]]
      # print
      cat('AUC score of ', method, ': ', ROC$auc[[1]], '.\n')
  }
  # Save the whole data
  save.image(file=paste(targetPath, 'prediction_environment.RData', sep = ""))
}



