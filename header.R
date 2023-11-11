# The RBGL package needs to be installed separately:
if (!('RBGL' %in% installed.packages()[,"Package"])) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("RBGL")
}
if (!('Rgraphviz' %in% installed.packages()[,"Package"])) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("Rgraphviz", version = "3.8")
}
# Then the rest of the pacakges
list.of.packages = c('gRbase','bnlearn','bnstruct','pROC','PRROC','psych','gRain','MASS', 
                     'xtable','lubridate','crayon', 'Rgraphviz', 'bnclassify', 'gtools',
                     'foreach', 'doParallel', 'mRMRe', 'gbm', 'class', 'caret', 'deepnet')
# Install missing packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
# Call for library
for (package in list.of.packages)
  eval(parse(text = paste('library(', package, ')', sep = "")))
################################################# function: check d-separation
dseparation <- function(z) {
  if (!(dsep(bestStructure, x = 'Infestation', y = y, z = z))) {
    cat("'Infestation' is d-separated from '", y, "'given", "\n")
  }
}
################################################# function: exclude
# exclude elements from a vector
exclude <- function(mainVector, elementsToBeExcluded) {
  for (element in elementsToBeExcluded) 
    mainVector = mainVector[!(mainVector %in% element)]
  return(mainVector)
}
################################################# function: what portion of the loop and how much time is left
howMuchTimeLeft <- function(iteration, sizeOfLoop, beginningTime, numberOfTimesVisitedHere, 
                            thePortionAfterWhichTimeLeftIsPrinted, printDay = FALSE) {
  # This function prints how much time is left until a loop finisihes.  
  # BeginningTime = Sys.time() and
  #  numberOfTimesVisitedHere <<- 0
  #  must be placed right before the loop, and the function itself inside the loop
  iterationsAfterWhichTimeLeftIsPrinted = round(sizeOfLoop*thePortionAfterWhichTimeLeftIsPrinted)
  if (iteration %% iterationsAfterWhichTimeLeftIsPrinted == 0){
    if (numberOfTimesVisitedHere == 0) {
      secondTime <<- Sys.time()
    }
    numberOfTimesVisitedHere <<- numberOfTimesVisitedHere + 1
    temp = seconds_to_period(round((1/thePortionAfterWhichTimeLeftIsPrinted - 
                                      numberOfTimesVisitedHere)*
                                     difftime(secondTime,beginningTime, units = "secs")))
    if (printDay == TRUE) {
      cat(sprintf('%00d%% complete. %02d days %02d:%02d:%02d left. ', 
                  round(100*iteration/sizeOfLoop,0),
                  day(temp), temp@hour, minute(temp), second(temp)))
    } else {
      cat(sprintf('%00d%% complete. %02d:%02d:%02d left. ', 
                  round(100*iteration/sizeOfLoop,0),
                  temp@hour, minute(temp), second(temp)))      
    }
  }
}  
################################################# function: checking for cycles 
doesMatrixContainCycles <- function(a){
  n = dim(a)[1]
  An = a 
  for (ii in 2:n) {
    An = An%*%a; # do not re-compute A^n from skratch
    if (tr(An) != 0) {
      #cat('got cycles\n')
      return(TRUE)
    }
  }
  return(FALSE)
}
################################################# function: checking undirected links
doesMatrixContainUndirectedLink <- function(a){
  for (i in dim(a)[1]) {
    for (j in dim(a)[2]){
      if ( (a[i,j] == 1) & (a[j,i] == 1)){
        return(TRUE)
      }
    }
  }
  return(FALSE)
}  
################################################# function: Visualize and save a bnlearn structure 
visualizeAndSaveAStructure <- function(bayesianStructure, layout, filePath, targetNode,
                                       main = "", shape = "ellipse", colorMarkovBlanket = TRUE,
                                       save = FALSE, visualize = TRUE) {
  # Adjust the title if it's hierarchical
  if (layout == 'dot')
    main = paste(main, '_hierarchical', sep = "")
  ## Visualize?  
  if (visualize & !colorMarkovBlanket)
    graphviz.plot(bayesianStructure, 
                  layout = layout, 
                  shape = shape, 
                  highlight = list(nodes = targetNode), 
                  main = main)
  if (visualize & colorMarkovBlanket) {
    if (length(mb(bayesianStructure, targetNode)) > 0) 
      graph = graphviz.plot(bayesianStructure, 
                      layout = layout, 
                      shape = shape, 
                      highlight = list(nodes = mb(bayesianStructure, targetNode), fill = "cyan"), 
                      main = main)
    else {
      graph = graphviz.plot(bayesianStructure, 
                            layout = layout, 
                            shape = shape, 
                            main = main)
    }
    nodeRenderInfo(graph)$col[match(mb(bayesianStructure, targetNode), nodes(bayesianStructure))] = "darkblue" 
    nodeRenderInfo(graph)$col[match(targetNode, nodes(bayesianStructure))] = "darkred" 
    nodeRenderInfo(graph)$fill[match(targetNode, nodes(bayesianStructure))] = "orangered" 
    renderGraph(graph)
  }
  ## Save?
  # Make sure the directory exists
  if (save) 
    dir.create(filePath, recursive = TRUE, showWarnings = FALSE)
  #
  if (save & !colorMarkovBlanket) {
    pdf(paste(filePath, main, ".pdf", sep = ""))
    graphviz.plot(bayesianStructure, 
                  layout = layout, 
                  shape = shape, 
                  highlight = list(nodes = targetNode), 
                  main = network)
    dev.off()
  }
  if (save & colorMarkovBlanket) {
    pdf(paste(filePath, main, ".pdf", sep = ""))
    if (length(mb(bayesianStructure, targetNode)) > 0) 
      graph = graphviz.plot(bayesianStructure, 
                            layout = layout, 
                            shape = shape, 
                            highlight = list(nodes = mb(bayesianStructure, targetNode), fill = "cyan"), 
                            main = main)
    else {
      graph = graphviz.plot(bayesianStructure, 
                            layout = layout, 
                            shape = shape, 
                            main = main)
    }
    nodeRenderInfo(graph)$col[match(mb(bayesianStructure, targetNode), nodes(bayesianStructure))] = "darkblue" 
    nodeRenderInfo(graph)$col[match(targetNode, nodes(bayesianStructure))] = "darkred" 
    nodeRenderInfo(graph)$fill[match(targetNode, nodes(bayesianStructure))] = "orangered" 
    renderGraph(graph)
    dev.off()
  }
}
################################################# function: maximum numnber of nodes that can be d-saparated
## Find the maximum number of nodes y that can be d-separated from a target node, with the minimum number of 
# known nodes z
maximumDSeparatedNodesWithMinimumKnownNodes <- function(bayesianNetwork, targetNode, filePath = "") {
  ## bayesianNetwork: an object of class bn  
  ## targetNode: string  
  # empty the file
  write("",filePath)
  numberOfNodes = length(nodes(bayesianNetwork))  
  # Loop from maximum number of nodes to minimum
  foundADSeparation = FALSE
  # -2: one to exclude `infestation', one to exclude z
  for (i in (numberOfNodes-2):1) {
    if (foundADSeparation == TRUE) break
    allSubsetsOfSizeIOfTheNodes = combn(exclude(nodes(bayesianNetwork), targetNode), i)
    # Loop through all possible i-size subsets
    for (j in 1:dim(allSubsetsOfSizeIOfTheNodes)[2]) {
      # Check if all nodes in the j-th column of `allSubsetsOfSizeIOfTheNodes' are d-separated
      # Y = allSubsetsOfSizeIOfTheNodes[, j]
      # check all possible z among all size-k subsets of the remaining nodes
      for (k in 1:(numberOfNodes-i-1)) {
        allSubsetsOfSizeKOFTheRemainingNondes = 
          combn(exclude(nodes(bayesianNetwork), c(targetNode, allSubsetsOfSizeIOfTheNodes[, j])), k)
        for (l in 1:dim(allSubsetsOfSizeKOFTheRemainingNondes)[2]) {
          z = allSubsetsOfSizeKOFTheRemainingNondes[, l]
          # Now loop over all nodes y of Y
          allNodesYAreDSeparated = TRUE 
          for (y in allSubsetsOfSizeIOfTheNodes[, j]) {
            if (!(dsep(bayesianNetwork, x = targetNode, y = y, z = z))) {
              allNodesYAreDSeparated = FALSE 
              break
            }
          }
          if (allNodesYAreDSeparated) {
            foundADSeparation = TRUE
            if (filePath == "") {
              cat(targetNode, " _||_ ", allSubsetsOfSizeIOfTheNodes[, j], ' | ', z, '\n') 
            } else {
              write(paste(targetNode, 
                          " _||_ ",
                          paste(allSubsetsOfSizeIOfTheNodes[, j], collapse = " "),
                          ' | ', 
                          paste(z, collapse = " ")),
                    file = filePath, 
                    append = TRUE)
            }
          }
        }
      }
    }
  }    
}
################################################# function: find all d-separable single nodes 
## List all nodes Y that can be d-separated from 'Infestation' given some set of nodes Z
singleDseparableNodes <- function(bayesianNetwork, targetNode, filePath) {
  ## bayesianNetwork: an object of class bn  
  ## targetNode: string  
  # Empty the write file
  write("", file = filePath)
  numberOfNodes = length(nodes(bayesianNetwork))  
  for (y in exclude(nodes(bayesianNetwork), targetNode)) {
    # When y gets separated by z, z will be added to blacklist. Checkin d-separation with more elements 
    # included in blackList is unnecessary.
    counter = 1
    blackList = list()
    ADseparationFoundInPreviousNumberOfVariables = FALSE
    # a loop over the number of variables needed to d-separate y from Infestation. 
    for (i in 1:(numberOfNodes)) {
      if ((numberOfNodes - length(blackList) - 2 < i) | ADseparationFoundInPreviousNumberOfVariables) break
      temp = combn(exclude(nodes(bayesianNetwork), c(targetNode, y)), i)
      # Check if it temp includes any of the blacklists
      areBlackListsExcludedFromZ = TRUE
      for (componentOfBlackList in blackList) 
        areBlackListsExcludedFromZ = areBlackListsExcludedFromZ & (all(componentOfBlackList %in% temp))
      if (areBlackListsExcludedFromZ) {
        for (j in 1:dim(temp)[2]) {
          z = temp[, j]
          if (dsep(bestStructure, x = targetNode, y = y, z = z)) {
            # Write them
            write(paste("targetNode _||_ ", y, ' | ', paste(z, collapse = " ")),
                  file = paste(targetPath, "singleNodeDSeparations.tex", sep = ""), 
                  append = TRUE)
            blackList[[counter]] = z
            counter = counter + 1
            ADseparationFoundInPreviousNumberOfVariables = TRUE
          }
        }
      }  
    }
  }
}
################################################# function: AUC
aucScoreFunction <- function(net, targetNode = 'Target', v, networkName = '', 
                             exactInference = TRUE, returnRandom = TRUE) {
  ## v: validation data set  
  ## net: a class bn (from bnlearn) object  
  if ((networkName =='bestStructure') | (networkName == 'NaiveBayes'))
    junction = compile(as_grain(net))
  else {
    junction = compile(as.grain(net)) 
    }
  # Determine the real(target) and learned (outputs) values
  eval(parse(text = paste('targets = v$', targetNode, sep = "")))
  outputs = rep(0, length(targets))
  var = head(names(v), -1)
  values = matrix(nrow = 1, ncol = dim(v)[2]-1)
  ## for time
  time1 = Sys.time()
  numberOfTimesVisitedHere <<- 0
  ##
  for (i in 1:length(targets)){
    #howMuchTimeLeft(i, length(targets), time1, numberOfTimesVisitedHere, 1/20)
    if (exactInference) {
      for (j in 1:(dim(v)[2]-1))
        values[j] = as.character(v[i, j])
      outputs[i] = querygrain(setEvidence(junction, 
                                          nodes = var, 
                                          states = values), 
                              nodes = targetNode, 
                              type = "joint")[2]
    } else {
      str = paste("(",  var, " == '",
                  sapply(v[i, 1:dim(v)[2]-1], as.character), "')",
                  sep = "", collapse = " & ")
      str2 = paste("(",  targetNode, " == '1')", sep = "")
      cmd = paste("cpquery(net, ", str2, ", ", str, ")", sep = "")
      outputs[i] = eval(parse(text = cmd))
    }
  }
  if (length(unique(outputs)) == 1)
    if (returnRandom) {
      return(c(.5,NaN,NaN))
      cat('\n')
    } else {
      return(c(.5,NaN))
      cat('\n')
    }  
  #roc(targets, outputs, plot = TRUE)
  ROC = roc.curve(scores.class0 = outputs[(targets == 1)], 
                 scores.class1 = outputs[(targets == 0)], 
                 curve = T, max.compute = T, min.compute = T, rand.compute = T)
  PR = pr.curve(scores.class0 = outputs[(targets == 1)], 
                 scores.class1 = outputs[(targets == 0)], 
                 curve = T, max.compute = T, min.compute = T, rand.compute = T)
  if (returnRandom) {
    plot(ROC)
    plot(PR) 
  }
  # export
  temp = c(ROC$auc, PR$auc.integral, PR$rand[2])
  #cat('AUPR for a random classifier is ', temp[[3]], '\n')
  if (returnRandom) 
    return(c(temp[[1]], temp[[2]], temp[[3]]))
    else return(c(temp[[1]], temp[[2]]))
}
################################################# function: FastAUCTrain
fastAucScoreFunctionOnTrain <- function(covariateIndicesOfV, 
                                 targetNode = 'Target') {
  v = t[,c(covariateIndicesOfV,dim(test)[2]), drop = F]
  ## v: validation data set  
  junction = compile(as.grain(net)) 
  # Determine the real(target) and learned (outputs) values
  eval(parse(text = paste('targets = v$', targetNode, sep = "")))
  outputs = rep(0, length(targets))
  var = head(names(v), -1)
  values = matrix(nrow = 1, ncol = dim(v)[2]-1)
  ## for time
  time1 = Sys.time()
  numberOfTimesVisitedHere <<- 0
  ##
  for (i in 1:length(targets)) {
    #howMuchTimeLeft(i, length(targets), time1, numberOfTimesVisitedHere, 1/20)
    for (j in 1:(dim(v)[2]-1))
      values[j] = as.character(v[i, j])
    outputs[i] = querygrain(setEvidence(junction, 
                                        nodes = var, 
                                        states = values), 
                            nodes = targetNode, 
                            type = "joint")[2]
  }
  if (length(unique(outputs)) == 1)
    return(c(.5,NaN))
  #roc(targets, outputs, plot = TRUE)
  ROC = roc.curve(scores.class0 = outputs[(targets == 1)], 
                  scores.class1 = outputs[(targets == 0)])
  return (ROC$auc)
}
################################################# function: FastAUCTest
fastAucScoreFunctionOnTest <- function(covariateIndicesOfV, 
                                        targetNode = 'Target') {
  v = test[,c(covariateIndicesOfV,dim(test)[2]), drop = F]
  ## v: validation data set  
  junction = compile(as.grain(net)) 
  # Determine the real(target) and learned (outputs) values
  eval(parse(text = paste('targets = v$', targetNode, sep = "")))
  outputs = rep(0, length(targets))
  var = head(names(v), -1)
  values = matrix(nrow = 1, ncol = dim(v)[2]-1)
  ##
  for (i in 1:length(targets)) {
    #howMuchTimeLeft(i, length(targets), time1, numberOfTimesVisitedHere, 1/20)
    for (j in 1:(dim(v)[2]-1))
      values[j] = as.character(v[i, j])
    outputs[i] = querygrain(setEvidence(junction, 
                                        nodes = var, 
                                        states = values), 
                            nodes = targetNode, 
                            type = "joint")[2]
  }
  if (length(unique(outputs)) == 1)
    return(c(.5,NaN))
  #roc(targets, outputs, plot = TRUE)
  ROC = roc.curve(scores.class0 = outputs[(targets == 1)], 
                  scores.class1 = outputs[(targets == 0)])
  return (ROC$auc)
}
################################################# function: plot barchart with y grids
plotBarChart <- function(height, weight, label, rotate, filePath, save = FALSE) {
  margin = c(8,5,2,2)+.1 #(7,... #bot, left, top, right)
  offset = 3
  ylim = NULL
  fontSize = 1.5
  if (max(height) < .1)
  {fontSize = 2; offset = 4; margin = c(11,6,4,2)+.1}
  if (max(height) > .1) 
  {ylim = c(0,.3);   offset = 4}
  if (max(height) > .3) 
    ylim = c(0,1)
  par(mar=margin)
  p <- barplot(height, weight, col = colours()[3], ylab = "Probability of bloom", cex.lab = fontSize,
               #col = colours()[1:6], 
               axisnames = FALSE, space = .05, cex.axis = fontSize, ylim = ylim)
  text(p, par("usr")[3], labels = label, srt = rotate, adj= 1, xpd = TRUE, cex= fontSize, pos = 1, 
       offset = offset )
  grid(nx = NA, ny=NULL)
  if (save) {
    pdf(filePath)  
    par(mar=margin)
    p <- barplot(height, weight, col = colours()[3], ylab = "Probability of bloom", cex.lab = fontSize,
                 #col = colours()[1:6], 
                 axisnames = FALSE, space = .05, cex.axis = fontSize, ylim = ylim)
    par(family = "serif")
    text(p, par("usr")[3], labels = label, srt = rotate, adj= 1, xpd = TRUE, cex=fontSize, pos = 1, offset = offset)
    grid(nx = NA, ny=NULL)
    dev.off()
  }
}
################################################# function: train and test for lakes
trainAndTestPartitionForLakes = function(data, 
                                         ratioOfTestToOriginalData = .15,
                                         discretizationMethod = "hartemink",
                                         discritizationLevels = 3,
                                         readTrainAndTest = F) {
  if (readTrainAndTest) {
    t <- read.csv(file = paste(targetPath, "t.csv", sep = ""),
                  header = TRUE,
                  sep = ",")
    test <- read.csv(file = paste(targetPath, "test.csv", sep = ""),
                     header = TRUE,
                     sep = ",")
    # Factorize them
    for (i in 1:(dim(t)[2])) {
      t[, i] = factor(t[, i])
      test[,i] = factor(test[, i])
    }
  }  else {
    # Select %15 of the data for test
    numberOfValidationSamples = round(ratioOfTestToOriginalData*dim(data)[1])
    ## The data on the last month(s) of the last year(s) a lake is a valid candidate for the test dataset.
    ## We have two goals:
    ## 1) if an instance at time t appears in the test, no instance with a time later than t appears in train,
    ## 2) there should be some lakes that only appear in the test.
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
        cat(lake,'\n')
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
    test = data[indicesOfV,]
    t = data[setdiff(1:dim(data)[1], indicesOfV),]
    # delete month = 10
    t = t[t$Month != 10,]
    ## Report the proeperties of the test dataset
    reportDetailsOfTestDataset(t, test, targetPath)
    # Backup the train and test
    tOriginal = t
    testOriginal = test
    #########################  Discretize them 
    # All covariates should be numbers
    for (i in 1:dim(t)[2]) 
      t[,i] = as.numeric(t[,i])
    # noncovariates should not be discretized and Month cannot be divided into 50 quantiles
    ignoreVariableIndices = c(1,2,3,8,20)
    # print them
    cat(names(t)[ignoreVariableIndices])
    ########### First Method: Use hartemink's algorithm
    if (discretizationMethod == 'hartemink') {
      # Quantize just based on the train
      t[,-ignoreVariableIndices] = 
        discretize(t[,-ignoreVariableIndices],
                   method = "hartemink",
                   breaks = discritizationLevels, 
                   ibreaks = 50, 
                   order = T,
                   idisc = "interval")
    }
    ########### Second Method: equal frequency
    if (discretizationMethod == 'frequency') {
      # Quantize just based on the train
      t[,-ignoreVariableIndices] = 
        discretize(t[,-ignoreVariableIndices],
                   method = "quantile",
                   breaks = discritizationLevels, 
                   order = T)
    }
    
    ########### Third Method: equal interval lengths
    if (discretizationMethod == 'interval') {
      # Quantize just based on the train
      t[,-ignoreVariableIndices] = 
        discretize(t[,-ignoreVariableIndices],
                   method = "interval",
                   breaks = discritizationLevels, 
                   order = T)
    }
    ## Factorize the test based on the train
    quantizedTest = testOriginal
    # For each variable,
    for (i in c(1:dim(t)[2])[-ignoreVariableIndices]) {
      # and each instance
      for (instance in 1:dim(test)[1]) {
        ## find the interval in t that it belongs to
        interval = matrix(nrow = length(levels(t[,i])), ncol = 2)
        # checkout http://www.endmemo.com/program/R/gsub.php
        interval[,1] = as.numeric(gsub(".(-?[0-9\\.?0-9?e\\+?0-9]+) ?, ?(-?[0-9\\.?0-9?e\\+?0-9]+).", "\\1",  levels(t[,i])))
        interval[,2] = as.numeric(gsub(".(-?[0-9\\.?0-9?e\\+?0-9]+) ?, ?(-?[-?0-9\\.?0-9?e\\+?0-9]+).", "\\2",  levels(t[,i])))
        # which interval?
        for (j in 1:(dim(interval)[1]-1)) 
          if (testOriginal[instance,i] <= interval[j,2]) {
            quantizedTest[instance,i] = levels(t[,i])[j]
            break
          }
        j = dim(interval)[1]
        if (testOriginal[instance,i] >= interval[j,1]) 
          quantizedTest[instance,i] = levels(t[,i])[j]
      }
    }
    test = quantizedTest
    # Factorize them
    for (i in 1:(dim(t)[2])) {
      t[, i] = factor(t[, i])
      test[,i] = factor(test[, i])
      # levels(test[,i]) = levels(t[,i])
    }
    if (!((0 %in% test$Target) & (1 %in% test$Target))) 
      print('The test dataset does not include both values for the Target.\n')
    ### The rest is to check for special cases:
    ## Check if some levels are not observed, (no need to check the last one, i.e., Infestation)
    for (i in 5:(dim(t)[2])) {
      while (nlevels(t[, i]) > length(unique(t[, i]))) {
        # which level is not in the data?
        levelPresence = levels(t[,i])  %in% unique(t[,i])
        # find the first two empty levels
        firstEmptyLevel = match(FALSE, levelPresence)
        # If this is the only empty and last level, then
        secondEmptyLevel = firstEmptyLevel
        # Otherwise,
        # Check if the next level is empty and is not the last level
        while (!(levelPresence[secondEmptyLevel+1]) & (secondEmptyLevel+1 < length(levelPresence))) 
          # if so, then go to next level
          secondEmptyLevel = secondEmptyLevel + 1
        ## Make the new range
        # If the firstEmptyLevel is not the first level
        if (firstEmptyLevel == 1) {
          newRange = paste(sapply(strsplit(levels(t[, i])[firstEmptyLevel], ","), "[", 1),
                           ",",
                           sapply(strsplit(levels(t[, i])[secondEmptyLevel+1], ","), "[", 2))
          levelsToBeChanged = firstEmptyLevel:(secondEmptyLevel+1)
        } else{
          newRange = paste(sapply(strsplit(levels(t[, i])[firstEmptyLevel-1], ","), "[", 1),
                           ",",
                           sapply(strsplit(levels(t[, i])[secondEmptyLevel], ","), "[", 2))
          levelsToBeChanged = (firstEmptyLevel-1):secondEmptyLevel
        }
        eval(parse(text = paste(c('t','test'), "[,i] = as.character(", c('t','test'), "[,i])", sep = "")))
        for (j in levelsToBeChanged) {
          t[(t[,i] == levels(data[, i])[j]),i] = newRange
          test[(test[,i] == levels(data[, i])[j]),i] = newRange
        }
        t[, i] = factor(t[, i])
        test[, i] = factor(test[, i])
      }
    }
    for (i in 1:(dim(t)[2]-1)) {
      t[, i] = factor(t[, i])
      test[, i] = factor(test[, i])
    }
    ## Check if evey variable in both v and t has at least 2 levels
    for (i in 4:(dim(dataBeforeFactorizing)[2])) {
      if ((nlevels(test[, i]) == 1) | (nlevels(t[, i]) == 1)) {
        m = round(min(dataBeforeFactorizing[,i]),2)
        M = round(max(dataBeforeFactorizing[,i]),2)
        # Get the numerical version of v[, i] and t[, i]
        counter = 0
        for (lake in unique(dataBeforeFactorizing$lake)) {
          if (counter > numberOfValidationSamples)
            break
          for (year in unique(dataBeforeFactorizing[data$lake==lake,2])) {
            if (counter > numberOfValidationSamples)
              break
            counter = counter + 1
            # Assumption: the lakes are sorted in ascending order of first, sampled year and second, sampled month.
            indicesOfV[counter] = max(intersect(which(dataBeforeFactorizing$lake %in% lake), which(dataBeforeFactorizing$Year %in% year)))
          }
        }
        v_numeric_i = dataBeforeFactorizing[indicesOfV,i]
        t_numeric_i = dataBeforeFactorizing[setdiff(1:dim(data)[1], indicesOfV),i]
        # Find the min and max of the numerics
        m_v = round(min(v_numeric_i),2)
        M_v = round(max(v_numeric_i),2)
        m_t = round(min(t_numeric_i),2)
        M_t = round(max(t_numeric_i),2)
        middleValue = (max(m_v, m_t)+min(M_v, M_t))/2
        # Modify test
        test[, i] = as.character(test[, i])
        test[(v_numeric_i <=middleValue), i] = paste('(', m, ',', (m+M)/2, ']', sep = "")
        test[(v_numeric_i > middleValue), i] = paste('(', (m+M)/2, ',', M, ']', sep = "")
        # Modify t
        t[, i] = as.character(t[, i])
        t[(t_numeric_i <=middleValue), i] = paste('(', m, ',', (m+M)/2, ']', sep = "")
        t[(t_numeric_i > middleValue), i] = paste('(', (m+M)/2, ',', M, ']', sep = "")
      }
    }
    ## Check if any variable has fewer levels in t than in v
    for (i in 1:(dim(t)[2]-1)) 
      if (nlevels(t[,i]) < nlevels(test[,i]))
        cat(green(names(t)[i] ," has fewer levels in the training dataset than in the validation dataset.\n"))
    # Delete the year and pixel name
    # Factorize them
    for (i in 1:(dim(t)[2])) {
      t[, i] = factor(t[, i])
      test[,i] = factor(test[, i])
    }
    t[,c(1:2)] <- NULL
    test[,c(1:2)] <- NULL
    # Write
    if (T) {
      write.csv(t,
                file = paste(targetPath, 
                             "t.csv", sep = ""),
                row.names = FALSE)    
      write.csv(test,
                file = paste(targetPath,
                             "test.csv", sep = ""),
                row.names = FALSE) 
    }
  }
  # Print the training and testing datasets
  return(list(t,test))
}
################################################# function: plot and print a conditional probability
plotAndPrintAConditionalProbability <- function(net, targetNode, yIndices, zIndices, zlevelIndices, t,
                                                plotFilePath = "", tableFilePath = "", 
                                                plotSave = FALSE, tableSave = FALSE, 
                                                printTable = FALSE) {
  ## t: training data set: just to get the levels
  ## net: a class bn (from bnlearn) object  
  junction = compile(as.grain(net))
  # Get y
  y = nodes(net)[yIndices]
  # Which ones should be known variables?
  # do it only if zIndices is not empty
  if (zIndices != "") {
    z = nodes(net)[zIndices]
    # What levels?
    values = matrix(nrow = 1, ncol = length(z))
    counter = 0
    for (j in zIndices){
      counter = counter + 1
      values[counter] = as.character(levels(t[,j])[zlevelIndices])
    }
  } else {
    z = ""
    values = ""
  }
  outputs = querygrain(setEvidence(junction, 
                                   nodes = z, 
                                   states = values), 
                       nodes = c(targetNode, y), 
                       type = "conditional")
  # Get the length of the intervals of the "conditioned-on variable
  intervalsOfY = rownames(outputs)
  interval = matrix(nrow = length(intervalsOfY), ncol = 3)
  interval[,1] = as.numeric(gsub(".(-?[0-9\\.?0-9?e\\+?0-9]+) ?, ?(-?[0-9\\.?0-9?e\\+?0-9]+).", "\\1", intervalsOfY))
  interval[,2] = as.numeric(gsub(".(-?[0-9\\.?0-9?e\\+?0-9]+) ?, ?(-?[0-9\\.?0-9?e\\+?0-9]+).", "\\2", intervalsOfY))
  interval[,3] = interval[,2]-interval[,1]
  # round up interval
  interval = round(interval,1)
  # sort the intervals of the "conditioned-on" variable
  outputs = outputs[order(interval[,1]),]
  ##
  # Save the outputs
  if (tableSave) 
    print(xtable(outputs, digits = 3), file = paste(tableFilePath,
                                                    'PcondtionedOn', paste(y, collapse = "_"),
                                                    'given', paste(z, collapse = "_"),
                                                    ".tex", sep = "")) 
  # sort the interval
  interval = interval[order(interval[,1]),]
  # adjust the labels
  for (i in 1:(length(intervalsOfY)))
    if (interval[i,3]==0) {
      interval[i,3] = 1
      intervalsOfY[i] = paste(interval[i,1])
    } else {
      intervalsOfY[i] = paste('[', interval[i,1], ',', interval[i,2], ')', sep = "")
    }
  # How much should the lables rotate?
  rotation = 45
  # Treat the binary cases differently
  if (length(levels(t[,yIndices])) == 2) 
    if (all(levels(t[,yIndices]) == c("0","1"))) {
      intervalsOfY = c("Low","High")
      interval[,3] = matrix(c(1,1), nrow = 2)
      rotation = 0
    }
  # Treat the trinary cases differently
  if (length(levels(t[,yIndices])) == 3) 
    if (all(levels(t[,yIndices]) == c("1","2","3"))) {
      intervalsOfY = c("Low","Medium","High")
      interval[,3] = matrix(rep(1,dim(interval)[1],1), nrow = 1)
      rotation = 0
    }
  # plot the resulting query
  #?options("scipen"=-2, "digits"=1)
  plotBarChart(outputs[,2], interval[,3], intervalsOfY, rotation,
               paste(plotFilePath,
                     'PcondtionedOn', paste(y, collapse = "_"),
                     'given', paste(z, collapse = "_"),
                     ".pdf", sep = ""), 
               save = plotSave)
  # print the table
  if(printTable)
    outputs
}
################################################# function: plot/print all single node conditional probabilities
plotAndPrintAllSingleNodeConditionedProbabilities <- function(net, targetNode, t,
                                                              plotFilePath = "", tableFilePath = "", 
                                                              plotSave = FALSE, tableSave = FALSE) {
  ## v: validation data set  
  ## net: a class bn (from bnlearn) object  
  junction = compile(as.grain(net))
  # for all y except the target node
  for (yIndices in 1:(length(nodes(net))-1)) {
    # no z
    zIndices = ""
    # no levels?
    zlevelIndices = ""
    plotAndPrintAConditionalProbability(net, targetNode, yIndices, zIndices, zlevelIndices, t,
                                        plotFilePath = plotFilePath, 
                                        tableFilePath = tableFilePath, 
                                        plotSave = plotSave, tableSave = tableSave)
  
  }  
}
################################################# function: report details of test dataset
reportDetailsOfTestDataset = function(
  train, 
  test,
  targetPath = "",
  printToConsole = T
) {
  # This is the main text file
  text = c()
  ########################### Portion of train test
  dataSize = dim(train)[1] + dim(test)[1]
  temp = paste0('The dataset is partitioned into ', 
                 round(100*dim(train)[1]/dataSize), '% train and ',
                 round(100*dim(test)[1]/dataSize), '% test.')
  text = c(text, temp)
  ########################### Number of lakes and years
  numberOfLakesInTest = length(unique(test$lake))
  numberOfYearsInTest = length(unique(test$Year))
  temp = paste0('The test dataset consists of ', numberOfLakesInTest, " lakes.")
  text = c(text, temp)
  ########################### Number of unique lakes in test
  numberOfLakesUniqueToTest = length(setdiff(test$lake,train$lake))
  # How many instances to they make?
  numberOfInstancesforTheUniqueLakesInTest = 0
  for (lake in setdiff(test$lake,train$lake))
    numberOfInstancesforTheUniqueLakesInTest = numberOfInstancesforTheUniqueLakesInTest + 
    sum(test$lake == lake)
  temp = paste0('A number of ', numberOfLakesUniqueToTest, 
                'lakes are exclusive to the test dataset and do not appear in the train. ', 
                'The instances corresponding to these lakes form ', 
                round(100*numberOfInstancesforTheUniqueLakesInTest/dim(test)[1]), 
                '% of the test dataset.')
  text = c(text, temp)
  ########################### distribution of Month over the test instances
  # First, make sure they are not factors. 
  if (is.factor(test$Month)) {
    cat('Month in the test dataset is in factor form. Change it to integer.')
  } else {
    distributionOfMonthOverTestInstances = 
      round(100*hist(test$Month, 
                     breaks=seq(min(test$Month)-0.5, max(test$Month)+0.5, by=1))$density)
    distributionMatrix = 
      matrix(distributionOfMonthOverTestInstances, 
        nrow = 1,
        dimnames = list("distribution", c(min(test$Month):max(test$Month))))
    temp1 = paste0("The distribution of 'month' over the test instances",
                  "is as in Table~\ref{table:distributionOfMonthOverTestInstances}.")
    temp2 = print(xtable(distributionMatrix, 
                   digits = 0, 
                   label = 'table:distributionOfMonthOverTestInstances'))
    text = c(text, temp1, temp2, '%')
  }
  ########################### Append the texts.
  text = paste(text, collapse = "\n")
  # Print?
  if (printToConsole == T)
    cat(text)
  ########################### Save. 
  cat(text, file = paste0(targetPath, 'detailsOnTestDataset.tex'))
}
################################################# 