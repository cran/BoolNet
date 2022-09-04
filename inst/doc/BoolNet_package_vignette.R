### R code from vignette source 'BoolNet_package_vignette.Snw'

###################################################
### code chunk number 1: BoolNet_package_vignette.Snw:81-82 (eval = FALSE)
###################################################
## install.packages("BoolNet")


###################################################
### code chunk number 2: BoolNet_package_vignette.Snw:86-87
###################################################
library(BoolNet)


###################################################
### code chunk number 3: BoolNet_package_vignette.Snw:166-167 (eval = FALSE)
###################################################
## cellcycle <- loadNetwork("cellcycle.txt")


###################################################
### code chunk number 4: BoolNet_package_vignette.Snw:172-173
###################################################
data(cellcycle)


###################################################
### code chunk number 5: BoolNet_package_vignette.Snw:222-223
###################################################
data(yeastTimeSeries)


###################################################
### code chunk number 6: BoolNet_package_vignette.Snw:230-231
###################################################
binSeries <- binarizeTimeSeries(yeastTimeSeries)


###################################################
### code chunk number 7: BoolNet_package_vignette.Snw:236-239
###################################################
net <- reconstructNetwork(binSeries$binarizedMeasurements, 
                          method="bestfit", 
                          maxK=4)


###################################################
### code chunk number 8: BoolNet_package_vignette.Snw:245-246
###################################################
net


###################################################
### code chunk number 9: BoolNet_package_vignette.Snw:253-254 (eval = FALSE)
###################################################
## plotNetworkWiring(net)


###################################################
### code chunk number 10: BoolNet_package_vignette.Snw:268-272
###################################################
net <- reconstructNetwork(binSeries$binarizedMeasurements, 
                          method="bestfit", 
                          maxK=4, 
                          excludedDependencies = list("Sic1" = c("Sic1", "Fkh2")))


###################################################
### code chunk number 11: BoolNet_package_vignette.Snw:281-287
###################################################
net <- reconstructNetwork(binSeries$binarizedMeasurements, 
method="bestfit", maxK=4)
functionIndices <- c(1,2,3,2) #select function index for each regulatory component
dontCareDefaults <- lapply(seq_along(net$interactions), function(idx) rep(F, sum(net$interactions[[idx]][[functionIndices[idx]]]$func == -1))) #determine number of don't care values for each selected function and set them to 0
names(dontCareDefaults) <- net$genes
singleNet <- chooseNetwork(net, functionIndices, dontCareValues = dontCareDefaults)


###################################################
### code chunk number 12: BoolNet_package_vignette.Snw:290-291
###################################################
singleNet


###################################################
### code chunk number 13: BoolNet_package_vignette.Snw:297-298
###################################################
set.seed(3176)


###################################################
### code chunk number 14: BoolNet_package_vignette.Snw:300-304
###################################################
series <- generateTimeSeries(cellcycle, 
                             numSeries=100, 
                             numMeasurements=10, 
                             noiseLevel=0.1)


###################################################
### code chunk number 15: BoolNet_package_vignette.Snw:309-312
###################################################
binSeries <- binarizeTimeSeries(series, method="kmeans")
net <- reconstructNetwork(binSeries$binarizedMeasurements, method="bestfit")
net


###################################################
### code chunk number 16: BoolNet_package_vignette.Snw:318-319
###################################################
set.seed(4463)


###################################################
### code chunk number 17: BoolNet_package_vignette.Snw:321-326
###################################################
series <- generateTimeSeries(cellcycle, 
                             numSeries=10, 
                             numMeasurements=10, 
                             perturbations=1,
                             noiseLevel=0.1)


###################################################
### code chunk number 18: <
###################################################
series$perturbations


###################################################
### code chunk number 19: BoolNet_package_vignette.Snw:337-341
###################################################
perturbations <- series$perturbations
series$perturbations <- NULL

binSeries <- binarizeTimeSeries(series, method="kmeans")


###################################################
### code chunk number 20: BoolNet_package_vignette.Snw:344-348
###################################################
net <- reconstructNetwork(binSeries$binarizedMeasurements, 
                          method="bestfit", 
                          perturbations=perturbations)
net


###################################################
### code chunk number 21: BoolNet_package_vignette.Snw:357-358
###################################################
net <- generateRandomNKNetwork(n=10, k=3)


###################################################
### code chunk number 22: BoolNet_package_vignette.Snw:361-362
###################################################
net <- generateRandomNKNetwork(n=10, k=c(1,2,3,1,3,2,3,2,1,1))


###################################################
### code chunk number 23: BoolNet_package_vignette.Snw:367-368
###################################################
net <- generateRandomNKNetwork(n=20, k=20, topology="scale_free")


###################################################
### code chunk number 24: BoolNet_package_vignette.Snw:372-373
###################################################
net <- generateRandomNKNetwork(n=10, k=3, linkage="lattice")


###################################################
### code chunk number 25: BoolNet_package_vignette.Snw:378-382
###################################################
net <- generateRandomNKNetwork(n=10, 
                               k=3, 
                               functionGeneration="biased", 
                               zeroBias=0.75)


###################################################
### code chunk number 26: BoolNet_package_vignette.Snw:389-397
###################################################
net1 <- generateRandomNKNetwork(n=10, 
                                k=3,
                                functionGeneration=generateCanalyzing, 
                                zeroBias=0.75)
net2 <- generateRandomNKNetwork(n=10, 
                                k=3, 
                                functionGeneration=generateNestedCanalyzing, 
                                zeroBias=0.75)


###################################################
### code chunk number 27: BoolNet_package_vignette.Snw:405-418
###################################################
isMonotone <- function(input, func)
{
  for (i in seq_len(ncol(input)))
  # check each input gene
  {
    groupResults <- split(func, input[,i])
    if (any(groupResults[[1]] < groupResults[[2]]) && 
        any(groupResults[[1]] > groupResults[[2]]))
      # the function is not monotone
      return(FALSE)
  }
  return(TRUE)
}


###################################################
### code chunk number 28: BoolNet_package_vignette.Snw:427-431
###################################################
net <- generateRandomNKNetwork(n=10,
                               k=3, 
                               validationFunction="isMonotone", 
                               failureIterations=1000)


###################################################
### code chunk number 29: BoolNet_package_vignette.Snw:442-444
###################################################
data(cellcycle)
knockedOut <- fixGenes(cellcycle, "CycD", 0)


###################################################
### code chunk number 30: BoolNet_package_vignette.Snw:447-448
###################################################
knockedOut <- fixGenes(cellcycle, 1, 0)


###################################################
### code chunk number 31: BoolNet_package_vignette.Snw:451-452
###################################################
overExpressed <- fixGenes(cellcycle, "CycD", 1)


###################################################
### code chunk number 32: BoolNet_package_vignette.Snw:455-456
###################################################
originalNet <- fixGenes(knockedOut, "CycD", -1)


###################################################
### code chunk number 33: BoolNet_package_vignette.Snw:461-462
###################################################
newNet <- fixGenes(cellcycle, c("CycD","CycE"), c(0,1))


###################################################
### code chunk number 34: BoolNet_package_vignette.Snw:473-475
###################################################
data(cellcycle)
stateTransition(cellcycle, rep(1,10))


###################################################
### code chunk number 35: BoolNet_package_vignette.Snw:479-481
###################################################
path <- getPathToAttractor(cellcycle, rep(0,10))
path


###################################################
### code chunk number 36: BoolNet_package_vignette.Snw:486-487 (eval = FALSE)
###################################################
## plotSequence(sequence=path)


###################################################
### code chunk number 37: BoolNet_package_vignette.Snw:500-501 (eval = FALSE)
###################################################
## sequenceToLaTeX(sequence=path, file="sequence.tex")


###################################################
### code chunk number 38: BoolNet_package_vignette.Snw:506-508
###################################################
startState <- generateState(cellcycle, specs=c("CycD"=1,"CycA"=1))
stateTransition(cellcycle,startState)


###################################################
### code chunk number 39: BoolNet_package_vignette.Snw:515-516
###################################################
data(igf)


###################################################
### code chunk number 40: BoolNet_package_vignette.Snw:521-523
###################################################
startState <- generateState(igf, specs=c("IGF"=1))
stateTransition(igf, startState)


###################################################
### code chunk number 41: BoolNet_package_vignette.Snw:526-527
###################################################
getPathToAttractor(network=igf,state=startState)


###################################################
### code chunk number 42: BoolNet_package_vignette.Snw:532-535
###################################################
startState <- generateState(igf, specs=list("IGF"=c(0,0,1)))

startState


###################################################
### code chunk number 43: BoolNet_package_vignette.Snw:545-546 (eval = FALSE)
###################################################
## plotSequence(network=igf, startState=startState)


###################################################
### code chunk number 44: BoolNet_package_vignette.Snw:552-553
###################################################
set.seed(54321)


###################################################
### code chunk number 45: BoolNet_package_vignette.Snw:555-556
###################################################
stateTransition(cellcycle, rep(1,10), type="asynchronous")


###################################################
### code chunk number 46: BoolNet_package_vignette.Snw:562-563
###################################################
set.seed(4321)


###################################################
### code chunk number 47: BoolNet_package_vignette.Snw:565-567
###################################################
stateTransition(cellcycle, rep(1,10), type="asynchronous", 
geneProbabilities=c(0.05,0.05,0.2,0.3,0.05,0.05,0.05,0.05,0.1,0.1))


###################################################
### code chunk number 48: BoolNet_package_vignette.Snw:574-576
###################################################
stateTransition(cellcycle, rep(1,10), type="asynchronous", 
chosenGene="CycE")


###################################################
### code chunk number 49: BoolNet_package_vignette.Snw:580-581
###################################################
set.seed(432)


###################################################
### code chunk number 50: BoolNet_package_vignette.Snw:583-585
###################################################
data(examplePBN)
stateTransition(examplePBN, c(0,1,1), type="probabilistic")


###################################################
### code chunk number 51: BoolNet_package_vignette.Snw:588-590
###################################################
stateTransition(examplePBN, c(0,1,1), type="probabilistic", 
chosenFunctions=c(2,1,2))


###################################################
### code chunk number 52: BoolNet_package_vignette.Snw:608-611 (eval = FALSE)
###################################################
## data(cellcycle)
## attr <- getAttractors(cellcycle)
## attr


###################################################
### code chunk number 53: BoolNet_package_vignette.Snw:614-616
###################################################
attr <- getAttractors(cellcycle)
attr


###################################################
### code chunk number 54: BoolNet_package_vignette.Snw:622-623 (eval = FALSE)
###################################################
## print(attr, activeOnly=TRUE)


###################################################
### code chunk number 55: BoolNet_package_vignette.Snw:626-627
###################################################
print(attr, activeOnly=TRUE)


###################################################
### code chunk number 56: BoolNet_package_vignette.Snw:636-637
###################################################
getAttractorSequence(attr, 2)


###################################################
### code chunk number 57: BoolNet_package_vignette.Snw:644-646 (eval = FALSE)
###################################################
## tt <- getTransitionTable(attr)
## tt


###################################################
### code chunk number 58: BoolNet_package_vignette.Snw:660-661 (eval = FALSE)
###################################################
## getBasinOfAttraction(attr, 1)


###################################################
### code chunk number 59: BoolNet_package_vignette.Snw:666-667 (eval = FALSE)
###################################################
## getStateSummary(attr, c(1,1,1,1,1,1,1,1,1,1))


###################################################
### code chunk number 60: BoolNet_package_vignette.Snw:679-680 (eval = FALSE)
###################################################
## plotStateGraph(attr)


###################################################
### code chunk number 61: BoolNet_package_vignette.Snw:686-687 (eval = FALSE)
###################################################
## plotStateGraph(attr, piecewise=TRUE)


###################################################
### code chunk number 62: BoolNet_package_vignette.Snw:702-703 (eval = FALSE)
###################################################
## attr <- getAttractors(cellcycle, method="random", startStates=100)


###################################################
### code chunk number 63: BoolNet_package_vignette.Snw:707-710 (eval = FALSE)
###################################################
## attr <- getAttractors(cellcycle, 
##                       method="chosen", 
##                       startStates=list(rep(0,10),rep(1,10)))


###################################################
### code chunk number 64: BoolNet_package_vignette.Snw:719-720 (eval = FALSE)
###################################################
## plotAttractors(attr, subset=2)


###################################################
### code chunk number 65: BoolNet_package_vignette.Snw:723-724 (eval = FALSE)
###################################################
## attractorsToLaTeX(attr, subset=2, file="attractors.tex")


###################################################
### code chunk number 66: BoolNet_package_vignette.Snw:738-742
###################################################
attr <- getAttractors(cellcycle, 
                      type="asynchronous",
                      method="random", 
                      startStates=500)


###################################################
### code chunk number 67: BoolNet_package_vignette.Snw:748-749 (eval = FALSE)
###################################################
## attr


###################################################
### code chunk number 68: BoolNet_package_vignette.Snw:775-780 (eval = FALSE)
###################################################
## attr <- getAttractors(cellcycle, 
##                       type="asynchronous",
##                       method="random", 
##                       startStates=500, 
##                       avoidSelfLoops=FALSE)


###################################################
### code chunk number 69: BoolNet_package_vignette.Snw:795-796 (eval = FALSE)
###################################################
## plotAttractors(attr, subset=2, mode="graph", drawLabels=FALSE)


###################################################
### code chunk number 70: BoolNet_package_vignette.Snw:804-806 (eval = FALSE)
###################################################
## sim <- simulateSymbolicModel(igf)
## sim


###################################################
### code chunk number 71: BoolNet_package_vignette.Snw:809-811
###################################################
sim <- simulateSymbolicModel(igf)
sim


###################################################
### code chunk number 72: BoolNet_package_vignette.Snw:825-826 (eval = FALSE)
###################################################
## plotAttractors(sim, subset=2)


###################################################
### code chunk number 73: BoolNet_package_vignette.Snw:829-830 (eval = FALSE)
###################################################
## plotStateGraph(sim)


###################################################
### code chunk number 74: BoolNet_package_vignette.Snw:844-845
###################################################
set.seed(43851)


###################################################
### code chunk number 75: BoolNet_package_vignette.Snw:847-848
###################################################
sim <- simulateSymbolicModel(igf, method="random", startStates=2)


###################################################
### code chunk number 76: BoolNet_package_vignette.Snw:851-852
###################################################
sim$sequences


###################################################
### code chunk number 77: BoolNet_package_vignette.Snw:867-870
###################################################
data(examplePBN)
sim <- markovSimulation(examplePBN)
sim


###################################################
### code chunk number 78: BoolNet_package_vignette.Snw:876-877 (eval = FALSE)
###################################################
## plotPBNTransitions(sim)


###################################################
### code chunk number 79: BoolNet_package_vignette.Snw:890-895
###################################################
data(cellcycle)
sim <- markovSimulation(cellcycle, 
                        numIterations=1024, 
                        returnTable=FALSE)
sim


###################################################
### code chunk number 80: BoolNet_package_vignette.Snw:904-909
###################################################
sim <- markovSimulation(cellcycle, 
                        numIterations=1024,
                        returnTable=FALSE, 
                        startStates=list(rep(1,10)))
sim


###################################################
### code chunk number 81: BoolNet_package_vignette.Snw:921-922
###################################################
set.seed(3361)


###################################################
### code chunk number 82: BoolNet_package_vignette.Snw:924-929
###################################################
data(cellcycle)
r <- perturbTrajectories(cellcycle, 
                         measure="hamming", 
                         numSamples=100, 
                         flipBits=1)


###################################################
### code chunk number 83: BoolNet_package_vignette.Snw:933-934
###################################################
r$value


###################################################
### code chunk number 84: BoolNet_package_vignette.Snw:938-944
###################################################
r <- perturbTrajectories(cellcycle, 
                         measure="sensitivity", 
                         numSamples=100, 
                         flipBits=1, 
                         gene="CycE")
r$value


###################################################
### code chunk number 85: BoolNet_package_vignette.Snw:949-954
###################################################
r <- perturbTrajectories(cellcycle, 
                         measure="attractor", 
                         numSamples=100, 
                         flipBits=1)
r$value


###################################################
### code chunk number 86: BoolNet_package_vignette.Snw:961-964
###################################################
perturbedNet <- perturbNetwork(cellcycle, 
                               perturb="functions", 
                               method="bitflip")


###################################################
### code chunk number 87: BoolNet_package_vignette.Snw:969-972
###################################################
perturbedNet <- perturbNetwork(cellcycle, 
                               perturb="functions", 
                               method="shuffle")


###################################################
### code chunk number 88: BoolNet_package_vignette.Snw:978-982
###################################################
perturbedNet <- perturbNetwork(cellcycle, 
                               perturb="transitions", 
                               method="bitflip", 
                               numStates=10)


###################################################
### code chunk number 89: BoolNet_package_vignette.Snw:990-1046 (eval = FALSE)
###################################################
## # Perform a robustness test on a network
## # by counting the numbers of perturbed networks
## # containing the attractors of the original net
## 
## library(BoolNet)
## 
## # load mammalian cell cycle network
## data(cellcycle)
## 
## # get attractors in original network
## attrs <- getAttractors(cellcycle, canonical=TRUE)
## 
## # create 1000 perturbed copies of the network and search for attractors
## perturbationResults <- sapply(1:1000, function(i)
## {
##   # perturb network and identify attractors
##   perturbedNet <- perturbNetwork(cellcycle, perturb="functions", method="bitflip")
##   perturbedAttrs <- getAttractors(perturbedNet, canonical=TRUE)
##   
##   # check whether the attractors in the original network exist in the perturbed network
##   attractorIndices <- sapply(attrs$attractors,function(attractor1)
##         {
##           index <- which(sapply(perturbedAttrs$attractors, function(attractor2)
##             {
##               identical(attractor1, attractor2)
##             }))
##           if (length(index) == 0)
##             NA
##           else
##             index
##         })
##   return(attractorIndices)
## })
## 
## # perturbationResults now contains a matrix
## # with the first 2 columns specifying the indices or the 
## # original attractors in the perturbed network 
## # (or NA if the attractor was not found) and the next 2 
## # columns counting the numbers of states
## # in the basin of attraction (or NA if the attractor was not found)
## 
## # measure the total numbers of occurrences of the original attractors in the perturbed copies
## numOccurrences <- apply(perturbationResults[seq_along(attrs$attractors),,drop=FALSE], 1,
##                       function(row)sum(!is.na(row)))
## 
## # print original attractors
## cat("Attractors in original network:\n")
## print(attrs)
## 
## # print information
## cat("Number of occurrences of the original attractors",
##   "in 1000 perturbed copies of the network:\n")
## for (i in 1:length(attrs$attractors))
## {
##   cat("Attractor ",i,": ",numOccurrences[i],"\n",sep="")
## }


###################################################
### code chunk number 90: BoolNet_package_vignette.Snw:1080-1085 (eval = FALSE)
###################################################
## data(cellcycle)
## res <- testNetworkProperties(cellcycle, 
##                       numRandomNets=100, 
##                       testFunction="testAttractorRobustness", 
##                       testFunctionParams = list(copies=100, perturb="functions"))


###################################################
### code chunk number 91: BoolNet_package_vignette.Snw:1109-1114 (eval = FALSE)
###################################################
## testNetworkProperties(cellcycle, 
##                       numRandomNets=100,
##                       testFunction="testTransitionRobustness",
##                       testFunctionParams=list(numSamples=100),
##                       alternative="less")  


###################################################
### code chunk number 92: BoolNet_package_vignette.Snw:1121-1124 (eval = FALSE)
###################################################
## testNetworkProperties(cellcycle, 
##                       numRandomNets=100, 
##                       testFunction="testIndegree")


###################################################
### code chunk number 93: BoolNet_package_vignette.Snw:1142-1146 (eval = FALSE)
###################################################
## testNetworkProperties(cellcycle, 
##                       numRandomNets=100, 
##                       testFunction="testIndegree", 
##                       accumulation="kullback_leibler")


###################################################
### code chunk number 94: BoolNet_package_vignette.Snw:1162-1174
###################################################
testBasinSizes <- function(network, accumulate=TRUE, params)
{
  attr <- getAttractors(network)
  basinSizes <- sapply(attr$attractors, function(a)
                      { 
                         a$basinSize
                      })
   if (accumulate)
     return(mean(basinSizes))
   else
     return(basinSizes)                  
}


###################################################
### code chunk number 95: BoolNet_package_vignette.Snw:1180-1184 (eval = FALSE)
###################################################
## testNetworkProperties(cellcycle, 
##                       numRandomNets=100,
##                       testFunction="testBasinSizes",
##                       xlab="Average size of basins of attraction")


###################################################
### code chunk number 96: BoolNet_package_vignette.Snw:1207-1208 (eval = FALSE)
###################################################
## saveNetwork(cellcycle, file="cellcycle.txt")


###################################################
### code chunk number 97: BoolNet_package_vignette.Snw:1213-1215
###################################################
net <- generateRandomNKNetwork(n=10, k=3, readableFunctions=FALSE)
saveNetwork(net, file="randomnet.txt", generateDNF=TRUE)


###################################################
### code chunk number 98: BoolNet_package_vignette.Snw:1230-1233 (eval = FALSE)
###################################################
## toSBML(cellcycle, file="cellcycle.sbml")
## sbml_cellcycle <- loadSBML("cellcycle.sbml")
## sbml_cellcycle


###################################################
### code chunk number 99: BoolNet_package_vignette.Snw:1252-1253 (eval = FALSE)
###################################################
## system.file("doc/example.btp", package="BoolNet")


###################################################
### code chunk number 100: BoolNet_package_vignette.Snw:1285-1286
###################################################
net <- loadBioTapestry(system.file("doc/example.btp", package="BoolNet"))


###################################################
### code chunk number 101: BoolNet_package_vignette.Snw:1288-1289 (eval = FALSE)
###################################################
## net <- loadBioTapestry("example.btp")


###################################################
### code chunk number 102: BoolNet_package_vignette.Snw:1294-1295 (eval = FALSE)
###################################################
## net


###################################################
### code chunk number 103: BoolNet_package_vignette.Snw:1318-1319 (eval = FALSE)
###################################################
## plotNetworkWiring(net)


###################################################
### code chunk number 104: BoolNet_package_vignette.Snw:1340-1343 (eval = FALSE)
###################################################
## data(cellcycle)
## attr <- getAttractors(cellcycle)
## toPajek(attr, file="cellcycle.net")


###################################################
### code chunk number 105: BoolNet_package_vignette.Snw:1347-1348 (eval = FALSE)
###################################################
## toPajek(attr, file="cellcycle.net", includeLabels=TRUE)


