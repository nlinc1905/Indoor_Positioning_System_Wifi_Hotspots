
setwd("C:/Users/Nick/Documents/R Directory/CRAWDAD IPS Data")

#------------------------------------------------------------------#
#----------------------Load and Clean Data-------------------------#
#------------------------------------------------------------------#

#Check for comment lines in the data
txt <- readLines("offline.final.trace.txt")
sum(substr(txt, 1, 1) == "#")
#So 5312 comments

head(txt)

#Build function to split the text on certain characters into rows, then piece together column names and rows of data
processLine <- function(x) {
  #Regex in the strsplit function to split data on ;=, characters (determined by looking at data)
  tokens = strsplit(x, "[;=,]")[[1]]
  #If no signals are recorded (token length is 10) then remove the row
  if (length(tokens) == 10)
    return(NULL)
  #For each signal recording, tokens 1, 3, 5, and 9 are column names
  tmp = matrix(tokens[-(1:10)], ncol=4, byrow=T)
  #Column bind column names with the other rows in matrix form
  cbind(matrix(tokens[c(2,4,6:8,10)], nrow=nrow(tmp), ncol=6, byrow=T), tmp)
}

#Run the processLine function over the entire dataset to build dataframe
lines <- txt[substr(txt, 1,1) != "#"]  #Removes comments from data
tmp <- lapply(lines, processLine)
train <- as.data.frame(do.call("rbind", tmp), stringsAsFactors=F)

#Check dimensions of training set
dim(train)
head(train)
str(train)

#Add variable names
names(train) <- c("time", "scanMac", "posX", "posY", "posZ", "orientation", "mac", "signal", "channel", "type")

#Convert time, position, orientation, and signal to numeric variables
varList <- c("time", "posX", "posY", "posZ", "orientation", "signal")
train[varList] <- lapply(train[varList], as.numeric)

#Only the hotspot measurements are wanted, and they are found by type=3
train <- train[train$type=="3",]
#Safe to remove the type variable since it has only 1 value now
train$type <- NULL

head(train)
dim(train)

#Convert time from milliseconds to seconds
train$time <- train$time/1000
class(train$time) <- c("POSIXt", "POSIXct")

summary(train[, varList])
summary(sapply(train[,c("mac", "channel", "scanMac")], as.factor))

#Since the scan device and the elevation (posZ) are constant, they can be dropped
train$scanMac <- NULL
train$posZ <- NULL
#Can also get rid of unneeded objects to clear up space
remove(lines)
remove(txt)

#Should have 8 values for orientation according to documentation, verify it
length(unique(train$orientation))
#203 different orientations but should only have 8, check distribution
plot(ecdf(train$orientation))
#Orientations are close to the 8 expected but are too precise, need to group them
roundOrientation <- function(angles) {
  refs = seq(0, by=45, length=9)
  q = sapply(angles, function(o) which.min(abs(o-refs)))
  c(refs[1:8], 0)[q]
}
#Create new orientations
train$angle <- roundOrientation(train$orientation)
plot(ecdf(train$angle))

#From summary data, it seems that number of MAC addreses = number of channels, so check if there's 1:1 mapping
c(length(unique(train$mac)), length(unique(train$channel)))
#Only 6 hotspots should appear but there are 12 access points and 8 channels, meaning there is noise
table(train$mac)
#5 of these mac addresses appear much less frequently than the rest - they're probably noise
goodMacs <- names(sort(table(train$mac), decreasing=T))[1:7] #keep the top 7 Macs
train <- train[train$mac %in% goodMacs,]

#Check remaining Macs to ensure there is a non-zero count for each combo
apply(with(train, table(mac, channel)), 1, function(x) sum(x >0))
#Since there is a 1:1 relationship between remaining Macs and channels, remove channel
train$channel <- NULL

#Remove nulls from positions
locDF <- with(train, by(train, list(posX, posY), function(x) x)) #Create data frame to store positions
locDF <- locDF[ !sapply(locDF, is.null)]
length(locDF)

#------------------------------------------------------------------#
#---------------------Exploratory Analysis-------------------------#
#------------------------------------------------------------------#

#Examine signal strength at different orientations and at different access points to see if the distributions are normal
library(lattice)
#Use lattice box plots to examine distributions for orienatations and locations
bwplot(signal ~ factor(angle) | mac, data=train, 
       subset = posX == 2 & posY == 12 & mac != "00:0f:a3:39:dd:cd",
       layout = c(2,3))
#See how signal strength varies with orientation?  Note that the 1 extra Mac address was dropped

#Check normality with density plots
densityplot( ~ signal | mac + factor(angle), data=train,
             subset = posX == 24 & posY ==4 & mac != "00:0f:a3:39:dd:cd",
             bw = 0.5, plot.points = F)

#Create a list of data frames for every combo of (x,y), angle, and access point
train$posXY <- paste(train$posX, train$posY, sep="-")  #Factor with all unique combos of x and y
byLocAngleAP <- with(train, by(train, list(posXY, angle, mac), function(x) x))
#Get summary stats on each frame
signalSummary <- lapply(byLocAngleAP,
                        function (oneLoc) {
                          ans = oneLoc[1, ]
                          ans$medSignal = median(oneLoc$signal)
                          ans$avgSignal = mean(oneLoc$signal)
                          ans$num = length(oneLoc$Signal)
                          ans$sdSignal = sd(oneLoc$signal)
                          ans$iqrSignal = IQR(oneLoc$signal)
                          ans
                        })
trainSummary <- do.call("rbind", signalSummary)

#By looking at SD, it appears SD increases with average signal strength
breaks <- seq(-90, -30, by = 5)
bwplot(sdSignal ~ cut(avgSignal, breaks = breaks), data=trainSummary,
       subset = mac != "00:0f:a3:39:dd:cd",
       xlab = "Mean Signal", ylab = "SD Signal")

#Look at skew by plotting difference (avg - med) against the number of observations in scatterplot with fitted line
with(trainSummary, smoothScatter((avgSignal - medSignal) ~ num,
                                 xlab = "Number of Observations",
                                 ylab = "mean - median"))
abline(h=0, col="#984ea3", lwd=2)
lo.obj <- with(trainSummary, loess(diff ~ num, 
                                   data = data.frame(diff = (avgSignal - medSignal),
                                                     num = num)))
lo.obj.pr <- predict(lo.obj, newdata = data.frame(num = (70:120)))
lines(x=70:120, y=lo.obj.pr, col="#4daf4a", lwd=2)

#View signal strength by location (pick an angle, say 0 degrees, and view a topographic heat map of signal strength)
oneAPAngle <- subset(trainSummary, mac == goodMacs[5] & angle == 0)
library(fields)  #Used for drawing heatmaps
smoothSS <- Tps(oneAPAngle[, c("posX", "posY")], oneAPAngle$avgSignal)
#Predict the value for a fitted surface at a grid of observed positions
vizSmooth <- predictSurface(smoothSS)
#Plot the predicted signal trength
plot.surface(vizSmooth, type = "C")
#Add locations where the measurements were taken
points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex = 0.5)

#Wrap all this into a function so that you can draw heat map for any angle and Mac address
surfaceSS <- function(d, m, a) {
  oneAPAngle = subset(d, mac == goodMacs[m] & angle == a)
  library(fields)
  smoothSS <- Tps(oneAPAngle[, c("posX", "posY")], oneAPAngle$avgSignal)
  vizSmooth <- predictSurface(smoothSS)
  plot.surface(vizSmooth, type = "C")
  points(oneAPAngle$posX, oneAPAngle$posY, pch=19, cex=0.5)
}

#Tell R to plot matrix
parCur <- par(mfrow=c(2,2), mar=rep(1,4))

#Call surfaceSS 4 times to draw 4 heatmaps
mapply(surfaceSS, m = rep(c(5, 1), each = 2), a=rep(c(0, 135), 2),
       d=list(data=trainSummary))
#Top left = mac 1 at angle 0, top right = mac 1 at angle 135, bottom is mac 2 at same angles

#Reset plotting parameters for future plots
par(parCur)

#Store the positions of all 6 access points (WiFi hotspots) since they are known
ap <- matrix( c(7.5, 6.3, 2.5, -0.8, 12.8, -2.8, 1, 14, 33.5, 9.3, 33.5, 2.8),
              ncol=2, byrow=T, dimnames = list(goodMacs[-2], c("x", "y")))
#Mac addresses are the row names
ap

#Remove the additional access point (the one that's not a hotspot) from trainSummary
trainSummary <- subset(trainSummary, mac != goodMacs[2])

#Calculate Euclidean distances from the device to the hotspots, starting with x coordinate and then moving on to y
diffs <- trainSummary[ , c("posX", "posY")] - ap[trainSummary$mac, ]
trainSummary$dist <- sqrt(diffs[,1]^2 + diffs[,2]^2)

#Scatter plot of hotspots
xyplot(signal ~ dist | factor(mac) + factor(angle), 
       data = trainSummary, pch=19, cex=0.3, xlab="distance")
#Note plot curvature - log transformation might help

#------------------------------------------------------------------#
#--------------------Load and Clean Test Data----------------------#
#------------------------------------------------------------------#

#Clean and load test data (the online dataset)

#Comments check
txt <- readLines("online.final.trace.txt")
sum(substr(txt, 1, 1) == "#")  #240 comments

processLine <- function(x) {
  #Regex in the strsplit function to split data on ;=, characters (determined by looking at data)
  tokens = strsplit(x, "[;=,]")[[1]]
  #If no signals are recorded (token length is 10) then remove the row
  if (length(tokens) == 10)
    return(NULL)
  #For each signal recording, tokens 1, 3, 5, and 9 are column names
  tmp = matrix(tokens[-(1:10)], ncol=4, byrow=T)
  #Column bind column names with the other rows in matrix form
  cbind(matrix(tokens[c(2,4,6:8,10)], nrow=nrow(tmp), ncol=6, byrow=T), tmp)
}

lines <- txt[substr(txt, 1,1) != "#"]  #Removes comments from data
tmp <- lapply(lines, processLine)
test <- as.data.frame(do.call("rbind", tmp), stringsAsFactors=F)
names(test) <- c("time", "scanMac", "posX", "posY", "posZ", "orientation", "mac", "signal", "channel", "type")

varList <- c("time", "posX", "posY", "posZ", "orientation", "signal")
test[varList] <- lapply(test[varList], as.numeric)

#Only the hotspot measurements are wanted, and they are found by type=3
test <- test[test$type=="3",]
#Safe to remove the type variable since it has only 1 value now
test$type <- NULL
test$time <- test$time/1000
class(test$time) <- c("POSIXt", "POSIXct")
test$scanMac <- NULL
test$posZ <- NULL
remove(lines)
remove(txt)

roundOrientation <- function(angles) {
  refs = seq(0, by=45, length=9)
  q = sapply(angles, function(o) which.min(abs(o-refs)))
  c(refs[1:8], 0)[q]
}
#Create new orientations
test$angle <- roundOrientation(test$orientation)

#Keep the same 6 Macs as in trainSummary
table(trainSummary$mac)
goodMacs <- names(sort(table(test$mac), decreasing=T))[2:7] #keep the sam
test <- test[test$mac %in% goodMacs,]

#Check remaining Macs to ensure there is a non-zero count for each combo
apply(with(test, table(mac, channel)), 1, function(x) sum(x >0))
#Since there is a 1:1 relationship between remaining Macs and channels, remove channel
test$channel <- NULL

#Remove nulls from positions
locDF <- with(test, by(test, list(posX, posY), function(x) x)) #Create data frame to store positions
locDF <- locDF[ !sapply(locDF, is.null)]
length(locDF)
head(test)

#This is where the similarities between loading offline and online sets end

#Create xy position
test$posXY = paste(test$posX, test$posY, sep="-")
length(unique(test$posXY))  #Confirmed 60 locations
#Tally signal strengths recorded at each location
tabtestXYA = table(test$posXY, test$angle)
tabtestXYA[1:6,]
#Signal strengths were recorded at 1 orientation per location

#Restructure test data to have 6 columns of signal strengths (one for each hotspot)
keepVars <- c("posXY", "posX", "posY", "orientation", "angle")
byLoc <- with(test,
              by(test, list(posXY),
                 function(x) {
                   ans = x[1, keepVars]
                   avgSS = tapply(x$signal, x$mac, mean)
                   y = matrix(avgSS, nrow=1, ncol=6,
                              dimnames = list(ans$posXY, names(avgSS)))
                   cbind(ans, y)
                 }))
testSummary <- do.call("rbind", byLoc)
dim(testSummary)  #Confirm 60 rows and 11 columns

#Since orientation DOES impact signal strength (refer to graph way up above), must match orienation in training and test sets

#For m angles, find the closest desired orientations to the new observations
m <- 3
angleNewObs <- 230
refs <- seq(0, by=45, length=8)
nearestAngle <- roundOrientation(angleNewObs)

if ( m%% 2 == 1) {
  angles = seq(-45*(m-1)/2, 45*(m-1)/2, length=m)
} else {
  m = m + 1
  angles = seq(-45*(m-1)/2, 45*(m-1)/2, length=m)
  if (sign(angleNewObs - nearestAngle) > -1)
    angles = angles[-1]
  else
    angles = angles[-m]
}

#Map the angles to values in refs (-45 maps to 335 and 405 maps to 45)
angles <- angles + nearestAngle
angles[angles < 0] <- angles [angles <0] + 360
angles[angles > 360] <- angles[angles >360] - 360

#Pick observations from trinSummary to analyze and build model with
trainSubset <- trainSummary[trainSummary$angle %in% angles, ]

#Aggregate the signal stengths from these angles and create a data structure like testSummary, using resapheSS() function
reshapeSS <- function (data, varSignal = "signal", keepVars = c("posXY", "posX", "posY")) {
  byLocation = with(data, by(data, list(posXY),
                             function(x) {
                               ans = x[1, keepVars]
                               avgSS = tapply(x[,varSignal], x$mac, mean)
                               y = matrix(avgSS, nrow=1, ncol=6, dimnames = list(ans$posXY, names(avgSS)))
                               cbind(ans, y)
                             }))
  newDataSS <- do.call("rbind", byLocation)
  return(newDataSS)
}

#Summarize and reshape trainSubset
trainSS <- reshapeSS(trainSubset, varSignal = "avgSignal")

#Wrap everything that was just done into a function and test it with angle 130
selectTrain = function(angleNewObs, signals = NULL, m = 1){
  # m is the number of angles to keep between 1 and 5
  refs = seq(0, by = 45, length  = 8)
  nearestAngle = roundOrientation(angleNewObs)
  
  if (m %% 2 == 1) 
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
  else {
    m = m + 1
    angles = seq(-45 * (m - 1) /2, 45 * (m - 1) /2, length = m)
    if (sign(angleNewObs - nearestAngle) > -1) 
      angles = angles[ -1 ]
    else 
      angles = angles[ -m ]
  }
  angles = angles + nearestAngle
  angles[angles < 0] = angles[ angles < 0 ] + 360
  angles[angles > 360] = angles[ angles > 360 ] - 360
  angles = sort(angles) 
  
  trainSubset = signals[ signals$angle %in% angles, ]
  reshapeSS(trainSubset, varSignal = "avgSignal")
}

#Test function with new obs angle 130
train130 <- selectTrain(130, trainSummary, m = 3)
head(train130)
length(train130[[1]])

#------------------------------------------------------------------#
#----------------------Location Estimation-------------------------#
#------------------------------------------------------------------#

#Build function findNN that finds the k nearest neighbors from any point
#Paramters are a numeric vector of 6 new signal strengths and the return value from selectTrain()
#Function returns the locations of the training observations in order of closeness to the new observation's signal strength
findNN <- function(newSignal, trainSubset) {
  diffs = apply(trainSubset[ , 4:9], 1, function(x) x - newSignal)
  dists = apply(diffs, 2, function(x) sqrt(sum(x^2)))
  closest = order(dists)
  return(trainSubset[closest, 1:3])
}

#Distance Weighted average the first k locations returned by findNN to estimate location of new observation
predXY = function(newSignals, newAngles, trainData, 
                  numAngles = 1, k = 3){
  
  closeXY = list(length = nrow(newSignals))
  
  for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    closeXY[[i]] = 
      findNN(newSignal = as.numeric(newSignals[i, ]), trainSS)
  }
  
  estXY = lapply(closeXY, 
                 function(x) sapply(x[ , 2:3], 
                                    function(x) mean(x[1:k])))
  estXY = do.call("rbind", estXY)
  return(estXY)
}

#Test function with 3 nearest neighbors and 3 orientations
estXYk3 <- predXY(newSignals = testSummary[ , 6:11], 
                 newAngles = testSummary[ , 4], 
                 trainSummary, numAngles = 3, k = 3)
#Test function with 1 nearest neighbors and 3 orientations
estXYk1 <- predXY(newSignals = testSummary[ , 6:11], 
                  newAngles = testSummary[ , 4], 
                  trainSummary, numAngles = 3, k = 1)

#Map of actual vs predicted locations of new signals shows model accuracy
#K=3 more accurate than k=1 because errors are shorter distance and less problematic because they follow the hallways
#Check website for plotting code

#Calculate sum of squared errors to measure fit
calcError <- function(estXY, actualXY) {
  sum( rowSums( (estXY - actualXY)^2) )
}
actualXY <- testSummary[ , c("posX", "posY")]
sapply(list(estXYk1, estXYk3), calcError, actualXY)
#SSE lower for k=3 model, reaffirming that it's better than k=1

#Use k-fold cross validation to find optimal k for KNN

#Perform k-fold cross validation with k=11 so each fold has 15 locations randomly selected
v <- 11
permuteLocs <- sample(unique(trainSummary$posXY))
permuteLocs <- matrix(permuteLocs, ncol = v, 
                      nrow = floor(length(permuteLocs)/v))

testFold <- subset(trainSummary, posXY %in% permuteLocs[ , 1])

reshapeSS <- function(data, varSignal = "signal", 
                      keepVars = c("posXY", "posX","posY"),
                      sampleAngle = FALSE, 
                      refs = seq(0, 315, by = 45)) {
  byLocation =
    with(data, by(data, list(posXY), 
                  function(x) {
                    if (sampleAngle) {
                      x = x[x$angle == sample(refs, size = 1), ]}
                    ans = x[1, keepVars]
                    avgSS = tapply(x[ , varSignal ], x$mac, mean)
                    y = matrix(avgSS, nrow = 1, ncol = 6,
                               dimnames = list(ans$posXY,
                                               names(avgSS)))
                    cbind(ans, y)
                  }))
  
  newDataSS = do.call("rbind", byLocation)
  return(newDataSS)
}

train <- train[ train$mac != "00:0f:a3:39:dd:cd", ]

keepVars <- c("posXY", "posX","posY", "orientation", "angle")

testCVSummary <- reshapeSS(train, keepVars = keepVars, 
                          sampleAngle = TRUE)

testFold <- subset(testCVSummary, 
                  posXY %in% permuteLocs[ , 1])

trainFold <- subset(trainSummary,
                   posXY %in% permuteLocs[ , -1])

estFold <- predXY(newSignals = testFold[ , 6:11], 
                 newAngles = testFold[ , 4], 
                 trainFold, numAngles = 3, k = 3)

actualFold <- testFold[ , c("posX", "posY")]
calcError(estFold, actualFold)

#Wrap the code above into loops over the folds and number of neighbors for K=20
K <- 20
err <- rep(0, K)

for (j in 1:v) {
  testFold = subset(testCVSummary, 
                    posXY %in% permuteLocs[ , j])
  trainFold = subset(trainSummary,
                     posXY %in% permuteLocs[ , -j])
  actualFold = testFold[ , c("posX", "posY")]
  
  for (k in 1:K) {
    estFold = predXY(newSignals = testFold[ , 6:11],
                     newAngles = testFold[ , 4], 
                     trainFold, numAngles = 3, k = k)
    err[k] = err[k] + calcError(estFold, actualFold)
  }
}

pdf(file = "Geo_CVChoiceOfK.pdf", width = 10, height = 6)
oldPar = par(mar = c(4, 3, 1, 1))
plot(y = err, x = (1:K),  type = "l", lwd= 2,
     ylim = c(1200, 2100),
     xlab = "Number of Neighbors",
     ylab = "Sum of Square Errors")

rmseMin <- min(err)
kMin <- which(err == rmseMin)[1]
segments(x0 = 0, x1 = kMin, y0 = rmseMin, col = gray(0.4), 
         lty = 2, lwd = 2)
segments(x0 = kMin, x1 = kMin, y0 = 1100,  y1 = rmseMin, 
         col = grey(0.4), lty = 2, lwd = 2)

#mtext(kMin, side = 1, line = 1, at = kMin, col = grey(0.4))
text(x = kMin - 2, y = rmseMin + 40, 
     label = as.character(round(rmseMin)), col = grey(0.4))
par(oldPar)
dev.off()

#Since errors level out at k=5 and higher, use k= and apply it to training and test data
estXYk5 <- predXY(newSignals = testSummary[ , 6:11], 
                 newAngles = testSummary[ , 4], 
                 trainSummary, numAngles = 3, k = 5)

#Tally the errors
calcError(estXYk5, actualXY)

#Below code is more optimized
predXY <- function(newSignals, newAngles, trainData, 
                  numAngles = 1, k = 3){
  
  closeXY = list(length = nrow(newSignals))
  
  for (i in 1:nrow(newSignals)) {
    trainSS = selectTrain(newAngles[i], trainData, m = numAngles)
    closeXY[[i]] = findNN(newSignal = as.numeric(newSignals[i, ]),
                          trainSS)
  }
  
  estXY = lapply(closeXY, function(x)
    sapply(x[ , 2:3], 
           function(x) mean(x[1:k])))
  estXY = do.call("rbind", estXY)
  return(estXY)
}
#Notice that number of angles is not optimized through cross validation, but it could be for a better model

#------------------------------------------------------------------#
#------------------------Drawing Results---------------------------#
#------------------------------------------------------------------#

floorErrorMap <- function(estXY, actualXY, trainPoints = NULL, AP = NULL){
  
  plot(0, 0, xlim = c(0, 35), ylim = c(-3, 15), type = "n",
       xlab = "", ylab = "", axes = FALSE)
  box()
  if ( !is.null(AP) ) points(AP, pch = 15)
  if ( !is.null(trainPoints) )
    points(trainPoints, pch = 19, col="grey", cex = 0.6)
  
  points(x = actualXY[, 1], y = actualXY[, 2], 
         pch = 19, cex = 0.8 )
  points(x = estXY[, 1], y = estXY[, 2], 
         pch = 8, cex = 0.8 )
  segments(x0 = estXY[, 1], y0 = estXY[, 2],
           x1 = actualXY[, 1], y1 = actualXY[ , 2],
           lwd = 2, col = "red")
}

trainPoints <- trainSummary[ trainSummary$angle == 0 & 
                              trainSummary$mac == "00:0f:a3:39:e1:c0" ,
                            c("posX", "posY")]

pdf(file="GEO_FloorPlanK3Errors.pdf", width = 10, height = 7)
oldPar <- par(mar = c(1, 1, 1, 1))
floorErrorMap(estXYk3, testSummary[ , c("posX","posY")], 
              trainPoints = trainPoints, AP = ap)
par(oldPar)
dev.off()

pdf(file="GEO_FloorPlanK5Errors.pdf", width = 10, height = 7)
oldPar <- par(mar = c(1, 1, 1, 1))
floorErrorMap(estXYk5, testSummary[ , c("posX","posY")], 
              trainPoints = trainPoints, AP = ap)
par(oldPar)
dev.off()