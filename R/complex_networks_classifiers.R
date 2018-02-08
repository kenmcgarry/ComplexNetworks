# complex_networks_classifiers.R
# Create SVM and RandomForest based classifiers
# ROC and PR curves from ROCR package
# https://www.r-bloggers.com/machine-learning-using-support-vector-machines/

library(ROCR)
library(kernlab)
library("e1071")
library(caret)

#  matrix mmt [150 x 12997] contains the raw data, needs to split. Row 150 contains "class" labels.
mmi <- mmt[,1:6000]
mmi <- data.frame(mmi)

mtest <- mmt[,6001:12000]
mtest <- data.frame(mtest)
#mmi$targets <- as.factor(mmi$targets)

for (i in 1:length(targettype)){
  if(targettype[i] == 1) {  
     targettype[i] <- runif(1, 7.0, 9.9)}}
for (i in 1:length(targettype)){
  if(targettype[i] == 0) {  
    targettype[i] <- runif(1,0.01,0.02)}}

# http://amunategui.github.io/smote/  # Supersampling Rare Events - build up smaller target class
# https://topepo.github.io/caret/index.html
# http://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
# fit a random forest model (using ranger package), re-runs the model over 25 bootstrap samples and 
# across 3 options of the tuning parameter
targets <- (mmi[150,]); targets <- t(targets); rownames(targets) <- NULL; colnames(targets) <- NULL
targets <- as.factor(targets)
testtargets <- mtest[150,];testtargets <- t(testtargets);rownames(testtargets) <- NULL; colnames(testtargets) <- NULL 
testtargets <- as.factor(testtargets)

rf_fit <- train(as.factor(targets) ~., data = t(mmi),method = "ranger")
targettype <- predict(rf_fit,t(mmi))

# https://www.r-bloggers.com/a-small-introduction-to-the-rocr-package/
pred <- prediction(targettype,mmi$targets)
roc1.perf <- performance(pred, "tpr", "fpr")
plot(roc1.perf)
roc3.perf <- performance(pred, "prec", "rec")
plot(roc3.perf)
roc4.perf <- performance(pred, "spec", "sens")
plot(roc4.perf)

# accuracy on training data
pred <- prediction(as.numeric(targettype),targets)
table(targettype,targets)
roc1.perf <- performance(pred, "tpr", "fpr")
plot(roc1.perf)

# accuracy on test data
targettype2 <- predict(rf_fit,t(mtest))
pred2 <- prediction(as.numeric(targettype2),testtargets)
table(targettype2,targets)
roc2.perf <- performance(pred2, "tpr", "fpr")
plot(roc2.perf)

# poor accuracy might be solved by using fewer non-target exemplars
mcrap <- t(mmt)
mcrap <- as.data.frame(mcrap)
positives <- mcrap[mcrap$targets == 1,]  # get all targets (1,443)
negatives <- mcrap[mcrap$targets == 0,]  # get all targets (11,554)
negatives <- sample_n(negatives, nrow(positives)) # only use 1,443 of them to match positives
balanced_dat <- rbind(positives,negatives)

## Prepare a training and a test set ##
ntrain <- round(nrow(balanced_dat)*0.8) # number of training examples
tindex <- sample(nrow(balanced_dat),ntrain) # indices of training samples
xtrain <- data.frame(balanced_dat[tindex,])
xtest <- data.frame(balanced_dat[-tindex,])

ytrain <- xtrain[,150]  # class labels for training data
ytest <- xtest[,150]   # class labels for test data
ytest <- as.factor(ytest); #rownames(ytest) <- NULL; colnames(ytest) <- NULL
#xtest <- t(xtest)
xtest <- data.frame(xtest)

# Ok, so retain Random Forest
rf_fit <- train(as.factor(ytest) ~., data=(xtest), method = "ranger")
#targettype <- predict(rf_fit,t(xtrain))








