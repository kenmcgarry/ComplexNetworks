# complex_networks_classifiers.R
# Create SVM and RandomForest based classifiers
# ROC and PR curves from ROCR package
# https://www.r-bloggers.com/machine-learning-using-support-vector-machines/

library(ROCR)
library(kernlab)
library("e1071")
library(caret)
library(data.table)
library(ranger)

#  matrix mmt [150 x 12997] contains the raw data, needs to split. Row 150 contains "class" labels.
mmi <- mmt[,1:1000]
mmi <- data.frame(mmi)

mtest <- mmt[,1001:2000]
mtest <- data.frame(mtest)


# http://amunategui.github.io/smote/  # Supersampling Rare Events - build up smaller target class
# https://topepo.github.io/caret/index.html
# http://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
# fit a random forest model (using ranger package), re-runs the model over 25 bootstrap samples and 
# across 3 options of the tuning parameter
targets <- (mmi[150,]); targets <- t(targets); rownames(targets) <- NULL; colnames(targets) <- NULL
targets <- as.factor(targets)
testtargets <- mtest[150,];testtargets <- t(testtargets);rownames(testtargets) <- NULL; colnames(testtargets) <- NULL 
testtargets <- as.factor(testtargets)

rf_fit <- train(as.factor(targets) ~., data = t(mmi[1:149,]),method = "ranger")
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

#===============
# poor accuracy might be solved by using fewer non-target exemplars
mcrap <- data.table::transpose(as.data.frame(mmt))
mcrap <- data.frame(mcrap)
colnames(mcrap) <- rownames(mmt)
rownames(mcrap) <- colnames(mmt)

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

xxtest <- data.table::transpose(xtest)
colnames(xxtest) <- rownames(xtest)
rownames(xxtest) <- colnames(xtest)

#############################################################

# Ok, so retain Random Forest
rf_fit <- train(as.factor(targets) ~., data=xtrain, method = "ranger",importance = "impurity_corrected")# 
vimport <- rf_fit[[11]]
vimport <- vimport$variable.importance
vimport <- data.frame(GOterm=names(vimport),importance=unname(vimport),stringsAsFactors = FALSE)
vimport <- vimport[order(vimport$importance,decreasing=TRUE),]
head(vimport,10)

targettype <- predict(rf_fit,xtrain)
targettype <- factor2int((targettype))
pred       <- prediction(targettype,xtrain$targets)
roc1.perf  <- performance(pred, "tpr", "fpr")
plot(roc1.perf,col="red", lwd=5)

targettype <- predict(rf_fit,xtest)
acc <- table(xtest$targets, targettype)
acc <- as.vector(acc); TN <- acc[1]; FN <- acc[2]; FP <- acc[3]; TP <- acc[4]  
cat("\naccuracy calculated by (TP + TN)/(TP + TN + FP + FN) = ",(TP + TN)/(TP + TN + FP + FN))
targettype <- factor2int(targettype)
pred       <- prediction(targettype,xtest$targets)
roc2.perf  <- performance(pred, "tpr", "fpr")
plot(roc2.perf,add = TRUE, col="blue",lwd=5)

perf <- performance(pred, measure = "auc")
cat("\nAUC: ", as.numeric(perf@y.values)) 
perf <- performance(pred, measure = "acc")
cat("\nAccuracy: ",max(perf@x.values[[1]]) )

ind = which.max(slot(perf,"x.values")[[1]] )
acc = slot(perf, "y.values")[[1]][ind]
cutoff = slot(perf, "x.values")[[1]][ind]
print(c(accuracy= acc, cutoff = cutoff))


#pred1 <- prediction(Avec.pred1, Avec)
#perf1 <- performance(pred1, "tpr", "fpr")
#plot(perfall, col="red", lwd=5)
#plot(perf1, add = TRUE, col="blue",lwd=5)
#abline(a=0, b= 1,lty=2 )
# CHUNK 33
#perf1.auc <- performance(pred1, "auc")
#slot(perf1.auc, "y.values")

