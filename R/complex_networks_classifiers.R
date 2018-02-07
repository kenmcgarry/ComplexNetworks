# complex_networks_classifiers.R
# Create SVM and RandomForest based classifiers
# ROC and PR curves from ROCR package
# https://www.r-bloggers.com/machine-learning-using-support-vector-machines/

library(ROCR)
library(kernlab)
library("e1071")
library(caret)

#  matrix mmt [150 x 12997] contains the raw data, needs to split. Row 150 contains "class" labels.
mmi <- t(mmt)
mmi <- data.frame(mmi)
mmi$targets <- as.factor(mmi$targets)

## Prepare a training and a test set ##
ntrain <- round(n*0.8) # number of training examples
tindex <- sample(n,ntrain) # indices of training samples
xtrain <- x[tindex,]
xtest <- x[-tindex,]
ytrain <- y[tindex]

#svm_model <- ksvm(targets~.,data=mmi,kernel="rbfdot",scaled=FALSE,kpar=list(sigma=0.05),C=1,cross=5)
svm_model <- ksvm(targets~.,data=mmi,type="C-svc",kernel="vanilladot",C=100,cross=5,scaled=c()) #cross=5

targettype <- predict(svm_model,mmi)
table(targettype,mmi$targets)

targettype <- as.numeric(targettype)  # change for prediction()

for (i in 1:length(targettype)){
  if(targettype[i] == 1) {  
     targettype[i] <- 0}}
for (i in 1:length(targettype)){
  if(targettype[i] == 2) {  
    targettype[i] <- 1}}


# https://www.r-bloggers.com/a-small-introduction-to-the-rocr-package/
pred <- prediction(targettype,mmi$targets)
roc1.perf <- performance(pred, "tpr", "fpr")
plot(roc1.perf)
roc3.perf <- performance(pred, "prec", "rec")
plot(roc3.perf)
roc4.perf <- performance(pred, "spec", "sens")
plot(roc4.perf)

# http://amunategui.github.io/smote/  # Supersampling Rare Events - build up smaller target class
# https://topepo.github.io/caret/index.html
# http://www.rebeccabarter.com/blog/2017-11-17-caret_tutorial/
# fit a random forest model (using ranger package), re-runs the model over 25 bootstrap samples and 
# across 3 options of the tuning parameter
rf_fit <- train(as.factor(mmi$targets) ~ ., data = mmi,method = "ranger")



