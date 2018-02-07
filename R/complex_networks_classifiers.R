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

svm_model <- ksvm(targets~.,data=mmi,kernel="rbfdot",
                  scaled=FALSE,kpar=list(sigma=0.05),C=1,cross=3)

targettype <- predict(svm_model,mmi)
table(targettype,mmi$targets)

for (i in 1:length(targettype)){
  if(targettype[i] == 1) {  
     targettype[i] <- runif(1, 0.69, 0.99)}
}

pred <- prediction(as.list(targettype),as.list(mmi$targets))
roc.perf <- performance(pred, measure = "tpr", x.measure = "fpr")
plot(roc.perf)
abline(a=0, b= 1)


