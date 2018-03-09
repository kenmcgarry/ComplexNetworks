# complex_networks_reviewers.R
# further experiments and analysis conducted as a result of reviewer feedback.

memory.limit(2010241024*1024) # use more RAM memory (20 GBs)
setwd("C:/R-files/complexnetworks")    # point to where my code lives
load("ComplexNets8thMarch2018.RData")
source("complex_networks_functions.R")  # load in the functions required for this work. 


# point 8 - retrain RandomForest using different numbers of trees
# Train Random Forest on data 
rf_model <-randomForest(as.factor(ytrain) ~.,data=xtrain[,1:149],proximity=TRUE,keep.forest=TRUE,
                        ntree=500,mtry=12)

predicted_rf_train <- predict(rf_model,newdata=xtrain[,1:149],type = "prob")
predicted_rf_test <- predict(rf_model,newdata=xtest[,1:149],type = "prob")  
pred_rf_test <- ROCR::prediction((predicted_rf_test[,2]),ytest)
rf.roc.test <- ROCR::performance(pred_rf_test, "tpr", "fpr")
rf.pr.test <- ROCR::performance(pred_rf_test, "prec", "rec")

pred_rf_train <- ROCR::prediction((predicted_rf_train[,2]),ytrain)
rf.roc.train <- ROCR::performance(pred_rf_train, "tpr", "fpr")
rf.pr.train <- ROCR::performance(pred_rf_train, "prec", "rec")

pred_rf <- stats::predict(rf_model,xtest)
confusionMatrix(data=pred_rf,reference=ytest,positive="1")
plot(rf.pr.test)

x <- as.numeric(unlist(rf.roc.test@x.values))
y <- as.numeric(unlist(rf.roc.test@y.values))
df1 <- data.frame(x,y)
x <- rf.roc.train@x.values
y <- rf.roc.train@y.values
y <- as.numeric(unlist(y[[1]]))
x <- as.numeric(unlist(x[[1]]))
df2 <- data.frame(x,y)

# ROC curves for train data and test data
ggplot(df1,aes(x,y))+
  geom_line(aes(color="Test data ROC"),size=1)+
  geom_line(data=df2,aes(color="Train data ROC"),size=1)+
  labs(x="False Positive Rate",y="True Positive Rate") +
  labs(color="Legend") +
  theme(legend.position = c(0.8, 0.2))

# PR (precision-recall) curves for train data and test data
x <- as.numeric(unlist(rf.pr.test@x.values))
y <- as.numeric(unlist(rf.pr.test@y.values))
df1 <- data.frame(x,y)
x <- rf.pr.train@x.values
y <- rf.pr.train@y.values
y <- as.numeric(unlist(y[[1]]))
x <- as.numeric(unlist(x[[1]]))
df2 <- data.frame(x,y)

ggplot(df1,aes(x,y))+
  geom_line(aes(color="Test data PR"),size=1)+
  geom_line(data=df2,aes(color="Train data PR"),size=1)+
  labs(x="Precision",y="Recall") +
  labs(color="Legend") +
  theme(legend.position = c(0.8, 0.2))


# point 6 - is k-coreness or GO-Slim the main factor in target prediction??

tablestuff <- tablestuff[order(core,decreasing = FALSE),]



