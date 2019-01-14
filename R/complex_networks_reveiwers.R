# complex_networks_reviewers.R
# further experiments and analysis conducted as a result of reviewer feedback.

#load("ComplexNets10thApril2018.RData")
source("complex_networks_functions.R")  # load in the functions required for this work. 
library(AUC)
library(randomForest)


## POINT 1: DATA IMBALANCE PROBLEM, WE HAVE 1,449 TARGETS AND 11,567 NON-TARGETS
# https://stats.stackexchange.com/questions/157714/r-package-for-weighted-random-forest-classwt-option/158030#158030
# https://stats.stackexchange.com/questions/163251/creating-a-test-set-with-imbalanced-data/163567#163567
#setwd("C:/R-files/NeuralNet")  
#load("NCA-27thMarch2018.RData")
## Convert matrix to dataframe and balance out the data by undersampling
mnew <- data.table::transpose(as.data.frame(mmt))
mnew <- data.frame(mnew)
colnames(mnew) <- rownames(mmt)
rownames(mnew) <- colnames(mmt)
positives <- mnew[mnew$targets == 1,]  # get all targets (1,449)
negatives <- mnew[mnew$targets == 0,]  # get all nontargets (11,567)

## Prepare a training and a test set / as per reviewer guidelines
nmin <- round(nrow(positives)/2)
mindex <- sample(nrow(positives),nmin) # indices of minority training samples
minTRAIN <- data.frame(positives[mindex,])  # populate the minority train dataframe
minTEST <- data.frame(positives[-mindex,])

nmaj <- round(nrow(negatives)/2)
majindex <- sample(nrow(negatives),nmaj) # indices of majority training samples
majTRAIN <- data.frame(negatives[majindex,]) # populate the minority train dataframe
majTEST <-  data.frame(negatives[-majindex,])

exindex1 <- sample(nrow(majTRAIN),500) # indices of maj
exindex2 <- sample(nrow(majTEST),500) # indices of maj

majEXPLORE1 <- data.frame(majTRAIN[exindex1,]) # EXPLORE data will be feed into RF
majEXPLORE2 <- data.frame(majTEST[exindex2,]) # EXPLORE data will be feed into RF
majEXPLORE <- rbind(majEXPLORE1[,1:149],majEXPLORE2[,1:149])  #

majTEST <- data.frame(majTEST[-exindex2,])
majTRAIN <-data.frame(majTRAIN[-exindex1,])

allTRAIN <- rbind(minTRAIN[,1:149],majTRAIN[,1:149])  # data for training
allTEST  <- rbind(minTEST[,1:149],majTEST[,1:149])   # data for testing
yTEST <- as.factor(c(minTEST[,150],majTEST[,150]))  # class labels testing
yTRAIN <- as.factor(c(minTRAIN[,150],majTRAIN[,150])) # class labels training

## train default RF and then with 2x 5x and 10x upsampling by stratification
rf1 <- randomForest(yTRAIN~.,allTRAIN, mtry=120, ntree=5000,nodesize=1)
rf2 <- randomForest(yTRAIN~.,allTRAIN, mtry=120, ntree=5000,nodesize=1,sampsize=c(100,200),strata=yTRAIN)
rf3 <- randomForest(yTRAIN~.,allTRAIN, mtry=120, ntree=5000,nodesize=1,sampsize=c(100,500),strata=yTRAIN)
rf4 <- randomForest(yTRAIN~.,allTRAIN, mtry=120, ntree=5000,nodesize=1,sampsize=c(300,724),strata=yTRAIN)
#rf4 <- randomForest(yTRAIN~.,allTRAIN, mtry=120, ntree=5000,nodesize=1,sampsize=c(724,724),strata=yTRAIN)

## plot ROC for training data votes
par(mfrow=c(1,1))
plot(roc(rf1$votes[,2],factor(1 * (rf1$y==1))),main="ROC curves for four models predicting class 1")
plot(roc(rf2$votes[,2],factor(1 * (rf2$y==1))),col=2,add=TRUE)
plot(roc(rf3$votes[,2],factor(1 * (rf3$y==1))),col=3,add=TRUE)
plot(roc(rf4$votes[,2],factor(1 * (rf4$y==1))),col=4,add=TRUE)

for (i in 1:4){
if(i==1) rf_model <- rf1
if(i==2) rf_model <- rf2
if(i==3) rf_model <- rf3
if(i==4) rf_model <- rf4

predicted_rf_train <- predict(rf_model,newdata=allTRAIN,type = "prob")
pred_rf_train <- ROCR::prediction((predicted_rf_train[,2]),yTRAIN)
rf.pr.train <- ROCR::performance(pred_rf_train, "prec", "rec")
#plot(rf.pr.train)

## Ok, so now test out RF on test data
predicted_rf_test <- predict(rf_model,allTEST) # cant use "prob" as confusionmatrix cant use it(error)
confusionMatrix(data=predicted_rf_test,reference=yTEST,positive="1") 
predicted_rf_test <- predict(rf_model,allTEST,type="prob")  # ROCR functions need type="prob"
pred_rf_test <- ROCR::prediction((predicted_rf_test[,2]),yTEST)
rf.roc.test <- ROCR::performance(pred_rf_test, "tpr", "fpr")
rf.pr.test <- ROCR::performance(pred_rf_test, "prec", "rec")

if(i==1){rf1.roc <- rf.roc.test; rf1.pr <- rf.pr.test}
if(i==2){rf2.roc <- rf.roc.test; rf2.pr <- rf.pr.test}
if(i==3){rf3.roc <- rf.roc.test; rf3.pr <- rf.pr.test}
if(i==4){rf4.roc <- rf.roc.test; rf4.pr <- rf.pr.test}
}

#rf.auc <- performance(pred_rf_test, measure = "auc")
#cat("\nAUC=",as.character(rf.auc@y.values[1]))
#plot(rf.pr.test)
#plot(rf.roc.test)

# plot ROC and PR using the prettier ggplot2 graphs
attributes(rf1.roc)$roc_name <- "No sampling"
attributes(rf2.roc)$roc_name <- "100:200"
attributes(rf3.roc)$roc_name <- "100:500"
attributes(rf4.roc)$roc_name <- "300:724"
roc_plot(rf1.roc,rf2.roc,rf3.roc,rf4.roc)

# PR plots of classifiers
attributes(rf1.pr)$pr_name <- "No sampling"
attributes(rf2.pr)$pr_name <- "100:200"
attributes(rf3.pr)$pr_name <- "100:500"
attributes(rf4.pr)$pr_name <- "300:724"
pr_plot(rf1.pr,rf2.pr,rf3.pr,rf4.pr)


# Now feed in majEXPLORE data, i.e. the 1,000 proteins reserved for detecting targets
candidates_rf <-  predict(rf_model,majEXPLORE,type="prob")
#gs <- get_gstatistics(ppi_net)
my_table <-  make_table(candidates_rf,gs[[2]],0.80)   # 0.5 = 178 targets; 0.8=36 targets
my_table$prob <- as.numeric(my_table$prob)  # convert from strings to numbers
xtable(my_table,digits=c(0,0,0,2,2,0))


## POINT 2: WHAT IS THE TRUE INFLUENCE OF K-CORENESS?
# build data structure to train Random Forest on coreness
kcore <- graph.coreness(as.undirected(ppi_net))
kcore <- as.data.frame(kcore)
kcore$proteins <- rownames(kcore)
rownames(coreness) <- c()
kcore$target <- vertex_attr(ppi_net,"target")
kcore$target <- as.factor(kcore$target)

rf5 <- randomForest(y=kcore$target, x=as.data.frame(kcore$kcore), 
                    ntree=10000,nodesize=1,type="classification",
                    sampsize=c(300,700),strata=kcore$target)
rf_model <- rf5
predicted_rf_kcore <- predict(rf_model,newdata=as.data.frame(kcore$kcore),type = "prob")
pred_rf_kcore <- ROCR::prediction((predicted_rf_kcore[,2]),kcore$target)
rf.pr.kcore <- ROCR::performance(pred_rf_kcore, "prec", "rec")
plot(rf.pr.kcore)
rf.roc.kcore <- ROCR::performance(pred_rf_kcore,"tpr", "fpr")
plot(rf.roc.kcore)

predicted_rf_kcore <- predict(rf_model,newdata=as.data.frame(kcore$kcore)) # cant use "prob" as confusionmatrix cant use it(error)
confusionMatrix(data=predicted_rf_kcore,reference=kcore$target,positive="1") 

boxplot(kcore~target,data=kcore, main="K-core Data",xlab="Targetness", ylab="k-coreness") 

#######################################################################################
  
## Prepare a training and a test set 
ntrain <- round(nrow(bdata)*0.8) # number of training examples
tindex <- sample(nrow(bdata),ntrain) # indices of training samples
xtrain <- data.frame(bdata[tindex,])
xtest <-  data.frame(bdata[-tindex,])
ytrain <- xtrain[,150]  # class labels for training data (column 150 is class label)
ytest <- xtest[,150]   # class labels for test data
ytest <- as.factor(ytest)
 

# POINT 8 - retrain RandomForest using different numbers of trees
# various numbers of trees were used, 1500, 1000, 600, 500, 100, 50, 10, 5 and 1
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
plot(rf.pr.train)

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
#y[1]<- 1 # sorry fudge to get rid of NaN
df1 <- data.frame(x,y)
x <- rf.pr.train@x.values
y <- rf.pr.train@y.values
y <- as.numeric(unlist(y[[1]]))
x <- as.numeric(unlist(x[[1]]))
#y[1] <-1 # sorry fudge to get rid of NaN
df2 <- data.frame(x,y)

ggplot(df1,aes(x,y))+
  geom_line(aes(color="Test data PR"),size=1)+
  geom_line(data=df2,aes(color="Train data PR"),size=1)+
  labs(x="Precision",y="Recall") +
  labs(color="Legend") +
  theme(legend.position = c(0.8, 0.2))

# ########## select proteins that are nontargets but not in train or test set #########
# use the trained classifiers on new data "unknown" for potential targets

allnontargets <- mcrap[mcrap$targets == 0,]
unknown <- negatives[!rownames(allnontargets) %in% rownames(negatives),]
unknown <- data.frame(unknown)

# shape up new data
uindex <- base::sample(nrow(unknown),10) # indices of training samples
candidates <- unknown[sample(1:nrow(unknown), 5000,replace=FALSE),] 
candidates <- unknown[1:2500,1:149]
candidates <- data.frame(candidates)
candidates_rf <-  predict(rf_model,candidates,type="prob")
table(candidates_rf)
candidates_rf <- data.frame(candidate = rownames(candidates), target=as.vector(candidates_rf),stringsAsFactors = FALSE)

candidates_rf <- filter(candidates_rf,target==1)

# get top ranking unseen proteins that may be candidates as protein targets
targettype <- predict(rf.model, candidates, type="prob")

proteins <- rownames(targettype)
targettype <- data.frame(targettype)
target <- as.vector(targettype$X1)
targettype <- data.frame(protein=proteins, prob=targettype$X1)
targettype <- filter(targettype, prob >= 0.8) # probabilities equal to greater than 0.8 classed as target

# ADD protein type to dataframe: protein_class for protein types
targettype[] <- lapply(targettype, as.character)
rlen <- nrow(targettype)
tablestuff <- data.frame(TargetClass="FakeType", Gene="RU12",stringsAsFactors = FALSE)
for (i in 1:rlen){
  tempstuff <- data.frame(filter(protein_class,Gene == targettype$protein[i]))
  cat("\nData is...",tempstuff$TargetClass)
  if(length(tempstuff$TargetClass)==0){
    tempstuff <- data.frame(TargetClass="Unknown",Gene=targettype$protein[i],stringsAsFactors = FALSE)}
  tablestuff <- rbind(tablestuff,tempstuff) 
}
# last entry is just to instantiate datastructure, so delete
tablestuff <- tablestuff[-1,] 

# ADD hubs to dataframe: hlist for protein names that are hubs
myhubs <- vector(mode="character", length=rlen)    
for (i in 1:rlen){
  temphubs <- hlist[grep(targettype$protein[i],hlist)]
  cat("\nFound a hub is...",temphubs)
  if(length(temphubs)==0){
    myhubs[i] <- "non-hub"}else{
      myhubs[i] <- "hub"}
}

tablestuff <- cbind(tablestuff,myhubs)
prob <- targettype$prob
tablestuff <- cbind(tablestuff,prob)

# ADD k-coreness to dataframe:
coreness <- graph.coreness(as.undirected(ppi_net))
Genes <- names(coreness)
Core <- as.vector(coreness)
coreness <- data.frame(coreness=Core, Gene=Genes, stringsAsFactors = FALSE)

core <- vector(mode="integer", length=rlen)    
for (i in 1:rlen){
  tempcore <- coreness[grep(tablestuff$Gene[i],coreness$Gene),]
  cat("\nCoreness..",tempcore$coreness)
  if(length(tempcore$coreness)==0){
    core[i] <- 0}else{
      core[i] <- tempcore$coreness}
}

tablestuff <- cbind(tablestuff,core)
rownames(tablestuff) <- NULL
tablestuff %>% mutate_if(is.factor, as.character) -> tablestuff  # convert nasty factors to strings
tablestuff <- arrange(tablestuff,desc(core))
xtable(tablestuff)

# POINT 6 - is k-coreness or GO-Slim the main factor in target prediction??
# so we have 2,500 unknown proteins, 
group_candidates <-  predict(rf_model,candidates,type="prob")
proteins <- rownames(group_candidates)
group_candidates <- data.frame(group_candidates)
target <- as.vector(group_candidates$X1)
group_candidates <- data.frame(protein=proteins, prob=group_candidates$X1)
# ADD protein type to dataframe: protein_class for protein types
group_candidates[] <- lapply(group_candidates, as.character)
rlen <- nrow(group_candidates)
tablestuff1 <- data.frame(TargetClass="FakeType", Gene="RU12",stringsAsFactors = FALSE)
for (i in 1:rlen){
  tempstuff <- data.frame(filter(protein_class,Gene == group_candidates$protein[i]))
  cat("\nData is...",tempstuff$TargetClass)
  if(length(tempstuff$TargetClass)==0){
    tempstuff <- data.frame(TargetClass="Unknown",Gene=targettype$protein[i],stringsAsFactors = FALSE)}
  tablestuff1 <- rbind(tablestuff1,tempstuff) 
}
# last entry is just to instantiate datastructure, so delete
tablestuff1 <- tablestuff1[-1,] 
# ADD hubs to dataframe: hlist for protein names that are hubs
myhubs <- vector(mode="character", length=rlen)    
for (i in 1:rlen){
  temphubs <- hlist[grep(group_candidates$protein[i],hlist)]
  cat("\nFound a hub is...",temphubs)
  if(length(temphubs)==0){
    myhubs[i] <- "non-hub"}else{
      myhubs[i] <- "hub"}
}

tablestuff1 <- cbind(tablestuff1,myhubs)
prob <- group_candidates$prob
tablestuff1 <- cbind(tablestuff1,prob)

# ADD k-coreness to dataframe:
coreness <- graph.coreness(as.undirected(ppi_net))
Genes <- names(coreness)
Core <- as.vector(coreness)
coreness <- data.frame(coreness=Core, Gene=Genes, stringsAsFactors = FALSE)

core <- vector(mode="integer", length=rlen)    
for (i in 1:rlen){
  tempcore <- coreness[grep(tablestuff1$Gene[i],coreness$Gene),]
  cat("\nCoreness..",tempcore$coreness)
  if(length(tempcore$coreness)==0){
    core[i] <- 0}else{
      core[i] <- tempcore$coreness}
}

tablestuff1 <- cbind(tablestuff1,core)
rownames(tablestuff1) <- NULL
tablestuff1 %>% mutate_if(is.factor, as.character) -> tablestuff1  # convert nasty factors to strings
tablestuff1 <- arrange(tablestuff1,desc(core))

tablestuff1 <- na.omit(tablestuff1)
not_targets <- filter(tablestuff1,prob < 0.5)
are_targets <- filter(tablestuff1,prob >=0.5)
are_targets %>% summarise(ave=median(core))
not_targets %>% summarise(ave=median(core))
are_targets %>% summarise(vari=IQR(core))
not_targets %>% summarise(vari=IQR(core))
# proves more or less that k-coreness is not discriminating enough to divide targets and nontargets


# POINT 3: will using the entire go vocab confer any accuracy benefits compared with go-slim?
all_go <- annotate_with_go(ppi_net)

# Now get all the GO terms for each protein, 
plist <- "NULL"
for (i in 1:length(all_go)){
  ptemp <- unlist(all_go[[i]])
  plist <- c(plist,ptemp)
}
plist <- plist[-1]  # remove rubbish 1st entry
plen <- length(unique(plist))  # how unique GO terms do we have?
plist <- unique(plist) # overwrite

# create the matrix for classification algorithms
mma <- matrix(0, plen, length(names(all_go)))  # Number of unique GO terms x Number of genes
colnames(mma) <- names(all_go)  # each colname is a protein
rownames(mma) <- plist          # each rowname is a GO term


# Load the huge mma matrix data, it takes 2.5 hours to calculate from raw data on my poor laptop! 
load("mma.RData")

mmt <- data.table::transpose(as.data.frame(mma))
mmt <- data.frame(mmt)
colnames(mmt) <- rownames(mma)
rownames(mmt) <- colnames(mma)

positives <- mmt[mmt$Target == 1,]  # get all targets (1,443)
negatives <- mmt[mmt$Target == 0,]  # get all nontargets (11,554)
allnegatives <- data.frame(negatives)
negatives <- sample_n(negatives, nrow(positives)) # only use 1,443 of them to match positives
balanced_dat <- rbind(positives,negatives)

## Prepare a training and a test set ##
ntrain <- round(nrow(balanced_dat)*0.8) # number of training examples
tindex <- sample(nrow(balanced_dat),ntrain) # indices of training samples
xtrain <- data.frame(balanced_dat[tindex,])
xtest <-  data.frame(balanced_dat[-tindex,])

ytrain <- xtrain[,ncol(xtrain)]  # class labels for training data
ytest <- xtest[,ncol(xtest)]   # class labels for test data
xtest <- xtest[,1:ncol(xtest)-1]  # miss out class label and update structure
ytest <- as.factor(ytest); #

# train RF on full GO terms :using 500 trees takes 1.5 hours with 56% test set accuracy,sensitivity=0.16,specificity=0.97
# train RF on full GO terms :using 1000 trees takes 3.5 hours with 57% test set accuracy,sensitivity=0.16,specificity=0.98

rf_model_go <-randomForest(as.factor(ytrain)~.,data=xtrain[,1:ncol(xtrain)-1],
                           proximity=TRUE,keep.forest=FALSE,ntree=1000,mtry=12)

pred_rf <- stats::predict(rf_model_go,xtest)
confusionMatrix(data=pred_rf,reference=ytest,positive="1")


# POINT 4 - problem of under-sampling in the "TESTING SET", 
# https://stats.stackexchange.com/questions/97926/suggestions-for-cost-sensitive-learning-in-a-highly-imbalanced-setting
# https://www.r-bloggers.com/handling-class-imbalance-with-r-and-caret-an-introduction/
# https://topepo.github.io/caret/measuring-performance.html
# https://svds.com/learning-imbalanced-classes/

library(DMwR)   # DMwr to balance the unbalanced class

balanced.data <- SMOTE(Class ~., dresstrain, perc.over = 4800, k = 5, perc.under = 1000)
as.data.frame(table(balanced.data$Class))




#################################################################################################
# POINT 3: will using the entire go vocab confer any accuracy benefits compared with go-slim?
#all_go <- annotate_with_go(ppi_net)
#
# Now get all the GO terms for each protein, 
#plist <- "NULL"
#for (i in 1:length(all_go)){
#  ptemp <- unlist(all_go[[i]])
#  plist <- c(plist,ptemp)
#}
#plist <- plist[-1]  # remove rubbish 1st entry
#plen <- length(unique(plist))  # how unique GO terms do we have?
#plist <- unique(plist) # overwrite
#
# create the matrix for classification algorithms
#mma <- matrix(0, plen, length(names(all_go)))  # Number of unique GO terms x Number of genes
#colnames(mma) <- names(all_go)  # each colname is a protein
#rownames(mma) <- plist          # each rowname is a GO term
#
# POINT 3: This is edited out because data matrix takes some time to compute
# Populate matrix with 1's where GO term(s) are present for that protein: TAKES AT LEAST 2.5 HOURS TO COMPUTE!!
#for(i in 1:ncol(mma)){    # for each protein (column) annotate matrix with GO terms allocated to it. 
#  ptemp <- all_go[[i]]
#  rowkeep <- which(rownames(mma) %in% ptemp) # recall GO terms are rownames,
#  mma[rowkeep,i] <- 1
#}
#
# Now add the target status to "mma" as the training label
#target <- matrix(0, ncol(mma), 1)
#colkeep <- which(colnames(mma) %in% drug_targets$Gene)
#target[colkeep] <- 1
#mmatemp <- rbind(mma,t(target))
#mma <- mmatemp
#rownames(mma)[nrow(mma)]<-"Target"
#save(mma,file="mma.RData")
###################################################################################################
