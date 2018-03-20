# complex_networks_reviewers.R
# further experiments and analysis conducted as a result of reviewer feedback.

memory.limit(2010241024*1024) # use more RAM memory (20 GBs)
setwd("C:/R-files/complexnetworks")    # point to where my code lives
load("ComplexNets9thMarch2018.RData")
source("complex_networks_functions.R")  # load in the functions required for this work. 

# POINT 4 - problem of under-sampling in the "TESTING SET", 
# https://stats.stackexchange.com/questions/97926/suggestions-for-cost-sensitive-learning-in-a-highly-imbalanced-setting



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
targettype <- filter(targettype, prob > 0.5)  # greater than 0.5 will be classed as target

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


# Load the huge mma matrix data, it takes 2.5 hours to calculate from raw data! 
load("mma.RData")
# Now add the target status as the training label


# This is edited out because it takes ages to compute
# Populate matrix with 1's where GO term(s) are present for that protein: TAKES AT LEAST 2.5 HOURS TO COMPUTE!!
#for(i in 1:ncol(mma)){    # for each protein (column) annotate matrix with GO terms allocated to it. 
#  ptemp <- all_go[[i]]
#  rowkeep <- which(rownames(mma) %in% ptemp) # recall GO terms are rownames,
#  mma[rowkeep,i] <- 1
#}






