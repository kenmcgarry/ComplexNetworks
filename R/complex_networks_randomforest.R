# randomforest package

library(randomForest)

rf.model <- randomForest(as.factor(ytrain) ~., data=xtrain[,1:149],importance=TRUE,proximity = TRUE,keep.forest=TRUE)
print(rf.model)
round(importance(rf.model), 2)
varImpPlot(rf.model,main="",type=2,color="black",pch=16) 

# get top ranking unseen proteins that may be candidates as protein targets
targettype <- predict(rf.model, unknown, type="prob")
proteins <- rownames(targettype)
targettype <- data.frame(targettype)
target <- as.vector(targettype$X1)
targettype <- data.frame(protein=proteins, prob=targettype$X1)
targettype <- filter(targettype, prob > 0.5)  # greater than 0.5 will be classed as target

targettype <- targettype[order(targettype$prob,decreasing = TRUE),] 
head(targettype,10)


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
# Combine prob + core into overall score for candidate protein ranking

# ADD evidence to dataframe: space for papers in literature supporting potential


xtable(tablestuff)


# A bit more on variable importance - must stick this code in appropriate R file!!!
# https://www.r-bloggers.com/variable-importance-plot-and-variable-selection/
# VI_F <- importance(rf.model)
library(clusterGeneration)
library(mnormt)
library(corrplot)

VIP <- rf.model$importance ; VIP <- data.frame(VIP) ;VIP <- VIP[order(VIP[4],decreasing = TRUE),]
names <-rownames(VIP)
names <- names[1:15]
par(mfrow=c(2, 3), xpd=NA)
for (name in names)
  partialPlot(rf.model, xtrain[,1:149], eval(name), main=name, xlab=name,which.class=1)


#plot(boruta_output, cex.axis=.7, las=2, xlab="", main="") 
# plot regression model on probabilities






