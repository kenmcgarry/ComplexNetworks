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
negatives <- mcrap[mcrap$targets == 0,]  # get all nontargets (11,554)
allnegatives <- data.frame(negatives)
negatives <- sample_n(negatives, nrow(positives)) # only use 1,443 of them to match positives
balanced_dat <- rbind(positives,negatives)

## Prepare a training and a test set ##
ntrain <- round(nrow(balanced_dat)*0.8) # number of training examples
tindex <- sample(nrow(balanced_dat),ntrain) # indices of training samples
xtrain <- data.frame(balanced_dat[tindex,])
xtest <-  data.frame(balanced_dat[-tindex,])

ytrain <- xtrain[,150]  # class labels for training data
ytest <- xtest[,150]   # class labels for test data
ytest <- as.factor(ytest); #rownames(ytest) <- NULL; colnames(ytest) <- NULL


#############################################################

# Ok, so retain Random Forest
rf_fit <- train(as.factor(targets) ~., data=xtrain, method = "ranger")# 
vimport <- rf_fit[[11]]
vimport <- vimport$variable.importance
vimport <- data.frame(GOterm=names(vimport),importance=unname(vimport),stringsAsFactors = FALSE)
vimport <- vimport[order(vimport$importance,decreasing=TRUE),]
goterms <- head(vimport,15)

targettype <- predict(rf_fit,xtrain)
targettype <- factor2int((targettype))
pred       <- prediction(targettype,xtrain$targets)
roc1.perf  <- performance(pred, "tpr", "fpr")
plot(roc1.perf,col="red", lwd=5)
roc2.perf  <- performance(pred, "sens", "spec")
plot(roc2.perf,col="red", lwd=5)
roc3.perf  <- performance(pred, "prec", "rec")
plot(roc3.perf,col="red", lwd=5)

x <- roc1.perf@x.values
y <- roc1.perf@y.values
y <- as.numeric(as.character(unlist(y[[1]])))
x <- as.numeric(as.character(unlist(x[[1]])))
#qplot(x,y, geom="point",ylab="True Positive Rate",xlab="False Positive Rate")
df1 <- data.frame(x,y)
x <- roc2.perf@x.values
y <- roc2.perf@y.values
y <- as.numeric(as.character(unlist(y[[1]])))
x <- as.numeric(as.character(unlist(x[[1]])))
df2 <- data.frame(x,y)

# Now do ROC graph
ggplot(df1,aes(x,y))+
  geom_line(aes(color="Train data ROC"),size=2)+
  geom_line(data=df2,aes(color="Test data ROC"),size=2)+
  labs(x="False Positive Rate",y="True Positive Rate") +
  labs(color="Legend") +
  theme(legend.position = c(0.8, 0.2))

# now do PR graph
targettype <- predict(rf_fit,xtrain)
targettype <- factor2int((targettype))
pred <- prediction(targettype,xtrain$targets)
roc3.perf  <- performance(pred, "prec", "rec")
x <- roc3.perf@x.values
y <- roc3.perf@y.values
y <- as.numeric(as.character(unlist(y[[1]])))
x <- as.numeric(as.character(unlist(x[[1]])))
df3 <- data.frame(y,x)
df3 <- na.omit(df3)

targettype <- predict(rf_fit,xtest)
targettype <- factor2int(targettype)
pred       <- prediction(targettype,xtest$targets)
roc4.perf  <- performance(pred, "prec", "rec")
x <- roc4.perf@x.values
y <- roc4.perf@y.values
y <- as.numeric(as.character(unlist(y[[1]])))
x <- as.numeric(as.character(unlist(x[[1]])))
df4 <- data.frame(x,y)
df4 <- na.omit(df4)
  
ggplot(df3,aes(x,y))+
  geom_line(aes(color="Train data PR"),size=2)+
  geom_line(data=df4,aes(color="Test data PR"),size=2)+
  labs(x="Precision",y="Recall") +
  labs(color="Legend") +
  theme(legend.position = c(0.8, 0.2))



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

# Examine freq counts of GO terms between targets and non-target proteins on top 15 terms
# truncate name to 10 letters
freqtarget <- filter(balanced_dat, targets ==1)
freqplain <- filter(balanced_dat, targets ==0)
importantGO <- goterms$GOterm; importantGO <- gsub("\\.", ":", importantGO); 
importantGO <- substr(importantGO, start = 1, stop = 10)

tempplain <- freqplain %>% dplyr:: select(starts_with(importantGO[1]))
for (i in 1:length(importantGO)){
  temp <- freqplain %>% dplyr:: select(starts_with(importantGO[i]))
  tempplain <- cbind(temp,tempplain)
}
tempplain <- tempplain[,1:15] # get rid of last entry

freqplain <- data.frame(term="",present=0,absent=0,stringsAsFactors = FALSE)
for (i in 1:length(importantGO)){
  freqplain[i,1] <- importantGO[i]; 
  temptable <- as.vector(table(tempplain[,i]))
  freqplain[i,2]  <- temptable[1]
  freqplain[i,3] <- temptable[2]
  }

temptarget <- freqtarget %>% dplyr:: select(starts_with(importantGO[1]))
for (i in 1:length(importantGO)){
  temp <- freqtarget %>% dplyr:: select(starts_with(importantGO[i]))
  temptarget <- cbind(temp,temptarget)
}
temptarget <- temptarget[,1:15] # get rid of last entry

freqtarget <- data.frame(term="",present=0,absent=0,stringsAsFactors = FALSE)
for (i in 1:length(importantGO)){
  freqtarget[i,1] <- importantGO[i]; 
  temptable <- as.vector(table(temptarget[,i]))
  freqtarget[i,2]  <- temptable[1]
  freqtarget[i,3] <- temptable[2]
}


freqtarget2 <- melt(freqtarget,id.var="term")
colnames(freqtarget2)[2] <- "GOterm"
ggplot(freqtarget2, aes(x = term, y = value, fill = GOterm)) + 
  theme(axis.text.x=element_text(face="bold",angle=40,hjust=1,size=12)) +
  theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
  ylab("Frequency counts of target terms") + 
  xlab("") +
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) +
  geom_bar(stat = "identity")

freqplain2 <- melt(freqplain,id.var="term")
freqplain2$term <- as.factor(freqplain2$term)
freqplain2$value <- as.integer(freqplain2$value)
colnames(freqplain2)[2] <- "GOterm"

ggplot(freqplain2, aes(term, value, fill = GOterm)) + 
  geom_bar(stat = "identity") +
  geom_rect(aes(xmin=0, xmax=16, ymin=500, ymax=1500), fill="white") +
  scale_y_continuous(limits=c(0,NA), breaks=(yticks), labels=yticks) +
  
  theme(axis.text.x=element_text(face="bold",angle=40,hjust=1,size=12)) +
  theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
  ylab("Frequency counts of non-target terms") + 
  xlab("")+
  theme(axis.title.y = element_text(color="black", size=14, face="bold")) 
  
###### break y-axis  ######
yticks <- c(0, 50, 100, 200, 500, 1000, 1500)
mycols <- c("red","red","red","red","red","red","red","red","red","red","red","red","red","red","red",
            "cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan","cyan")

gap.barplot(freqplain2$value,gap=c(500,520),ytics=c(0,50,100,200,250,300,400,500,                                                                 900,1000,1100,1200,1300,1400,1500),
            ylab="Frequency counts of non-target terms",main="", xlab="", xaxlab=character(30),
            col=mycols)
staxlab(1,1:30,freqplain2$term,srt=45)
#######
gap.barplot(freqtarget2$value,gap=c(450,500),ytics=c(0,100,200,300,400,450,500,600,700,800,900,1000,                                                                 900,1000,1100,1200,1300,1400,1500),
            ylab="Frequency counts of non-target terms",main="", xlab="", xaxlab=character(30),
            col=mycols)
staxlab(1,1:30,freqplain2$term,srt=45)

###############  select proteins that are nontargets but not in train or test set
# use trained RandomForest on these candidates for potential targets

allnontargets <- mcrap[mcrap$targets == 0,]
unknown <- negatives[!rownames(allnontargets) %in% rownames(negatives),]
unknown <- data.frame(unknown)

uindex <- base::sample(nrow(unknown),10) # indices of training samples
candidates <- unknown[sample(1:nrow(unknown), 100,replace=FALSE),] 
candidates <- unknown[1:2000,]
candidates <- data.frame(candidates)
targettype <- predict(rf_fit,unknown)
candidates <- data.frame(protein=rownames(unknown),target=targettype)















