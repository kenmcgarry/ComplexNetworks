# complex_networks_functions.R
# helper functions and calls in libraries
# started: 5/1/18
library(sand)
library(igraph)
library(GO.db)
library(GOstats)
library(org.Sc.sgd.db)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(igraph)
library(NCBI2R)
library(rentrez)
library(ggplot2)
library(biomaRt)
library(scales)
library(grid)
library(RColorBrewer)
library(xtable)


# Makes first letter of string uppercase
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  return(x)
}

# Drug targets of the various drugs - usefully contains protein type (e.g. GPCR) as well.
load_drugtargets <- function(){
  drug_targets <- file.path('C://R-files//disease','drug.target.interaction.tsv.gz') %>% read.delim(na.strings='',header =TRUE,stringsAsFactors = FALSE) 
  names(drug_targets)[names(drug_targets)=="DRUG_NAME"] <- "DrugName"
  names(drug_targets)[names(drug_targets)=="TARGET_CLASS"] <- "TargetClass"
  names(drug_targets)[names(drug_targets)=="GENE"] <- "Gene"
  drug_targets <- drug_targets[,c(1,4,6)]  # Keep only need three variables
  drug_targets$DrugName <- firstup(drug_targets$DrugName)   # convert first letter to uppercase to match existing data
  drug_targets <- na.omit(drug_targets)# remove NA's
  
  return(drug_targets)
}

# R provides a tail and head command to view last six and first six elements, so why not the middle six?
middle <- function(mydata) {
  len <- nrow(mydata)
  startpoint <- round(len/2)
  endpoint <- startpoint+5
  mydata[startpoint:endpoint,]
  
}

# For each protein in the list findout how many publiactions it has been mentioned in
get_pubs <- function(protein_list){
  npubs <- count_articles(protein_list)
  return(npubs)
}


# See how many research articles are written about our proteins. Uses rentrez package.
count_articles <- function (protein_list){
  for (i in 1:length(protein_list)){
    
    print(protein_list[i])
    pname <- paste(protein_list[i],'[GENE]) AND (Homo sapiens[ORGN])',sep="")
    ids<-entrez_search(db="pubmed", term=pname,retmax=40000)
    atemp <- cbind(protein_list[i],length(ids$ids))
    
    if(i!=1){
      articles <- rbind(articles,atemp)} 
    else{
      articles <- atemp}
  }
  return(articles)
}

# Uses defeunct NCBI apckage, For each protein get the proteins they interact with.
get_interactions <- function(protein_list){
  
  for (i in 1:length(protein_list)){
    print(protein_list[i])
    pname <- paste(protein_list[i],'[sym]',sep="")
    ids<-GetIDs(pname)
    plist<-GetInteractions(ids)
    plist<-unique(plist[13])
    
    cvd<-rep(protein_list[i],nrow(plist))
    ptemp <- cbind(cvd,plist)
    
    if(i!=1){
      ppi <- rbind(ppi,ptemp)} 
    else{
      ppi <- ptemp}
  }
  return(ppi)
}

# Uses defunct NCBI package, Count how many interactions each protein has.
count_interactions <- function(protein_list) {
  for (i in 1:length(protein_list)){
    print(protein_list[i])
    pname <- paste(protein_list[i],'[sym]',sep="")
    ids<-GetIDs(pname)
    plist<-GetInteractions(ids[1])
    plist<-unique(plist[13])
    ptemp <- cbind(protein_list[i],nrow(plist))
    if(i!=1){
      ppi <- rbind(ppi,ptemp)} 
    else{
      ppi <- ptemp}
  }
  
  return(ppi)
}


# intial code from Kolaczyk and Csardi book, good for intro to GO annoations
example_from_kolaczyk <- function(){
  oldw <- getOption("warn")
  options(warn = -1)
  set.seed(42)
  load("ppi_CC.RData")
  upgrade_graph(ppi.CC)
  #save(ppi.CC, file = "ppi_CC.RData")
  
  summary(ppi.CC)
  V(ppi.CC)[ICSC == 1]$color <- "yellow"
  V(ppi.CC)[ICSC == 0]$color <- "lightblue"
  
  plot(ppi.CC, vertex.size=5, vertex.label=NA,
       main="Network of interacting proteins responsible for cell communications: \n yellow=involved in ICSC, blue=not involved")
  
  clu<- clusters(ppi.CC)
  ppi.CC.gc <- induced.subgraph(ppi.CC,clu$membership==which.max(clu$csize))
  nn.ave <- sapply(V(ppi.CC.gc), function(x) mean(V(ppi.CC.gc)[nei(x)]$ICSC))
  
  par(mfrow=c(2,1))
  hist(nn.ave[V(ppi.CC.gc)$ICSC ==1], col="yellow", ylim=c(0,30),xlab="Proportion neighbours w/ ICSC",
       main="Egos w/ ICSC")
  hist(nn.ave[V(ppi.CC.gc)$ICSC ==0], col="lightblue", ylim=c(0,30),xlab="Proportion neighbours without ICSC",
       main="Egos w/out ICSC")
  
  nn.pred <- as.numeric(nn.ave > 0.5)
  mean(as.numeric(nn.pred) != V(ppi.CC.gc)$ICSC)
  
  # now do the GO stuff, as per page 140 of Kolaczyk
  x <- as.list(org.Sc.sgdGO2ALLORFS)
  current.icst <- x[names(x) == "GO:0035556"]
  ev.code <- names(current.icst[[1]])
  icst.ida <- current.icst[[1]][ev.code == "IDA"]
  orig.icsc <- V(ppi.CC.gc)[ICSC == 1]$name
  candidates <- intersect(icst.ida, V(ppi.CC.gc)$name)
  new.icsc <- setdiff(candidates, orig.icsc)
  new.icsc
  
  options(warn = oldw)
  
}










