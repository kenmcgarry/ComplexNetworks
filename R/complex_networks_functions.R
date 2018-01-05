# complex_networks_functions.R
# helper functions and calls in libraries
# started: 5/1/18

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

get_pubs <- function(protein_list){
  npubs <- count_articles(protein_list)
  #nacts <- count_interactions(c("APP","TP53"))
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














