# GO_slim_functions.R
# 12/04/2018
library(ontologySimilarity)
library(ontologyIndex)
library(ontologyPlot)
library(GSEABase)

data(go)
data(gene_GO_terms)
data(GO_IC)


# The following functions are to avoid the errors kicked out by goSlim (and crashing R) 
# when it cant find any terms
error_go_mf <- function(){
  return("MF: Error detected!")
  #go_default <- "error"
}
error_go_cc <- function(){
  return("CC: Error detected!")
  #go_default <- "error"
}
error_go_bp <- function(){
  return("BP: Error detected!")
  #go_default <- "error"
}

# receives a single goslim ontology (there are several) and returns a list of column names, used to annotate 
# the training matrix. Will comprise a GO term number, database type (MF,BF, CC), and short text description
getnames_goslim <- function(gslim){
  
 #gslim <- get_ontology("C://R-files//GOslim//goslim_generic.obo")
 mf <- get_descendants(gslim, c("GO:0003674")) # MF terms
 bp <- get_descendants(gslim, c("GO:0008150")) # BP terms
 cc <- get_descendants(gslim, c("GO:0005575")) # CC terms
 cat("\nFound", length(mf) + length(bp) + length(cc), "terms.")  
 
 mf_prop <- vector(mode="character", length=length(mf))
 bp_prop <- vector(mode="character", length=length(bp))
 cc_prop <- vector(mode="character", length=length(cc))
 
 for (i in 1:length(mf)){
   temp <- get_term_property(ontology=gslim, property="ancestors", term=mf[i], as_names=TRUE)
   if(length(temp)==1)
     mf_prop[i] <- temp[1]
   if(length(temp)==2)
     mf_prop[i] <- temp[2]
   if(length(temp)==3)
     mf_prop[i] <- temp[3]
   if(length(temp)==4)
     mf_prop[i] <- temp[4]
   if(length(temp)==5)
     mf_prop[i] <- temp[5]
   if(length(temp)==6)
     mf_prop[i] <- temp[6]
   if(length(temp)==7)
     mf_prop[i] <- temp[7]
   if(length(temp)==8)
     mf_prop[i] <- temp[8]
   if(length(temp)==9)
     mf_prop[i] <- temp[9]
 }
 
 for (i in 1:length(bp)){
   temp <- get_term_property(ontology=gslim, property="ancestors", term=bp[i], as_names=TRUE)
   if(length(temp)==1)
     bp_prop[i] <- temp[1]
   if(length(temp)==2)
     bp_prop[i] <- temp[2]
   if(length(temp)==3)
     bp_prop[i] <- temp[3]
   if(length(temp)==4)
     bp_prop[i] <- temp[4]
   if(length(temp)==5)
     bp_prop[i] <- temp[5]
   if(length(temp)==6)
     bp_prop[i] <- temp[6]
   if(length(temp)==7)
     bp_prop[i] <- temp[7]
   if(length(temp)==8)
     bp_prop[i] <- temp[8]
   if(length(temp)==9)
     bp_prop[i] <- temp[9]
   if(length(temp)==10)
     bp_prop[i] <- temp[10]
   if(length(temp)==11)
     bp_prop[i] <- temp[11]
   if(length(temp)==12)
     bp_prop[i] <- temp[12]
   if(length(temp)==13)
     bp_prop[i] <- temp[13]
   if(length(temp)==14)
     bp_prop[i] <- temp[14]
 }
 
 for (i in 1:length(cc)){
   temp <- get_term_property(ontology=gslim, property="ancestors", term=cc[i], as_names=TRUE)
   if(length(temp)==1)
     cc_prop[i] <- temp[1]
   if(length(temp)==2)
     cc_prop[i] <- temp[2]
   if(length(temp)==3)
     cc_prop[i] <- temp[3]
   if(length(temp)==4)
     cc_prop[i] <- temp[4]
   if(length(temp)==5)
     cc_prop[i] <- temp[5]
   if(length(temp)==6)
     cc_prop[i] <- temp[6]
   if(length(temp)==7)
     cc_prop[i] <- temp[7]
   if(length(temp)==8)
     cc_prop[i] <- temp[8]
   if(length(temp)==9)
     cc_prop[i] <- temp[9]
   if(length(temp)==10)
     cc_prop[i] <- temp[10]
   if(length(temp)==11)
     cc_prop[i] <- temp[11]
   if(length(temp)==12)
     cc_prop[i] <- temp[12]
    }
 
 mf <- paste(mf,":MF:",mf_prop,sep="")
 bp <- paste(bp,":BP:",bp_prop,sep="")
 cc <- paste(cc,":CC:",cc_prop,sep="")
 thenames <- c(mf,bp,cc)
 return(thenames) 
 
}


# annotate_go() will receive a list of proteins and annotate with GO terms it will return
# a matrix of terms and proteins. Uses the GO data by Daniel Greene.
# This fucntion is modified from complexnetworks.
annotate_with_go <- function(plist){
  category <- c("MF","BP","CC")
  nproteins <- length(plist)
  
  tempgo <- gene_GO_terms[plist]  # Annotate!!
  tempgo <- tempgo[!sapply(tempgo, is.null)]  # Not all proteins have GO annotations so remove them.
  go_proteins <- names(tempgo) # Unfortunately, we are left with only 13,417 proteins.
  # sort GO terms by the three categories, breakdown maybe useful at later date for summary statistics
  cc <- go$id[go$name == "cellular_component"]
  bp <- go$id[go$name == "biological_process"]
  mf <- go$id[go$name == "molecular_function"] 
  temp_cc <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=cc, x))
  temp_bp <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=bp, x))
  temp_mf <- lapply(tempgo, function(x) intersection_with_descendants(go, roots=mf, x))
  
  return(tempgo)
}

# go_slim_annotation() receives a list of protein names and the goslim ontology to use.
# Purpsoe with goslim is to reduce the complexities of numerous GO annoatations into a few key terms.
# http://www.geneontology.org/page/go-slim-and-subset-guide#On_the_web. This uses the GSEABase package.
# WARNING: this stage may take a looong time to compute....the more genes the longer it takes.
go_slim_annotation <- function(mylist,obo){
  gostuff <- annotate_with_go(mylist)
  gostuff <-gostuff[names(which(lapply(gostuff, length) >1))]  # keep only genes with 1 or more GO annotations
  assign("go_mf", TRUE, env=globalenv())  # global variables to enable recovery from errors generated by goSlim()
  assign("go_bp", TRUE, env=globalenv())
  assign("go_cc", TRUE, env=globalenv())
  
  # load in the goslim ontology of your choice
  gslim_ont <- get_ontology(obo)
  tempnames <- getnames_goslim(gslim_ont)
  
  # create the matrix for classification algorithms, if a GO term is present mark it with by "1" in that column
  mm <- matrix(0,length(tempnames), length(names(gostuff)))  # mm=Number of GO terms in GoSlim (CC+BP+MF) x Number of genes
  colnames(mm) <- names(gostuff)  # protein names
  rownames(mm) <- tempnames
  
  slim <- getOBOCollection(obo,evidenceCode="ANY")
  
  for (i in 1:length(names(gostuff))){
    myCollection <- GOCollection(gostuff[[i]])
    genename <- names(gostuff[i])
    go_mf <- tryCatch(goSlim(myCollection, slim, "MF"),error=function(e) {go_mf <- error_go_mf()})
    go_bp <- tryCatch(goSlim(myCollection, slim, "BP"),error=function(e) {go_bp <- error_go_bp()})
    go_cc <- tryCatch(goSlim(myCollection, slim, "CC"),error=function(e) {go_cc <- error_go_cc()})
    
    # get lengths of MF, BP and CC
    mf <- get_descendants(gslim_ont, c("GO:0003674")) # MF terms
    bp <- get_descendants(gslim_ont, c("GO:0008150")) # BP terms
    cc <- get_descendants(gslim_ont, c("GO:0005575")) # CC terms
    
    if(length(go_mf) ==1) {mm[1:length(mf),i]  <- as.vector(matrix(0,ncol=length(mf)))} else{
      go_mf[go_mf$Count != 0,]$Count <- 1; # convert non-zero numbers into 1's
      mm[1:length(mf),i] <- go_mf$Count}            # found MF annotations, assign to matrix
    
    if(length(go_bp)==1) {mm[(length(mf)+1):(length(bp)+length(mf)),i] <- as.vector(matrix(0,ncol=length(bp)))} else{
      go_bp[go_bp$Count != 0,]$Count <- 1
      mm[(length(mf)+1):(length(bp)+length(mf)),i] <- go_bp$Count}
   
    if(length(go_cc) == 1) {mm[(length(bp)+length(mf)+1):(length(tempnames)),i] <-as.vector(matrix(0,ncol=length(cc)))} else{
      go_cc[go_cc$Count != 0,]$Count <- 1
      mm[(length(bp)+length(mf)+1):(length(tempnames)),i] <- go_cc$Count}
  }
  return(mm)  # return matrix  
}


# go_slim_annotation1() receives a list of protein names and the goslim ontology to use.
# Purpose with goslim is to reduce the complexities of numerous GO annoatations into a few key terms.
# http://www.geneontology.org/page/go-slim-and-subset-guide#On_the_web. This uses the GSEABase package.
# WARNING: this takes a looong time to compute..... approx 3 hours on my laptop
go_slim_annotation1 <- function(mylist,obo){
  gostuff <- annotate_with_go(mylist)
  gostuff <-gostuff[names(which(lapply(gostuff, length) >1))]  # keep only genes with 1 or more GO annotations
  assign("go_mf", TRUE, env=globalenv())  # global variables to enable recovery from errors generated by goSlim()
  assign("go_bp", TRUE, env=globalenv())
  assign("go_cc", TRUE, env=globalenv())
  
  # load in the goslim ontology of your choice
  gslim_ont <- get_ontology(obo)
  tempnames <- getnames_goslim(gslim_ont)
  
  # create the matrix for classification algorithms, if a GO term is present mark it with by "1" in that column
  mm <- matrix(0,length(tempnames), length(names(gostuff)))  # mm=Number of GO terms in GoSlim (CC+BP+MF) x Number of genes
  colnames(mm) <- names(gostuff)  # protein names
  rownames(mm) <- tempnames
  
  slim <- getOBOCollection(obo,evidenceCode="ANY")
  
  for (i in 1:length(names(gostuff))){
    myCollection <- GOCollection(gostuff[[i]])
    genename <- names(gostuff[i])
    go_mf <- tryCatch(goSlim(myCollection, slim, "MF"),error=function(e) {go_mf <- error_go_mf()})
    go_bp <- tryCatch(goSlim(myCollection, slim, "BP"),error=function(e) {go_bp <- error_go_bp()})
    go_cc <- tryCatch(goSlim(myCollection, slim, "CC"),error=function(e) {go_cc <- error_go_cc()})
    
    # get lengths of MF, BP and CC
    mf <- get_descendants(gslim_ont, c("GO:0003674")) # MF terms
    bp <- get_descendants(gslim_ont, c("GO:0008150")) # BP terms
    cc <- get_descendants(gslim_ont, c("GO:0005575")) # CC terms
    
    if(length(go_mf) ==1) {mm[1:length(mf),i]  <- as.vector(matrix(0,ncol=length(mf)))} else{
      go_mf[go_mf$Count != 0,]$Count <- 1; # convert non-zero numbers into 1's
      mm[1:length(mf),i] <- go_mf$Count}            # found MF annotations, assign to matrix
    
    if(length(go_bp)==1) {mm[(length(mf)+1):(length(bp)+length(mf)),i] <- as.vector(matrix(0,ncol=length(bp)))} else{
      go_bp[go_bp$Count != 0,]$Count <- 1
      mm[(length(mf)+1):(length(bp)+length(mf)),i] <- go_bp$Count}
    
    if(length(go_cc) == 1) {mm[(length(bp)+length(mf)+1):(length(tempnames)),i] <-as.vector(matrix(0,ncol=length(cc)))} else{
      go_cc[go_cc$Count != 0,]$Count <- 1
      mm[(length(bp)+length(mf)+1):(length(tempnames)),i] <- go_cc$Count}
  }
  return(mm)  # return matrix  
}




