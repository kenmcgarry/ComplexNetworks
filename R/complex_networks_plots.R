# complex_networks_plots.R
# plot the pretty pictures used in the paper.


# Show the various types of protein targets and their freqs
bar_plot_drugtargets <- function(){
  # Will need cropping when saved as PDF - use https://pdfresizer.com/crop
  mar.default <- c(7,4,4,2) + 0.1
  par(mar = mar.default + c(0, 4, 0, 0)) 
  
  tempnames <- sort(table(drug_targets$TargetClass),decreasing = TRUE)
  mp <- barplot(sort(table(drug_targets$TargetClass),decreasing = TRUE),col="blue",space=1,
                ylab="Number of drug targets (proteins)",ylim=c(0,6000),main="",xlab = "",xaxt = "n")
                
  mytable <- length(table(drug_targets$TargetClass))
  end_point <- 0.5 + (mytable) + (mytable)-1
  
  text(seq(1.5,end_point,by=2), par("usr")[1],#-0.25, # -0.25
       srt = 45, adj= 1, xpd = TRUE,
       labels = paste(names(tempnames)),cex=0.9)
}


# Do not plot the entire graph as its way too large - plot a small subset only
plot_small <- function(){
  return(null)
ad <- get.adjacency(g)
nodesize=igraph::degree(g)
for (i in 1:length(nodesize))
  if(nodesize[i] <= 30){nodesize[i] <-10} else {nodesize[i]<-15}

nodecolor=character(ncol(ad))  # create a character for every column in adjaceny matrix
for (i in 1:length(nodecolor))
  if(nodesize[i] <= 10){nodecolor[i] <- "lightblue"} else {nodecolor[i]<-"darkgray"}
nodelabel<-V(g)$name
plot(gsmall, edge.color="darkgray", 
     vertex.color=nodecolor,
     vertex.label=nodelabel,
     vertex.label.cex=0.6, 
     vertex.label.font=0.5, 
     vertex.frame.color="darkgreen",
     vertex.size=18,
     vertex.label.color="black", 
     vertex.label.family = "sans", 
     layout=layout.kamada.kawai(g))
}


bar_plot_drugtargets()  # barplot of the protein types of drug targets

