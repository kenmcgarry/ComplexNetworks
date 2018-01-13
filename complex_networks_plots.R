# complex_networks_plots.R
# plot the pretty pictures used in the paper.

# ggplot2 vewrsion of barplot - prettier?
bar_plot_gg2 <- function(){
  tempnames <- sort(table(drug_targets$TargetClass),decreasing = TRUE)
  Frequency <- as.vector(tempnames)
  pnames <- names(tempnames)
  df <- data.frame(Frequency,pnames)
  df$Frequency <- as.numeric(df$Frequency)

  ggplot(df, aes(x = pnames,y=(Frequency))) + scale_x_discrete(limits = pnames) + 
    geom_bar(stat="identity",fill="red")+
    theme(axis.text.x=element_text(face="bold",angle=40,hjust=1,size=12)) +
    theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
    ylab("Frequency count of proteins") + 
    xlab("")+
    scale_y_continuous(expand = c(0,0),breaks = c(0,500,1000,1500,2000,2500,3000,3500,4000))+
    theme(axis.title.y = element_text(color="black", size=14, face="bold"))
   
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

# produce a power law graph of degree for paper
plot_power <- function(){
  m = displ$new(gs$degree)
  ##Estimate the cut-off
  estimate_xmin(m)
  m$setXmin(105); m$setPars(2.644)
  
  temp <- plot(m,xlab="Degree",ylab="CDF") # only use this to grab coordinates for ggplot2
  
  ggplot(data=temp,aes(x=x,y=y))+
    geom_point()+
    geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)))+
    labs(title="",x="Degree distribution",y="CDF")+
    theme(axis.text.x=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.title.y = element_text(color="black", size=14, face="bold"))+
    theme(axis.title.x = element_text(color="black", size=14, face="bold"))
  
}

