# complex_networks_plots.R
# plot the pretty pictures used in the paper.

# ggplot2 vewrsion of barplot - prettier?
bar_plot_gg2 <- function(dt,br,mycolor){
  tempnames <- sort(table(dt$TargetClass),decreasing = TRUE)
  Frequency <- as.vector(tempnames)  # counts for each protein type
  pnames <- names(tempnames)
  df <- data.frame(Frequency,pnames)
  df$Frequency <- as.numeric(df$Frequency)

  if(br==1){mybreaks = c(0,500,1000,1500,2000,2500,3000,3500,4000);limits<-max(mybreaks)}
  if(br==2){mybreaks = c(0,200,400,600,800,1000);limits<-max(mybreaks)}
  
  ggplot(df, aes(x = pnames,y=(Frequency))) + scale_x_discrete(limits = pnames) + 
    geom_bar(stat="identity",fill=mycolor)+
    theme(axis.text.x=element_text(face="bold",angle=40,hjust=1,size=12)) +
    theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
    ylab("Frequency count of proteins") + 
    xlab("")+
    scale_y_continuous(expand = c(0,0),breaks = mybreaks,limits = c(0,limits)) +
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
plot_power <- function(gs){
  m = displ$new(gs$degree)
  ##Estimate the cut-off
  estimate_xmin(m)
  m$setXmin(105); m$setPars(2.644)
  
  temp <- plot(m,xlab="Degree",ylab="Percentiles") # only use this to grab coordinates for ggplot2
  
  ggplot(data=temp,aes(x=x,y=y*100))+
    geom_point()+
    geom_smooth(se = FALSE, method = "gam", formula = y ~ s(log(x)))+
    labs(title="",x="Degree distribution",y="Percentiles")+
    theme(axis.text.x=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
    theme(axis.title.y = element_text(color="black", size=14, face="bold"))+
    theme(axis.title.x = element_text(color="black", size=14, face="bold"))+
    geom_vline(xintercept = 5, linetype="dashed", color = "red", size=1)
  
  cat("\nPercentiles",quantile(gs$degree, c(.70, .80, .90)) )

}


# ggplot2 version of power law plot
plot_power2 <- function(gppi){
  
G.degrees <- degree(gppi)
G.degree.histogram <- as.data.frame(table(G.degrees))
G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
  geom_point(alpha = 0.5, color = "blue",size = 3)+
  #labs(title="",x="Degree distribution",y="Percentiles")+
  theme(axis.text.x=element_text(face="bold",angle=0,hjust=1,size=12)) +
  theme(axis.text.y=element_text(face="bold",angle=0,hjust=1,size=12)) +
  theme(axis.title.y = element_text(color="black", size=14, face="bold"))+
  theme(axis.title.x = element_text(color="black", size=14, face="bold"))+
  
  scale_x_continuous("Degree",expand = c(0,0),
                     breaks = c(1, 3, 10, 30, 100, 200, 500),
                     trans = "log10") +
  scale_y_continuous("Protein Frequency",expand=c(0,0),
                     breaks = c(1, 3, 10, 30, 100, 300, 1000, 2000),
                     trans = "log10") +
  ggtitle("Degree Distribution (log-log)") +
  geom_segment(aes(x = 1, y = 300, xend = 15, yend = 300),color="red",linetype="dashed",size=1)  + # Horiz
  geom_segment(aes(x = 15, y = 1, xend = 15, yend = 300),color="red",linetype="dashed",size=1)    # Vert
}



