# build a igraph object and get its statistics
# Do not plot the entire graph as its way too large - plot a small subset only
ghint <- graph.data.frame(ppi_hint)
ghint<-as.undirected(ghint); 
gs<-get_gstatistics(ghint)
head(gs)


# Try and produce a power law graph of degree
m = displ$new(gs$degree)
##Estimate the cut-off
estimate_xmin(m)
m$setXmin(105); m$setPars(2.644)
plot(m,xlab="Degree",ylab="CCDF" ,panel.first = grid(lty = 6, col = "cornsilk2"))
lines(m, col=2)
abline(h=c(10,20,50), v=c(1,2,5,10,20,50,100,200,500), col="gray", lty=3)

# link salience - might solve hub idenification problem
# https://www.nature.com/articles/ncomms1847
# https://github.com/csgillespie/poweRlaw
# Calculate power law for our protein networks - any diff between hubs and non-hubs?
# also graph coreness https://jcasasr.wordpress.com/2015/02/03/plotting-the-coreness-of-a-network-with-r-and-igraph/
# get subgraphs based on coreness measure
subgraph <- make_ego_graph(ghint, order=1, c("RPSA", "RPS2", "RPL7"))
plot(subgraph[[3]])

coreness = graph.coreness(as.undirected(subgraph[[1]]))
head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness));

# get subgraphs based on top ranking coreness measure
subgraph <- make_ego_graph(ghint, order=1, c("RPSA", "RPS2", "RPL7","RPL10","TP53"))
plot(subgraph[[5]])

# create layout
ll <- corenessLayout(subgraph[[5]]);
# plot
plot(subgraph[[5]], layout=ll, vertex.size=15, vertex.color=colbar[coreness], vertex.frame.color=colbar[coreness], main='Coreness');

# Try different layouts for best looking graph
net <- subgraph[[1]]
# create network and plot it
layouts <- grep("^layout_",ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <-do.call(layout,list(net))
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }




