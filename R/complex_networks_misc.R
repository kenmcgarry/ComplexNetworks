# complex_networks_misc.R
# part completed functions and code that has yet to be assigned.
# Not called up as yet, as will cause errors.

######################### CORENESS OF NETWORKS ######################################
# link salience - might solve hub idenification problem
# https://www.nature.com/articles/ncomms1847
# https://github.com/csgillespie/poweRlaw
# Calculate power law for our protein networks - any diff between hubs and non-hubs?
# also graph coreness https://jcasasr.wordpress.com/2015/02/03/plotting-the-coreness-of-a-network-with-r-and-igraph/
# get subgraphs based on coreness measure

subgraph_ppi <- make_ego_graph(ppi_net, order=1, c("RPSA", "PNMA1", "MZT2B","BRCA2","TUBGCP4"))
subgraph_ppi <- make_ego_graph(ppi_net, order=1, c("PSMA3"))


sg <- 4;
coreness <- graph.coreness(as.undirected(ppi_net))
coreness <- graph.coreness(as.undirected(subgraph_ppi[[sg]]))
coreness <- graph.coreness(as.undirected(subgraph_ppi))

head(sort(coreness, decreasing=TRUE))
colbar <- rainbow(max(coreness));
ll <- corenessLayout(subgraph_ppi[[sg]]);
plot(subgraph_ppi[[sg]], layout=ll, vertex.size=15, vertex.color=colbar[coreness], 
     vertex.label.color= "black", vertex.frame.color=colbar[coreness], main='');


# Try different layouts for best looking graph
net <- subgraph_ppi[[sg]]
# create network and plot it
layouts <- grep("^layout_",ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(4,4), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <-do.call(layout,list(net))
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }


