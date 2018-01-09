# build a igraph object and get its statistics
# Do not plot the entire graph as its way too large - plot a small subset only
ghint <- graph.data.frame(ppi_hint)
ghint<-as.undirected(ghint); 
gs<-get_gstatistics(ghint)
head(gs)

# link salience - might solve hub idenification problem
# https://www.nature.com/articles/ncomms1847
# https://github.com/csgillespie/poweRlaw
# Calculate power law for our protein networks - any diff between hubs and non-hubs?
# also graph coreness https://jcasasr.wordpress.com/2015/02/03/plotting-the-coreness-of-a-network-with-r-and-igraph/
cores = graph.coreness(as.undirected(ghint))
head(sort(cores, decreasing=TRUE), 3)

# create network and plot it
layouts <- grep("^layout_",ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <-do.call(layout,list(net))
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }
