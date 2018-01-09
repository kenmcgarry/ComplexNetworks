# create network and plot it
layouts <- grep("^layout_",ls("package:igraph"), value=TRUE)[-1]
# Remove layouts that do not apply to our graph.
layouts <- layouts[!grepl("bipartite|merge|norm|sugiyama|tree", layouts)]
par(mfrow=c(3,3), mar=c(1,1,1,1))
for (layout in layouts) {
  print(layout)
  l <-do.call(layout,list(net))
  plot(net, edge.arrow.mode=0, layout=l, main=layout) }
