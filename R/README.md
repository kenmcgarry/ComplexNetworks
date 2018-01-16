The following files:

`complex_networks_main.R` is the starting point. It loads in R libraries and calls in the other R files, loading in functions and data files.

`complex_networks_functions.R` contains the workhorse functions.

`complex_networks_plots.R` produces plots that appear in the paper.

`complex_networks_data.R` reads the HINT data and target data in. Proprocess data and removes duplicate entries.

`complex_networks_buildnets.R` creates the igraph network structures and conducts statistical analysis on them. Looks for hub proteins and those proteins that are hubs AND targets.
