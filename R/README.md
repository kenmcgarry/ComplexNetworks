The following files:

`complex_networks_main.R` is the starting point, it calls in the other R files, loading in functions and data files.

`complex_networks_functions.R` contains the workhorse functions and loads in R libraries.

`complex_networks_reviewers.R` The reviewers had suggested several improvements, we implemented these in full so best to examine this file first.

`complex_networks_plots.R` produces plots that appear in the paper.

`complex_networks_tables.R` produces tables that appear in the paper.

`complex_networks_data.R` reads the HINT data and target data in. Proprocess data and removes duplicate entries.

`complex_networks_buildnets.R` creates the igraph network structures and conducts statistical analysis on them. Looks for hub proteins and those proteins that are hubs AND targets.

`complex_networks_misc.R` contains unassigned functions and snippets of code that will be placed in 'complex_networks_functions.R` when fully debugged and working. 
