### There are a number of issues with some changes to either package updates or more than likely some lack of robustness in my own code. I have tracked the issues down to the ontology annotation process. In order to build Random Forest classifers we need each gene to be annotated with the goslim obo file. 

### In the meantime (as a fix) I have uploaded an old RDATA workspace (15thJanuary2019.RData). Load this in and the relevant data structure is "mt". From this point use the code in `complex_networks_reviewers.R` this should allow you train Random Forests etc on the data. 

### I will try and figure out why the annoation process has failed! Any help will be apreciated!

### The following files:

`complex_networks_main.R` is the starting point, it calls in the other R files, loading in functions and data files.

`complex_networks_functions.R` contains the workhorse functions and loads in R libraries.

`GO_slim_functions.R` I have modified the goslim functions to be more generic and not be hard coded for a specific obo file.

`complex_networks_reviewers.R` The reviewers had suggested several improvements, we implemented these in full so best to examine this file at some point.

`complex_networks_plots.R` produces plots that appear in the paper.

`complex_networks_tables.R` produces tables that appear in the paper.

`complex_networks_data.R` reads the HINT data and target data in. Proprocess data and removes duplicate entries.

`complex_networks_buildnets.R` creates the igraph network structures and conducts statistical analysis on them. Looks for hub proteins and those proteins that are hubs AND targets.

`complex_networks_misc.R` contains unassigned functions and snippets of code that will be placed in 'complex_networks_functions.R` when fully debugged and working. 
