# complex_networks_data.R
# Load in data and preprocess where necessary

# HINT (High-quality INTeractomes) is a curated compilation of high-quality protein-protein 
# interactions from 8 interactome resources (BioGRID, MINT, iRefWeb, DIP, IntAct, HPRD, MIPS 
# and the PDB). Contains 12,429 unique proteins with 59,128 interactions between them. http://hint.yulab.org/

ppi_hint <- read.csv("c:\\R-files\\proteins\\HINT-2017.csv", header=TRUE,stringsAsFactors = FALSE)

drug_targets <- load_drugtargets()  # loads in and preprocesses the drug.targets file.


# Covert gene names all to uppercase
ppi_hint <- mutate_all(ppi_hint, funs(toupper))

# get rid of duplicate A and B columns - NB for some reason duplicated() function cannot find them!!!
ppi_hint <- ppi_hint %>% 
  dplyr::filter(Gene_A != Gene_B)

# add bioplex ppi data to hint ppi data
# http://bioplex.hms.harvard.edu/

bioplex <- read.csv(file="C://R-files//proteins//BioPlex.csv", header=TRUE, sep=",")
# Need to combine "drug_targets" ppi with "ppi_hint"
dt_genes <- unique(drug_targets$Gene)
ppigenes <- c(ppi_hint[,1],ppi_hint[,2])
ppigenes <- unique(ppigenes)

ppigenes <- c(ppigenes,bioplex[,1],bioplex[,2])

shite <- ppigenes[ppigenes %in% dt_genes]













