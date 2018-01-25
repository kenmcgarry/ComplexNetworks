# complex_networks_data.R
# Load in data and preprocess where necessary

# HINT (High-quality INTeractomes) is a curated compilation of high-quality protein-protein 
# interactions from 8 interactome resources (BioGRID, MINT, iRefWeb, DIP, IntAct, HPRD, MIPS 
# and the PDB). Contains 12,429 unique proteins with 59,128 interactions between them. http://hint.yulab.org/

ppi_hint <- read.csv("c:\\R-files\\proteins\\HINT-2017.csv", header=TRUE,stringsAsFactors = FALSE)

drug_targets <- load_drugtargets()  # loads in and preprocesses the drug.targets file.
drug_targets <- drug_targets[!(duplicated(drug_targets[c("DrugName","Gene")]) | duplicated(drug_targets[c("DrugName","Gene")], fromLast = TRUE)), ]


# Covert gene names all to uppercase
ppi_hint <- mutate_all(ppi_hint, funs(toupper))

# get rid of duplicate A and B columns - NB for some reason duplicated() function cannot find them!!!
ppi_hint <- ppi_hint %>% 
  dplyr::filter(Gene_A != Gene_B)

# add bioplex ppi data to hint ppi data
# http://bioplex.hms.harvard.edu/

bioplex <- read.csv(file="C://R-files//proteins//BioPlex.csv", header=TRUE, sep=",")
ppi <- rbind(ppi_hint,bioplex) 

ppi <- ppi %>% 
  dplyr::filter(Gene_A != Gene_B)

rm(bioplex,ppi_hint)












