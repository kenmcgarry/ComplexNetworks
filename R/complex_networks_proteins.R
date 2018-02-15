# complex_networks_proteins.R

# Now analyze main protein classes comprising candidate targets
# Examine drug_targets dataframe

GPCR <- filter(drug_targets,TargetClass == "GPCR")
enzyme <- filter(drug_targets,TargetClass == "Enzyme")
IC <- filter(drug_targets,TargetClass == "Ion channel")
transporter <- filter(drug_targets,TargetClass == "Transporter")


