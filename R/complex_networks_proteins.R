# complex_networks_proteins.R

# Now analyze main protein classes comprising candidate targets
# Examine drug_targets dataframe

p_GPCR <- filter(drug_targets,TargetClass == "GPCR")
p_enzyme <- filter(drug_targets,TargetClass == "Enzyme")
p_IC <- filter(drug_targets,TargetClass == "Ion channel")
p_transporter <- filter(drug_targets,TargetClass == "Transporter")

# Alternatively, assess the drugs targetiing the k-core neighborhood.
