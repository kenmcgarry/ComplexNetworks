# ComplexNetworks
## Updated 14th May 2018

In this work we use complex network theory to provide a statistical model of the connectivity patterns of human proteins and their interaction partners. Our intention is to identify important proteins that may be predisposed to be potential candidates as drug targets for therapeutic interventions. Target proteins usually have more interaction partners than non-target proteins, but there are no hard-and-fast rules for defining the actual number of interactions. We devise a statistical measure for identifying hub proteins, we score our target proteins with gene ontology annotations. The important drugable protein targets are likely to have similar biological functions that can be assessed for their potential therapeutic value. Our system provides a statistical analysis of the local and distant neighborhood protein interactions of the potential targets using complex network measures. This approach builds a more accurate model of drug-to-target activity and therefore the likely impact on treating diseases. We integrate high quality protein interaction data from the HINT database and disease associated proteins from the DrugTarget database. Other sources include biological knowledge from Gene Ontology and drug information from DrugBank. The problem is a very challenging one since the data is highly imbalanced between target proteins and the more numerous nontargets. We use undersampling on the training data and build Random Forest classifier models which are used to identify previously unclassified target proteins. We validate and corroborate these findings from the available literature.


If you find the source code useful please cite our paper:

K. McGarry and S. Mcdonald (2018) Complex network theory for the identification and assessment of candidate protein targets. Computers in Biology and Medicine, 97. pp. 113-123.
