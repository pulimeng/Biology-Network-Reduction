# Biology Network Reduction for Graph-Based Machine Learning

This repo provides algorithms and examples for the network reduction for different biological networks (protein-protein interaction, DNA, drug-target interaction). The reduction is necessary for the machine learning applications since most of these networks are large and sparse (both topological and feature-wise). The size of the networks makes them difficult for certain applications (graph classification). And the sparsity in features also poses trouble for the machine learning models since the useful information can be well scattered in the network. Imagine looking at two images with only handful of pixels differen and try to find the underlying difference. Thus it is important to reduce the networks before they can be used for any graph-based machine learning.

This README file is written by Limeng Pu.

# Dependencies

Change line 12 of the file accordingly then run `conda env create -f environment.yml`. 

# Biology Networks

1. Protein-protein interaction (PPI) network

    This network describes the interaction between the proteins. It comes from STRING database (https://string-db.org/). Version 11 is      currently in use. It is established by literature mining. Each interaction has a score, which is the confidence score derived by the    authors. The higher the score, more likely the interaction exists in reality.
    - Number of nodes (proteins): 19,354.
    - Number of edges (interactions): 11,759,454.
    - Format: protein1, protein2, score.
    - Direction: undirected graph
