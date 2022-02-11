---
layout: default
title: Output
nav_order: 5
has_children: true
---

# Output

This contains the finalized PDB files, plots, snapshots, and presentations I've given to share my results.
Files are formatted as follows:

## PDBs
Files named: '[pdb_id]_[n]_domains.pdb' contain the results of the Domain Subdivision for n domains.
The labels of the domains are assigned to the b-factor values of the pdb for convenient visualization.

## Plots
Files named: '[pdb_id]_[n]_domains.png' are plots of the scores of each attempted subdivision for each pdb id.
The n in this case will be that of the highest scoring subdivision.

## Final Clustering
Files named: '[pdb_id]\_[n]\_domains_optimal.png' are snapshots of the highest scoring subdivision of the pdb
visualized in ChimeraX. An Icosahedral structure is overlaid to indicate the choice of subdivision. The ChimeraX sessions
used to create these screenshots are saved under Subdivisions/Sessions with the format: '[pdb_id]_[n]_domains.cxs'


# Summary Of HK97 Results

<object data="./Presentations/HK97_Results.pdf" width="1000" height="1000" type='application/pdf'></object>

