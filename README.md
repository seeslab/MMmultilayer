# MMmultilayer
Mixed membership EM algorithms for multi-layered data from the manuscript Tensorial and bipartite block models for link prediction in layered networks and temporal networks (Tarrés-Deulofeu, Godoy-Lorite, Guimerà and Sales-Pardo)

In the manuscript we analyze two types of datasets: 1) a drug-drug interaction dataset in different cancer cell lines; 2) time-resolved e-mail communication network dataset.

For 1) links (drug-drug interactions) can be of three different types (0,1,2)

The code for the tensorial model applied to multi-layered data with three types of interactions is:
nodesdrugsinference.py
Usage: 

The code for the bipartite model applied to multi-layered data with three types of interactions is:
linksdrugsinference.py

Usage:

For 2) links (email communicatio between pairs of nodes) can be either active or inactive

The code for the tensorial model applied to multi-layered data with two types of interactions is:
timeinference.py

Usage:

The code for the bipartite model applied to multi-layered data with two types of interactions is:
linkstimeinference.py

Usage:

