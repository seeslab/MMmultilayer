# MMmultilayer
Mixed membership EM algorithms for multi-layered data from the manuscript Tensorial and bipartite block models for link prediction in layered networks and temporal networks (Tarrés-Deulofeu, Godoy-Lorite, Guimerà and Sales-Pardo)

The code for the tensorial model applied to multi-layered data with R types of interactions is:
TensorialMMSBM.py
Usage: pypy TensorialMMSBM.py traindata testdata #nodes #layers #groupsnodes #groupslayers #labels #initializations output

The code for the bipartite model applied to multi-layered data with R types of interactions is:
BipartiteMMSBM.py
Usage: pypy BipartiteMMSBM.py traindata testdata #nodes #layers #groupsnodes(K) #groupslayers(L) #labels(R) #initializations(for the mean) output

traindata - file with training data (see XXXtrainXXX.dat in the repository for an example and explanations below)
testdata - file with test data (see XXXtestXXX.dat in the repository for an example and explanations below)
#nodes - number of nodes in dataset
#layers - number of layers in dataset
#groupsnodes - number of groups of nodes K considered in the MMSBM
#groupslayers - number of groups of layers L considered in the dataset
#labels - number of possible interaction types, R
#initializations - number of maximization runs. The final output is the average over results for each run.
output - if 1, the code pirnts the model parameters for each one of the runs; if 0, it does not print any model parameters

Train/Test datafile specifications:
Nodes(links) and layers should be numbered consecutively [0,...,#nodes] , [0,...,#layers]
The format of the files should be a 4 column file:
node1 node2 layer interaction_type

Columns must be separated with simple spaces. The files drugstrain0.dat and drigstest0.dat provide an example of input file.

Code output: 
Probability of an link being of type 0,1,2 ...R for each link in the training dataset. For each link, the output probability is the average probability over the #initializations runs.


