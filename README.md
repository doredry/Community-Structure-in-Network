In complex networks, a network is said to have community structure if the nodes of the
network can be grouped into groups of nodes with dense connections internally, and sparser
connections between groups.

In this project, I have implemented an algorithm for detecting community structures (or
clusters) in a network. The ability to detect such groups is of significant importance.

For example, partitioning a protein-protein interaction network into clusters can provide a
modular view of the

Input file format:
The First value represents the number of nodes in the network, n = |V|.
The second value represents the number of edges of the first node, i.e., k1. It is followed by
the k1 indices of its neighbors, in increasing order.
The next value is k2, followed by the k2 indices of the neighbors of the second node, then k3
and its k3 neighbors, and so on until node n.

Output file format: 
The First value represents the number of groups in the division.
The second value represents the number of nodes in the first group, followed by the indices
of the nodes in the group, in increasing order.
The next value is the number of nodes in the second group, followed by the indices of the
nodes in group, then the number of nodes and indices of nodes in the third group, and so
on until the last group.

Example: cluster [input file name] [output file name]

Read more about the algorithm I used at: https://doi.org/10.1073/pnas.0601602103
