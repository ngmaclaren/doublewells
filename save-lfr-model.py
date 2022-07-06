import numpy as np
import networkx as nx
# import matplotlib.pyplot as plt

from networkx.generators.community import LFR_benchmark_graph as LFR

# These parameters are chosen to match the number of expected nodes (100) and expected density (0.06)
# used for the other model networks in this study. Additionally, some paramters, such as max and min
# community size, were chosen so that the model would converge.

# number of nodes
n = 100
# Power law exponent for the degree distribution
tau1 = 1.8
# Power law exponent for the size of the clusters/communities
tau2 = 1.04
# The probability of an edge crossing between communities
mu = 0.1
# The average degree
avg_degree = 5
# The maximum degree
max_degree = int(n/4)
# Maximum community size
max_community = max_degree
# Minimum community size
min_community = int(n/10)

# Generate the network with the above settings. 
g = LFR(
    n
    , tau1
    , tau2
    , mu
    , seed = 12345
    , max_degree = max_degree
    , average_degree = avg_degree
    , max_community = max_community
    , min_community = min_community
)

# Remove any self-loops.
g.remove_edges_from(nx.selfloop_edges(g))
print(nx.density(g))

# Visualize the network.
nx.draw(g)
plt.show()

# And save the adjacency matrix.
A = nx.adjacency_matrix(g)
with open("./data/LFR-NetworkX.txt", "wb") as f:
    for line in A.todense():
        np.savetxt(f, line, fmt='%.0f')
