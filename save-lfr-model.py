import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from networkx.generators.community import LFR_benchmark_graph as LFR

n = 100
tau1 = 1.8
tau2 = 1.04
mu = 0.1
avg_degree = 5
max_degree = int(n/4)
max_community = max_degree
min_community = int(n/10)

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

g.remove_edges_from(nx.selfloop_edges(g))
print(nx.density(g))

nx.draw(g)
plt.show()

A = nx.adjacency_matrix(g)

with open("./data/LFR-NetworkX.txt", "wb") as f:
    for line in A.todense():
        np.savetxt(f, line, fmt='%.0f')
