import numpy as np
import networkx as nx
import math

n = 100
p = 2*math.log(n)/n;

G = nx.erdos_renyi_graph(n, p) # Undirected ER Graph
ADJ = nx.to_numpy_array(G) # Adjacency Matrix
L = np.diag(ADJ.sum(axis=1)) - ADJ # Laplace of the Graph
A = np.eye(n) - (1/n)*(L) # Dynamics
B = np.eye(n) # Input Matrix

print(ADJ)
