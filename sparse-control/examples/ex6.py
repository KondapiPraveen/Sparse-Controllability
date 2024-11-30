import numpy as np
import networkx as nx
import math

n = 20
# p = 2*math.log(n)/n;
p = 0.89

G = nx.erdos_renyi_graph(n, p) # Undirected ER Graph
ADJ = nx.to_numpy_array(G) # Adjacency Matrix
L = np.diag(ADJ.sum(axis=1)) - ADJ # Laplace of the Graph
A = np.eye(n) - (1/n)*(L) # Dynamics
A = 2*A
B = 10*np.eye(n) # Input Matrix

#print(ADJ)
