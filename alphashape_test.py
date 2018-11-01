import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import math as m


def add_edge(edges, i, j, outer_only=True):
    """
    Add a line between the i-th and j-th points,
    if not in the list already
    """
    if (i, j) in edges or (j, i) in edges:  # already added
        if outer_only: edges.remove((j, i))
        return

    edges.add((i, j))

tri = Delaunay(X)
tri_set = set(tuple(row) for row in tri.simplices)
edges = set()
edge_points = []

alpha = 4

for ia, ib, ic in tri.simplices:

    a = np.linalg.norm(X[ia]-X[ib])
    b = np.linalg.norm(X[ib]-X[ic])
    c = np.linalg.norm(X[ic]-X[ia])

    s = (a + b + c) / 2.0

    area = m.sqrt(s*(s-a)*(s-b)*(s-c))
    circum_r = a*b*c/(4.0*area)

    if circum_r < 1.0 / alpha:
        add_edge(edges, ia, ib)
        add_edge(edges, ib, ic)
        add_edge(edges, ic, ia)

boundary_set = set()
for edge in edges:
    boundary_set.add(tuple(X[edge[0], :]))
    boundary_set.add(tuple(X[edge[1], :]))

for edge in edges:
    plt.plot(X[edge, 0], X[edge, 1], '0.6')
plt.plot(X[:, 0], X[:, 1], 'o')
for point in boundary_set:
    plt.plot(point[0], point[1], 'ro')

plt.xlim((-1, 1))
plt.ylim((-1, 1))
plt.show()

