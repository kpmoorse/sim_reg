import numpy as np
from scipy.spatial import Delaunay
import math as m


def alphashape_2D(X, alpha, outform="pt_idx"):

    tri = Delaunay(X)
    edges = set()

    for ia, ib, ic in tri.simplices:

        a = np.linalg.norm(X[ia] - X[ib])
        b = np.linalg.norm(X[ib] - X[ic])
        c = np.linalg.norm(X[ic] - X[ia])

        s = (a + b + c) / 2.0

        area = m.sqrt(s * (s - a) * (s - b) * (s - c))
        circum_r = a * b * c / (4.0 * area)

        if circum_r < 1.0 / alpha:
            add_edge(edges, ia, ib)
            add_edge(edges, ib, ic)
            add_edge(edges, ic, ia)

    if outform == "pt_idx":
        output = set()
        for edge in edges:
            output.add(edge[0])
            output.add(edge[1])
    elif outform == "pt_val":
        output = set()
        for edge in edges:
            output.add(tuple(X[edge[0], :]))
            output.add(tuple(X[edge[1], :]))

    return output


def add_edge(edges, i, j, outer_only=True):
    """
    Add a line between the i-th and j-th points,
    if not in the list already
    """
    if (i, j) in edges or (j, i) in edges:  # already added
        if outer_only: edges.remove((j, i))
        return

    edges.add((i, j))