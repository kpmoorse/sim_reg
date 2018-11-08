from alphashape import *
import numpy as np
import matplotlib.pyplot as plt
import math as m

# Generate concave data set
X = []
for _ in range(500):
    flag = 0
    while not flag:
        x1 = np.random.randint(1, 4)
        x2 = np.random.randint(1, 4)
        if x1 != 2 and x2 == 2: continue
        else: flag = 1

    x1 += np.random.random() - 0.5
    x2 += np.random.random() - 0.5
    X.append([x1, x2])
X = np.array(X)
X = (X-2)/2  # Rescale to fit on [-1, 1]

boundary_set, edges = alphashape_2D(X, 8)

for edge in edges:
    plt.plot(X[edge, 0], X[edge, 1], '0.6')
plt.plot(X[:, 0], X[:, 1], 'o')
for point in boundary_set:
    plt.plot(point[0], point[1], 'ro')
plt.xlim((-1, 1))
plt.ylim((-1, 1))
plt.show()