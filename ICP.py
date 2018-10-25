import numpy as np
import matplotlib.pyplot as plt
import math as m
import random
from sklearn.neighbors import NearestNeighbors
import cv2

# Generate ideal target set
n_pts = 10
# targetSet = np.zeros((n_pts, 3))
# for i in range(n_pts):
#     th = random.random()*2*m.pi
#     targetSet[i, :] = np.array([m.cos(th), m.sin(th), 1])
targetSet = np.concatenate((np.random.random((n_pts, 2))*2-1, np.ones((n_pts, 1))), axis=1)

# Define transformation matrices
T_r = lambda th: np.array([[m.cos(th), -m.sin(th), 0],
                          [m.sin(th), m.cos(th), 0],
                          [0, 0, 1]]).T

T_t = lambda x, y: np.array([[1, 0, x],
                             [0, 1, y],
                             [0, 0, 1]]).T

# Define data set as a rigid transformation of the target set
offSet = targetSet.dot(T_t(0.1, 0.1)).dot(T_r(0.05))

# plt.plot(targetSet[0, :], targetSet[1, :], 'o')
# plt.plot(offSet[0, :], offSet[1, :], 'o')

# Calculate nearest neighbors
nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(targetSet)
distances, indices = nbrs.kneighbors(offSet)
indices = np.ndarray.flatten(indices)

T_est = cv2.estimateRigidTransform(targetSet.tolist(), offSet.tolist(), False)

# Plot results
plt.plot(targetSet[:, 0], targetSet[:, 1], 'o')
plt.plot(offSet[:, 0], offSet[:, 1], 'o')
for i in range(n_pts):
    plt.plot([offSet[i, 0], targetSet[indices[i], 0]], [offSet[i, 1], targetSet[indices[i], 1]], 'k')
    # plt.plot([targetSet[i, 0], offSet[indices[i], 0]], [targetSet[i, 1], offSet[indices[i], 1]], 'k')

plt.ylim((-1.2, 1.2))
plt.xlim((-1.2, 1.2))
plt.show()
