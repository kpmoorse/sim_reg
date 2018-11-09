import numpy as np
import matplotlib.pyplot as plt
import math as m
from sklearn.neighbors import NearestNeighbors
import random

# Generate ideal target set
n_pts = 25
targetSet = np.concatenate((np.random.random((n_pts, 2))*2-1, np.ones((n_pts, 1))), axis=1)

#Initialize transform matrix
T_est = np.identity(3)

# Define error transformation matrices
T_r = lambda th: np.array([[m.cos(th), -m.sin(th), 0],
                          [m.sin(th), m.cos(th), 0],
                          [0, 0, 1]]).T

T_t = lambda x, y: np.array([[1, 0, x],
                             [0, 1, y],
                             [0, 0, 1]]).T

#Initialize error ranges
shift = 0.2 #units
rotate = 0.1 #radians
noise = 0.05 #units

# Define data set as a rigid transformation of the target set
T_init = T_t(random.random()*2*shift-shift, random.random()*2*shift-shift).dot(
    T_r(random.random()*2*rotate-rotate)
)
offSet = targetSet.dot(T_init)
offSet[:, :2] += np.random.random((n_pts, 2))*2*noise-noise
np.random.shuffle(offSet)

# plt.plot(targetSet[0, :], targetSet[1, :], 'o')
# plt.plot(offSet[0, :], offSet[1, :], 'o')

for i in range(10):
    offSet_current = offSet.dot(T_est)

    # Calculate nearest neighbors
    nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(targetSet)
    distances, indices = nbrs.kneighbors(offSet_current)
    indices = np.ndarray.flatten(indices)

    # T_est = cv2.estimateRigidTransform(targetSet.tolist(), offSet.tolist(), False)
    T_est = T_est.dot(np.linalg.pinv(offSet_current).dot(targetSet[indices, :])) #NOTE: transform is NOT affine

# Plot results
plt.plot(targetSet[:, 0], targetSet[:, 1], 'o')
plt.plot(offSet[:, 0], offSet[:, 1], 'o')
plt.plot(offSet.dot(T_est)[:, 0], offSet.dot(T_est)[:,1], '.')
# for i in range(n_pts):
#     plt.plot([offSet[i, 0], targetSet[indices[i], 0]], [offSet[i, 1], targetSet[indices[i], 1]], 'k')
#     # plt.plot([targetSet[i, 0], offSet[indices[i], 0]], [targetSet[i, 1], offSet[indices[i], 1]], 'k')
plt.ylim((-1.2, 1.2))
plt.xlim((-1.2, 1.2))
plt.legend(['Scene', 'Initial Model', 'Final Model'])
plt.show()
