from ICP import *
import numpy as np
import random
import matplotlib.pyplot as plt
import math as m

# n_pts = 100
# S = append_unity(np.random.random((n_pts, 2))*2-1)

# Generate hexgrid of scene points
hex_ratio = m.sqrt(3)/2
grid = [10, m.ceil(10/hex_ratio)]
spacing = [2/(grid[0]-1), 2/(grid[1]-1)]
n_pts_S = grid[0]*grid[1]
S = np.zeros((n_pts_S, 2))
for i in range(grid[0]):
    for j in range(grid[1]):
        S[i*grid[1]+j, 0] = -1 + spacing[1]*j * hex_ratio
        S[i*grid[1]+j, 1] = -1 + spacing[0]*i + 0.5*spacing[0]*(j % 2)
S = append_unity(S)
n_pts = 66

# Define error transformation matrices
T_r = lambda th: np.array([[m.cos(th), -m.sin(th), 0],
                          [m.sin(th), m.cos(th), 0],
                          [0, 0, 1]]).T

T_t = lambda x, y: np.array([[1, 0, x],
                             [0, 1, y],
                             [0, 0, 1]]).T

# Initialize error ranges
shift = 0.05  # units
rotate = 0.1  # radians
noise = 0.05  # units

# Generate specific errors
shift_th = random.random()*2*m.pi
shift_x = shift*m.cos(shift_th)
shift_y = shift*m.sin(shift_th)
rotate_th = random.random()*2*rotate-rotate

# Define data set as a rigid transformation of the target set
T_init = T_t(shift_x, shift_y).dot(T_r(rotate_th))
M = S.dot(T_init)
M[:, :2] += np.random.random((M.shape[0], 2))*2*noise-noise
ix_init = np.random.permutation(n_pts_S)[:n_pts]
M = M[ix_init, :]

# Run ICP on generated data
T_est, indices = ICP(S, M, 15)

# Plot results
plt.plot(S[:, 0], S[:, 1], 'o')
plt.plot(M[:, 0], M[:, 1], 'o')
plt.plot(M.dot(T_est)[:, 0], M.dot(T_est)[:, 1], '.')

M_T = M.dot(T_est)
success = n_pts
for i in range(n_pts):
    plt.plot([M[i, 0], S[ix_init[i], 0]], [M[i, 1], S[ix_init[i], 1]], '0.5')
    plt.plot([M_T[i, 0], S[ix_init[i], 0]], [M_T[i, 1], S[ix_init[i], 1]], '0.5')
    if indices[i] != ix_init[i]:
        plt.plot([M_T[i, 0], S[indices[i], 0]], [M_T[i, 1], S[indices[i], 1]], 'r')
        success -= 1
    # plt.plot([targetSet[i, 0], offSet[indices[i], 0]], [targetSet[i, 1], offSet[indices[i], 1]], 'k')
success = success/n_pts

plt.ylim((-1.5, 1.5))
plt.xlim((-1.5, 1.5))
plt.legend(['Scene', 'Initial Model', 'Final Model'])
plt.show()

print("[shift = (%.03f, %.03f), rotate = %.4f]" % (shift_x, shift_y, rotate_th))
print("Initial error (mse) = ", np.sum((M - S[indices, :])**2)/n_pts)
print("Final error (mse) = ", np.sum((M_T - S[indices, :])**2)/n_pts)
print("Success (%) = ", success*100)
