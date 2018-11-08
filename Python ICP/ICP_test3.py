#Test ICP algorithm with centroid-distance-squared (cd-square) weighting vs uniform

from ICP import *
import numpy as np
import random
import matplotlib.pyplot as plt
import math as m
from tqdm import tqdm

shift_range = np.arange(0, 0.1, 0.01)
n_samples = 100

ols_dist = []
wls_dist = []

# Loop over many shift values
for shift in tqdm(shift_range):

    # Loop over many samples
    ols_sample = []
    wls_sample = []
    for _ in range(n_samples):

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
        n_pts = m.floor(n_pts_S*0.67)

        # Define error transformation matrices
        T_r = lambda th: np.array([[m.cos(th), -m.sin(th), 0],
                                  [m.sin(th), m.cos(th), 0],
                                  [0, 0, 1]]).T

        T_t = lambda x, y: np.array([[1, 0, x],
                                     [0, 1, y],
                                     [0, 0, 1]]).T

        # Initialize error ranges
        # shift = 0.05  # units
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

        # Calculate weights based on distance from centroid
        # (More distant points are weighted more heavily)
        N = M.shape[0]
        m_bar = np.mean(M, axis=0)
        w_list = np.linalg.norm(M - np.tile(m_bar, (N, 1)), axis=1)
        # w_list = np.random.random(N)
        W = np.diag(w_list)

        # Run ICP on generated data
        T_ols, indices = ICP_2D(S, M, 15)
        success_ols = sum(ix_init == indices) / n_pts
        ols_sample.append(success_ols)
        T_wls, indices = ICP_2D(S, M, 15, W=W)
        success_wls = sum(ix_init == indices) / n_pts
        wls_sample.append(success_wls)

    ols_dist.append(sum(np.array(ols_sample)==1)/n_samples)
    wls_dist.append(sum(np.array(wls_sample) == 1) / n_samples)

# bins = np.linspace(0, 1, 25)
# plt.hist(ols_sample, bins, alpha=0.5, edgecolor='k')
# plt.hist(wls_sample, bins, alpha=0.5, edgecolor='k')
# plt.title('ICP proportions of successful correspondence')
# plt.legend(('OLS', 'WLS'))
# plt.show()
plt.title('ICP proportions of successful correspondence')
plt.plot(shift_range, ols_dist, shift_range, wls_dist)
plt.xlabel('Shift Magnitude')
plt.legend(('OLS', 'WLS'))
plt.show()
