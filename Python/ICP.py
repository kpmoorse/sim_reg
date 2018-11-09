# ICP_2D: Apply the Iterative Closest Point (ICP) algorithm to a model to generate T_est,
# such that model * T_est ~= scene
#
#   model is an Nx2 matrix or an Nx3 matrix with column 3 all equal to unity;
#       each row corresponds to a point in a 2D point cloud
#
#   scene is an Nx2 matrix or an Nx3 matrix with column 3 all equal to unity;
#       each row corresponds to a point in a 2D point cloud
#
#   num_iters is an integer specifying the number of iterations of ICP to apply
#       before returning T_est
#
# T_est is a 3x3 transform matrix. T_est seems to be always invertible, which would
# make it an affine transform: a combination of translation, rotation, scale, and
# shear transforms.
#
# Kellan Moorse 2018-10-26

import numpy as np
from sklearn.neighbors import NearestNeighbors
import matplotlib.pyplot as plt


def ICP_2D(scene, model, num_iters, W=[]):

    scene = np.array(scene)
    model = np.array(model)

    # Assert input type and format requirements
    assert all([len(x) == len(scene[0]) for x in scene]), "scene must be a 2-dimensional matrix"
    assert scene.shape[1] in [2, 3], "the rows of scene must correspond to 2-dimensional points"
    assert all([len(x) == len(model[0]) for x in model]), "model must be a 2-dimensional matrix"
    assert model.shape[1] in [2, 3], "the rows of model must correspond to 2-dimensional points"
    assert num_iters % 1 == 0, "num_iters must be an integer"

    # Append column of ones to each input if it does not exist
    if scene.shape[1] < 3:
        scene = append_unity(scene)
    else:
        assert all(scene[:, -1] == 1), "All elements in the third column of scene must be unity"
    if model.shape[1] < 3:
        model = append_unity(model)
    else:
        assert all(model[:, -1] == 1), "All elements in the third column of model must be unity"

    # Initialize transform matrix
    T_est = np.identity(3)

    for i in range(num_iters):

        model_current = model.dot(T_est)

        # Calculate nearest neighbors
        nbrs = NearestNeighbors(n_neighbors=1, algorithm='auto').fit(scene)
        distances, indices = nbrs.kneighbors(model_current)
        indices = np.ndarray.flatten(indices)

        # Approximate transform based on nearest neighbors
        T_est_new = T_est.dot(linreg(model_current, scene[indices, :], W))

        # Break if there is no change in T_est
        if (T_est_new == T_est).all:
            break
        else:
            T_est = T_est_new

        # plt.plot(scene[:, 0], scene[:, 1], 'o')
        # plt.plot(model[:, 0], model[:, 1], 'o')
        # plt.plot(model.dot(T_est)[:, 0], model.dot(T_est)[:, 1], '.')
        # plt.show()
        # input()

    return T_est, indices


def ICP_3D(scene, model, num_iters, W=[]):

    # Assert input type and format requirements
    assert all([len(x) == len(scene[0]) for x in scene]), "scene must be a 2-dimensional matrix"
    assert scene.shape[1] in [3, 4], "the rows of scene must correspond to 2-dimensional points"
    assert all([len(x) == len(model[0]) for x in model]), "model must be a 2-dimensional matrix"
    assert model.shape[1] in [3, 4], "the rows of model must correspond to 2-dimensional points"
    assert num_iters % 1 == 0, "num_iters must be an integer"

    scene = np.array(scene)
    model = np.array(model)

    # Append column of ones to each input if it does not exist
    if scene.shape[1] < 3:
        scene = append_unity(scene)
    else:
        assert all(scene[:, -1] == 1), "All elements in the third column of scene must be unity"
    if model.shape[1] < 3:
        model = append_unity(model)
    else:
        assert all(model[:, -1] == 1), "All elements in the third column of model must be unity"


# Append a column (or row) of ones to a matrix
def append_unity(matrix, ax=1):

    return np.concatenate((matrix, np.ones((np.shape(matrix)[0], 1))), axis=ax)

def linreg(X, Y, W=[]):

    if len(W) == 0:
        W = np.identity(X.shape[0])

    return np.linalg.inv(X.T.dot(W).dot(X)).dot(X.T).dot(W).dot(Y)