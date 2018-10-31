import numpy as np
import matplotlib.pyplot as plt


def prepend_unity(matrix, ax=1):

    return np.concatenate((np.ones((np.shape(matrix)[0], 1)), matrix), axis=ax)


f = lambda x: np.sin(np.pi*x)
K = lambda x, h: (1-(np.abs(x)/h)**3)**3

N = 50

xx = prepend_unity(np.arange(-1, 1.01, 0.01).reshape((-1, 1)))
X = prepend_unity(np.random.random((N, 1))*2-1)
y = f(X[:, 1:2]) + np.random.normal(0, 0.1, (N, 1))

B_ols = np.linalg.inv(X.T.dot(X)).dot(X.T).dot(y)

W = np.identity(N)
B_wls0 = np.linalg.inv(X.T.dot(W).dot(X)).dot(X.T).dot(W).dot(y)

W = np.diag(K(X[:, 1], 1))
B_wls1 = np.linalg.inv(X.T.dot(W).dot(X)).dot(X.T).dot(W).dot(y)


plt.plot(xx[:, 1], f(xx[:, 1]))
plt.plot(X[:, 1:2], y, 'o')
plt.plot(xx[:, 1], xx.dot(B_ols))
plt.plot(xx[:, 1], xx.dot(B_wls0), '--')
plt.plot(xx[:, 1], xx.dot(B_wls1))
plt.ylim((-1.5, 1.5))

plt.legend(('True function', 'Noisy Data', 'OLS Fit', 'WLS Fit (uniform)', 'WLS Fit (center-weighted)'))

plt.show()
