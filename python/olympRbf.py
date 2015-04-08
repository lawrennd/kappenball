import scipy.io
import numpy as np
import matplotlib.pyplot as plt

olympicData = scipy.io.loadmat('/home/neil/public_html/olympics.mat')['male100']

years = olympicData[:, 0]
times = olympicData[:, 1]

t = times[:, None]

for numBasis in range(27, 28):
    mu = np.r_[1896:2012:(numBasis*1j)]
    Phi = np.zeros((years.shape[0], numBasis))
    # Set the width of the basis functions
    width = 2*(years.max() - years.min())/numBasis
    # precompute the scale of the exponential
    scale = -1./(2.*width*width)
    for i in range(0, numBasis):
        Phi[:, i] = np.exp(scale*(years-mu[i])**2)
    
    PhiTPhi = np.dot(Phi.T, Phi)
    PhiTt = np.dot(Phi.T, t)

    #wStar = np.dot(np.linalg.inv(PhiTPhi), PhiTt)
    wStar = np.linalg.solve(PhiTPhi, PhiTt)

    yearsPred = np.r_[1896:2012:1]
    PhiPred = np.zeros((yearsPred.shape[0], numBasis))
    for i in range(0, numBasis):
        PhiPred[:, i] = np.exp(width*(yearsPred-mu[i])**2)

    tPred = np.dot(PhiPred, wStar)

    plt.figure()
    plt.axis([1892, 2012, 9, 13])
    plt.plot(yearsPred, tPred)
    plt.plot(years, t, 'rx')
    plt.show()
