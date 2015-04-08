import scipy.io
import numpy as np
import matplotlib.pyplot as plt

olympicData = scipy.io.loadmat('/home/neil/public_html/olympics.mat')['male100']

years = olympicData[:, 0]
times = olympicData[:, 1]

t = times[:, None]
stdt = np.sqrt(np.var(t))
meant = np.mean(t)
t = (t-meant)/stdt
for order in range(1, 27):
    yearsSpan = years[-1]-years[0]
    Phi = np.zeros((years.shape[0], order+1))
    for i in range(0, order+1):
        Phi[:, i] = (1*(years-years.min())/yearsSpan-.5)**i
    
    PhiTPhi = np.dot(Phi.T, Phi)
    PhiTt = np.dot(Phi.T, t)

    wStar = np.linalg.solve(PhiTPhi, PhiTt)

    yearsPred = np.r_[1896:2012:1]
    PhiPred = np.zeros((yearsPred.shape[0], order+1))
    for i in range(0, order+1):
        PhiPred[:, i] = (1*(yearsPred-years.min())/yearsSpan-.5)**i

    tPred = np.dot(PhiPred, wStar)

    plt.axis([1892, 2012, 9, 13])
    plt.plot(yearsPred, tPred*stdt + meant)
    plt.plot(years, t*stdt+meant, 'rx')
    plt.show()
