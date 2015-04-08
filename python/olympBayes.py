import scipy.io
import numpy as np
import matplotlib.pyplot as plt

olympic_data = scipy.io.loadmat('/home/neil/public_html/olympics.mat')['male100']
years = olympic_data[:, 0]
times = olympic_data[:, 1]
num_data = years.shape[0]

alpha = 100
sigma2 = 0.09

t = times[:, None]
log_marginal = []
#basis_range = range(1, 28)
basis_span = (1896, 2012)
basis_range = range(5, 100)
#basis_span = (1850, 2050)
add_bias = False
add_linear = False

if add_bias:
    alpha = 2
for num_basis in basis_range:
    mu = np.linspace(basis_span[0], basis_span[1], num_basis)

    # Set the width of the basis functions
    width = 2*(years.max() - years.min())/num_basis
    # precompute the scale of the exponential
    scale = -1./(2.*width*width)
    Phi = np.zeros((years.shape[0], num_basis))
    for i in range(0, num_basis):
        Phi[:, i] = np.exp(scale*(years-mu[i])**2)

    if add_bias:
        Phi = np.hstack((10*np.ones((num_data, 1)), Phi))

    if add_linear:
        Phi = np.hstack(((years[:, None]-1950)/100, Phi))

    PhiPhiT = np.dot(Phi, Phi.T)
    K = alpha*PhiPhiT + np.eye(num_data)*sigma2

    log_marginal.append(-num_data/2.*np.log(2*np.pi) - 0.5*np.log(np.linalg.det(K)) - 0.5*np.dot(np.dot(t.T, K), t))
    

plt.close('all')
plt.plot(np.array(log_marginal).flatten())
plt.show()
best_basis = basis_range[np.argmax(log_marginal)]
print("Best basis contains " + str(best_basis) + " basis functions")
raw_input("Press return key\n")
plt.close()

Phi = np.zeros((years.shape[0], num_basis))
mu = np.linspace(basis_span[0], basis_span[1], best_basis)

# Set the width of the basis function
width = 2*(years.max() - years.min())/best_basis
scale = -1./(2.*width*width)
for i in range(0, best_basis):
    Phi[:, i] = np.exp(scale*(years-mu[i])**2)


if add_bias:
    Phi = np.hstack((10*np.ones((num_data, 1)), Phi))

if add_linear:
    Phi = np.hstack(((years[:, None]-1950)/100, Phi))

PhiTPhi = np.dot(Phi.T, Phi)
PhiTt = np.dot(Phi.T, t)

C_w = np.linalg.inv(1./sigma2*PhiTPhi + 1./alpha*np.eye(Phi.shape[1]))
mu_w = np.dot(C_w/sigma2, PhiTt)


num_pred = 201
years_pred = np.linspace(1850, 2050, num_pred)
Phi_pred = np.zeros((years_pred.shape[0], num_basis))
for i in range(0, num_basis):
    Phi_pred[:, i] = np.exp(scale*(years_pred-mu[i])**2)

if add_bias:
    Phi_pred = np.hstack((10*np.ones((num_pred, 1)), Phi_pred))

if add_linear:
    Phi_pred = np.hstack(((years_pred[:, None]-1950)/100, Phi_pred))

# Compute mean and one standard deviation
t_pred = np.dot(Phi_pred, mu_w).flatten()
t_var = (np.dot(Phi_pred, C_w)*Phi_pred).sum(axis=1)
t_std = np.sqrt(t_var)

plt.plot(years_pred, t_pred,color="#204a87",linewidth=2)
plt.fill(np.hstack((years_pred,years_pred[::-1])),np.hstack((t_pred+t_std,t_pred[::-1]-t_std[::-1])),color="#729fcf",linewidth=0.5,alpha=0.3)
plt.plot(years, t, 'kx',mew=1.5)
plt.axis([1850, 2050, 0, 13])
plt.show()
raw_input("Press return key\n")
plt.close()

w_samps = np.random.multivariate_normal(mu_w.flatten(), C_w, 10)
t_samp = np.dot(Phi_pred, w_samps.T)
plt.plot(years_pred, t_samp,linewidth=2)
plt.plot(years, t, 'kx',mew=1.5)
plt.axis([1850, 2050, 0, 13])
plt.show()
raw_input("Press return key\n")
plt.close()
