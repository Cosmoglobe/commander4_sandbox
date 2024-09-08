import numpy as np
import matplotlib.pyplot as plt


def inversion_sampler(lnlike, data, sigma, nu, theta_min, theta_max, num_steps):
    theta_arr = np.linspace(theta_min, theta_max, num_steps)
    ln_prob = np.zeros(num_steps)
    for i in range(num_steps):
        ln_prob[i] = lnlike(theta_arr[i], data, sigma, nu)
    ln_prob -= ln_prob.max()

    prob = np.exp(ln_prob)
    F = np.cumsum(prob)/prob.sum()
    eta = np.random.uniform()
    theta_hat = theta_arr[np.argmin(abs(eta -F))]
    return theta_hat

def signal(beta, nu):
    return 100*(nu/30)**beta

def lnlike(beta, d, sigma, nu):
    return -0.5*(((signal(beta, nu) - d)/sigma)**2).sum()

freq = np.array([30,44,70,100,353])
sigma = np.array([20, 10, 5, 8, 30])/2

d = signal(-3.0, freq) + np.random.randn(len(freq))*sigma

plt.figure()
plt.errorbar(freq, d, yerr=sigma)


freq_tot = np.geomspace(30, 353)

beta_min = -10
beta_max = -1
num_steps = 1000

num_samps = 1000
beta_samps = np.zeros(num_samps)
for i in range(num_samps):
    betahat = inversion_sampler(lnlike, d, sigma, freq, \
            beta_min, beta_max, num_steps)
    if i < 10:
        plt.plot(freq_tot, signal(betahat, freq_tot))
        print(betahat)
    beta_samps[i] = betahat
plt.xscale('log')
plt.figure()
plt.hist(beta_samps)
plt.show()
