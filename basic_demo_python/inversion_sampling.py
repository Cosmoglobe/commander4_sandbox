import numpy as np
import matplotlib.pyplot as plt


freq = np.array([30,44,70,100,353])
sigma = np.array([20, 10, 5, 8, 30])

def signal(beta):
    return 100*(freq/30)**beta


d = 100*(freq/30)**-3 + np.random.randn(len(freq))*sigma



beta = np.linspace(-10,-1,1000)
chi2 = np.zeros_like(beta)
for i in range(len(beta)):
    chi2[i] = (((signal(beta[i]) - d)/sigma)**2).sum()
chi2 -= chi2.min()
plt.plot(beta, chi2)
plt.figure()
P = np.exp(-chi2/2)/np.sum(np.exp(-chi2/2))
plt.plot(beta, P)

plt.figure()
plt.plot(beta, np.cumsum(P))

F = np.cumsum(P)


plt.figure()
plt.errorbar(freq, d, yerr=sigma)
freq = np.geomspace(30, 353)
for i in range(10):
    eta = np.random.uniform()
    betahat = beta[F > eta][0]
    plt.plot(freq, signal(betahat))
plt.xscale('log')
plt.show()
