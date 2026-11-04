#To solve SDE dX = -(1/\tau)*X(t)dt + sqrt(c)*dW(t)

import numpy as np
import matplotlib.pyplot as plt
import os

#Parameters
tau     = 10                #Relaxation time
T       = 10*tau
c       = 1                 #Diffusion Constant
X0      = 1
N       = 1000
dt      = T/N
#Generate Brownian path
rng     = np.random.default_rng(seed=10)
dW      = np.sqrt(dt)*rng.standard_normal(N)
W       = np.concatenate([[0.0],np.cumsum(dW)])

#Integrate over the path
X = np.zeros(N+1)
X[0]=X0
for i in range(N):
    X[i+1]=X[i]-(X[i]/tau)*dt+np.sqrt(c)*dW[i] 

#Calculate analytic solution
t       = np.linspace(0,T,N+1)
I       = np.concatenate([[0.0], np.cumsum(np.exp(t[:-1]/tau)*dW)])   # \int e^{s/τ} dW
Xtrue   = np.exp(-t/tau)*(X0+np.sqrt(c)*I)
# t_em    = np.linspace(0,1,N_em+1)

#Plot
plt.plot(t, X, 'r--*', label = "Numerical")
plt.plot(t, Xtrue, 'b', label="Exact")
plt.xlabel("t")
plt.ylabel("X(t)")
plt.legend()
script_dir  = os.path.dirname(os.path.abspath(__file__))
filename    = os.path.join(script_dir,"ou_process.png")
plt.savefig(filename)