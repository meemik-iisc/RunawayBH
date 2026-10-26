#To solve SDE dX = Lambda*X(t)dt + mu*X(t)*dW(t)

import numpy as np
import matplotlib.pyplot as plt
import os

#Parameters
Lambda  = 2
mu      = 1
X0      = 1
delt    = 2**(-8)
R       = 4
dt      = R*delt
N       = int(1/delt)
N_em    = int(1/dt)
#Generate Brownian path
rng     = np.random.default_rng(seed=10)
dW      = np.sqrt(delt)*rng.standard_normal(N)
W       = np.concatenate([[0.0],np.cumsum(dW)])

#Integrate over the path
X = np.zeros(N_em+1)
X[0]=X0
for i in range(N_em):
    dWj = 0
    for j in range(i*R, (i+1)*R):
       dWj+=dW[j]
    X[i+1]=X[i]+Lambda*X[i]*dt+mu*X[i]*dWj 

#Calculate analytic solution
t       = np.linspace(0,1,N+1)
Xtrue   = X0*np.exp((Lambda-0.5*mu**2)*t+mu*W)
t_em    = np.linspace(0,1,N_em+1)

#Plot
plt.plot(t_em, X, 'r--*', label = "EM method")
plt.plot(t, Xtrue, 'b', label="Exact")
plt.legend()
script_dir  = os.path.dirname(os.path.abspath(__file__))
filename    = os.path.join(script_dir,"euler_murayama.png")
plt.savefig(filename)