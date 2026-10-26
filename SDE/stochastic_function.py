import numpy as np
import matplotlib.pyplot as plt
import os
T   =1
N   =1000
M   = 100
dt  = T/N

rng = np.random.default_rng(seed=10)

t   = np.linspace(0,T,N+1)
dW  = np.sqrt(dt)*rng.standard_normal((M,N))
W   = np.cumsum(dW, axis=1)
W   = np.hstack([np.zeros((M,1)),W])
tM  = np.tile(t,(M,1))
u   = np.exp(tM+(0.5)*W)

plt.plot(t, u[:5].T,'r', lw=0.5)
plt.plot(t, u.mean(axis=0),'b', label="mean")
plt.plot(t, np.exp(9*t/8), 'b', linestyle="dashed", label = "theoretical mean")
plt.legend()
script_dir  = os.path.dirname(os.path.abspath(__file__))
filename    = os.path.join(script_dir,"stochastic_function.png")
plt.savefig(filename)