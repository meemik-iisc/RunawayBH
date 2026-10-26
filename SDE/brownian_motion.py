import numpy as np
import matplotlib.pyplot as plt
import os
T=1
N=1000
dt = T/N

rng = np.random.default_rng(seed=40)

W=np.zeros(N+1)
t=np.zeros(N+1)
for i in range(N):
    x = rng.standard_normal()
    W[i+1]=W[i]+np.sqrt(dt)*x
    t[i+1]=t[i]+dt 

script_dir  = os.path.dirname(os.path.abspath(__file__))
filename    = os.path.join(script_dir,"brownian_motion.png")
plt.plot(t,W)
plt.savefig(filename)