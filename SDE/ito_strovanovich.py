import numpy as np
import matplotlib.pyplot as plt
import os
T   =1
N   =1000
dt  = T/N

rng = np.random.default_rng(seed=40)

t   = np.linspace(0,T,N+1)
dW  = np.sqrt(dt)*rng.standard_normal(N)
W   = np.concatenate([[0.0], np.cumsum(dW)])

ito = 0
stro = 0
ito_analytic = 0.5*(W[-1]**2-T)
stro_analytic = 0.5*W[-1]**2
for i in range(N):
    ito+=W[i]*dW[i]
    stro+= 0.5*(W[i+1]+W[i])*dW[i]

print(f"Ito             = {ito:.6f},    Analytic = {ito_analytic:.6f},  Error = {np.abs(ito-ito_analytic):.6f}")
print(f"Strovanovich    = {stro:.6f},   Analytic = {stro_analytic:.6f}, Error = {np.abs(stro-stro_analytic):.6f}")
print(f"Difference      = {np.abs(ito-stro):.6f}")
