import numpy as np
from scipy.interpolate import interp1d
import os
import matplotlib.pyplot as plt
#Constants
G       = 6.67e-8       #cm^3/g/s^2
kB      = 1.38e-16      #erg/K
mp      = 1.67e-24      #g
kpc     = 3.086e21      #cm
pc      = 3.086e18      #cm
Msun    = 1.989e33      #g
mu      = 1.0
s_yr    = 3.154e7       #s
# K       = 2.91e28       #g^(-2/3)cm^4s^(-2)
gamma   = 5.0/3.0       #Polytropic index
gm1     = gamma-1
kms     = 1e5

# ============================================================================
# COOLING TABLE
# ============================================================================
script_dir = os.path.dirname(os.path.abspath(__file__))
cool_dir = os.path.join(script_dir,"cooltable.dat")
data = np.loadtxt(cool_dir)
T_tab = data[:, 0]
Lambda_tab = data[:, 1]
Lambda_interp = interp1d(T_tab, Lambda_tab, kind='linear', 
                         bounds_error=False, fill_value='extrapolate')

def cool_lambda(T):
    """Safe cooling function with clamping"""
    T_safe = np.clip(T, T_tab.min() * 0.95, T_tab.max() * 1.05)
    return Lambda_interp(T_safe)

def hydro_dt(Tamb, Tclump):
    cs_amb      = np.sqrt((kB*T_amb)/(mu*mp))    #cm/s
    cs_clump    = np.sqrt((kB*T_clump)/(mu*mp))
    return      CFL*dx/np.max(cs_amb,cs_clump)
def cool_dt(n, T):
    num = gamma*kB*T
    den = gm1*n*cool_lambda(T)
    return num/den
CFL = 0.3
Lx  = 8*pc             #cm
nx  = 64
dx  = Lx/nx

n_clump     = np.logspace(-2,3,1000)
n_amb       = 0.1
T_amb       = 1e6                #K
T_clump     = 1e4                #K
cool_dt_amb = cool_dt(n_amb, T_amb)
for n in n_clump:
    cool_dt_clump = cool_dt(n, T_clump)
    tcool_min = np.min(cool_dt_clump,cool_dt_amb)
    print(f"n_clump = {n:>8g}cm^-3")
# print(f"For n = {n:.1e} cm^-3: ")
# print(f"T = 1e6K    : dt_hydro = {hydro_dt(T1)/s_yr:.2e}yr    : dt_cool = {cool_dt(n, T1)/s_yr:.2e}yr")
# print(f"T = 1e4K    : dt_hydro = {hydro_dt(T2)/s_yr:.2e}yr    : dt_cool = {cool_dt(n, T2)/s_yr:.2e}yr")
