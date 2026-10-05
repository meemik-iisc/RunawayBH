import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------
# Constants (cgs) -- same values AthenaK uses
# ---------------------------------------------------------------------------
kB     = 1.380649e-16      # erg/K
amu    = 1.660538921e-24   # g   (AthenaK uses the atomic mass unit, not m_p)
pc     = 3.086e18          # cm
s_yr   = 3.15576e7         # s
mu     = 1.0               # AthenaK default when <units>/mu is not set
gamma  = 5.0/3.0
gm1    = gamma - 1.0

script_dir = os.path.dirname(os.path.abspath(__file__))
# ---------------------------------------------------------------------------
# Cooling function
# USE_TABLE = False : exact port of ISMCoolFn() in src/srcterms/ismcooling.hpp
# USE_TABLE = True  : your cooltable.dat (T [K], Lambda [erg cm^3/s]),
#                     interpolated in log-log space
# ---------------------------------------------------------------------------
USE_TABLE = False

_LHD = np.array([
  -22.5977, -21.9689, -21.5972, -21.4615, -21.4789, -21.5497, -21.6211, -21.6595,
  -21.6426, -21.5688, -21.4771, -21.3755, -21.2693, -21.1644, -21.0658, -20.9778,
  -20.8986, -20.8281, -20.7700, -20.7223, -20.6888, -20.6739, -20.6815, -20.7051,
  -20.7229, -20.7208, -20.7058, -20.6896, -20.6797, -20.6749, -20.6709, -20.6748,
  -20.7089, -20.8031, -20.9647, -21.1482, -21.2932, -21.3767, -21.4129, -21.4291,
  -21.4538, -21.5055, -21.5740, -21.6300, -21.6615, -21.6766, -21.6886, -21.7073,
  -21.7304, -21.7491, -21.7607, -21.7701, -21.7877, -21.8243, -21.8875, -21.9738,
  -22.0671, -22.1537, -22.2265, -22.2821, -22.3213, -22.3462, -22.3587, -22.3622,
  -22.3590, -22.3512, -22.3420, -22.3342, -22.3312, -22.3346, -22.3445, -22.3595,
  -22.3780, -22.4007, -22.4289, -22.4625, -22.4995, -22.5353, -22.5659, -22.5895,
  -22.6059, -22.6161, -22.6208, -22.6213, -22.6184, -22.6126, -22.6045, -22.5945,
  -22.5831, -22.5707, -22.5573, -22.5434, -22.5287, -22.5140, -22.4992, -22.4844,
  -22.4695, -22.4543, -22.4392, -22.4237, -22.4087, -22.3928])

def ism_cool_fn(T):
    """Vectorized port of AthenaK's ISMCoolFn (erg cm^3/s)."""
    T = np.asarray(T, dtype=float)
    logt = np.log10(T)
    ki   = 2.0e-19*np.exp(-1.184e5/(T + 1.0e3)) + 2.8e-28*np.sqrt(T)*np.exp(-92.0/T)
    hot  = 10.0**(0.45*logt - 26.065)
    ipps = np.clip((25.0*logt).astype(int) - 103, 0, 100)
    x0   = 4.12 + 0.04*ipps
    dx_  = logt - x0
    spex = 10.0**((_LHD[ipps+1]*dx_ - _LHD[ipps]*(dx_ - 0.04))*25.0)
    return np.where(logt <= 4.2, ki, np.where(logt > 8.15, hot, spex))

if USE_TABLE:
    data = np.loadtxt(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "cooltable.dat"))
    logT_tab, logL_tab = np.log10(data[:, 0]), np.log10(data[:, 1])
    def cool_lambda(T):
        logT = np.clip(np.log10(T), logT_tab.min(), logT_tab.max())
        return 10.0**np.interp(logT, logT_tab, logL_tab)
else:
    cool_lambda = ism_cool_fn

# ---------------------------------------------------------------------------
# Timescales
# ---------------------------------------------------------------------------
def sound_speed(T):
    return np.sqrt(gamma*kB*T/(mu*amu))                  # cm/s

def hydro_dt(T_amb, T_clump, dx):
    return CFL*dx/(np.maximum(sound_speed(T_amb), sound_speed(T_clump)))

def cool_time(n, T):
    """t_cool = e/(n^2 Lambda) = kB T / ((gamma-1) n Lambda)  [s]"""
    return kB*T/(gm1*n*cool_lambda(T))

# ---------------------------------------------------------------------------
# Setup (edit to match your input)
# ---------------------------------------------------------------------------
CFL     = 0.3
Lx      = 8*pc
res_arr = np.array([32,64,128,256])
n_amb, T_amb, T_clump = 0.1, 1e6, 1e4
cool_cfl, cool_nsub_max = 0.5, 100

n_clump   = 1e3
tc_clump  = cool_time(n_clump, T_clump)
tc_amb    = cool_time(n_amb, T_amb)
tcool_min = np.minimum(tc_clump, tc_amb)

print(f"n_clump = {n_clump:.1e}\t tcool_amb = {tc_amb/s_yr:.2e}yr\t tcool_clump = {tc_clump/s_yr:.2e}yr\t tcool_min = {tcool_min/s_yr:.2e}yr")
print(f"{'res'}\t{'dx[pc]'}\t\t{'t_hydro[yr]'}")
for res in res_arr:
    dx = Lx/res
    t_hydro = hydro_dt(T_amb, T_clump, dx)
    print(f"{res}\t{dx/pc:.2e}\t{t_hydro/s_yr:.2e}")

