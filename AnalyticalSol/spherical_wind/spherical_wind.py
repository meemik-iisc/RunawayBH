import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os
#Constants
G       = 6.67e-8       #cm^3/g/s^2
kB      = 1.38e-16      #erg/K
mp      = 1.67e-24      #g
kpc     = 3.086e21      #cm
pc      = 3.086e18      #cm
Msun    = 1.989e33      #g
mu      = 1.0
s_yr    = 3.154e7       #s
s_Myr   = 1e6*s_yr
# K       = 2.91e28       #g^(-2/3)cm^4s^(-2)
gamma   = 5.0/3.0       #Polytropic index
gm1     = gamma-1
kms     = 1e5

#Parameters
v_bh    = 1000*1e5                  #cm/s
M_bh    = 2e7*Msun                  #g
R_bondi = 2*G*M_bh/(v_bh**2)        #cm
epsilon = 0.1*R_bondi               #cm
# rho_vir = mp                        #g/cm^3
n_cgm   = 1                      #g/cm^3
rho_cgm = n_cgm*mp*mu
T_cgm   = 1e6                       #K
cs_cgm  = np.sqrt(kB*T_cgm/(mu*mp))
P_cgm   = n_cgm*kB*T_cgm

#Outflow Parameters
v_wind  = 1e2*kms       #cm/s
R0      = 0.25*kpc
t_sys   = 70*1e6*s_yr   #s
eta     = 0.1
M_dot_w = ((v_bh**2)*4*np.pi*(R0**2)*rho_cgm)/(v_wind)
# M_dot_w = 2*Msun/s_yr
print(f"Wind Mass Flow Rate = {M_dot_w*s_yr/Msun:.2f} Msun/yr")
Rinj    = 50*pc

#Code Units
L_code  = 1.0*kpc        #cm
# L_code  = 1.0*pc         #cm
v_code  = 1e3*kms            #cm/s
rho_code= 1.0*mp         #g/cm^3

M_code  = rho_code*L_code**3
t_code  = L_code/v_code
P_code  = rho_code*v_code**2
T_code  = (mu*mp/kB)*(v_code**2)

# ============================================================================
# COOLING TABLE
# ============================================================================
data = np.loadtxt('cooltable.dat')
T_tab = data[:, 0]
Lambda_tab = data[:, 1]
Lambda_interp = interp1d(T_tab, Lambda_tab, kind='linear', 
                         bounds_error=False, fill_value='extrapolate')

def cool_lambda(T):
    """Safe cooling function with clamping"""
    T_safe = np.clip(T, T_tab.min() * 0.95, T_tab.max() * 1.05)
    return Lambda_interp(T_safe)


K = (np.sqrt(4*np.pi)/(2**(0.2)))*np.pow(rho_cgm*v_wind**3/M_dot_w, 0.3)
print(f"Kappa = {K:.2e}")

#Solve IVP
def rhs(t,r_TS):
    return v_wind-(K*r_TS*(t**(-0.4)))

#Initial condition
t0 = np.sqrt((M_dot_w)/(rho_cgm*(v_wind**3)))
# t0 = 1*s_Myr
print(f"t01 = {t0/s_Myr:.2f} Myr")
t1 = 1/(K*v_wind)
print(f"t1 = {t1/s_Myr:.2e} Myr")
r_TS0 = (v_wind*t0**(0.4))/K
print(f"r_TS(t0) = {r_TS0/kpc:.2f} kpc")
t_end = 1e2*t0
t_eval = np.linspace(t0, t_end, 1000)

# Solve the IVP
sol = solve_ivp(rhs, (t0, t_end), [r_TS0], t_eval=t_eval, method='RK45',
                rtol=1e-6, atol=1e-9)

# Check that the solver succeeded
if not sol.success:
    print("Solver failed:", sol.message)
r_TS = sol.y[0]
t_arr = sol.t

r_OS = ((M_dot_w*(v_wind**2)*(t_arr**3))/(rho_cgm))**(0.2)

output_dir = os.path.dirname(os.path.abspath(__file__))
# Plot the numerical solution
fig1,ax1 = plt.subplots(figsize=(8,6))
# ax1.plot(np.log10(t_arr/t0), np.log10(r_TS/kpc), color = 'red', label = r'$r_{TS}$')
# ax1.plot(np.log10(t_arr/t0), np.log10(r_OS/kpc), color = 'blue', label = r'$r_{OS}$')
ax1.plot(t_arr/t0, r_TS/kpc, color = 'red', label = r'$r_{TS}$')
ax1.plot(t_arr/t0, r_OS/kpc, color = 'blue', label = r'$r_{OS}$')
ax1.set_xlabel(r'$t/t_0$')
ax1.set_ylabel('r (kpc)')
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.grid(True, which='both', linestyle='-', alpha=0.8)
ax1.minorticks_on()
ax1.set_title(f"r vs t for $t_0$ = {t0/s_Myr:.0f} Myr")
plt.tight_layout()
plt.legend(loc='best')
save_filename1 = os.path.join(output_dir,f"r_vs_t.png")
plt.savefig(save_filename1)
plt.close()

# Compute v_TS(t) = dr_TS/dt using the ODE and the solution r_TS(t)


v_TS = v_wind - K * r_TS * (t_arr**(-0.4))
v_OS =  ((M_dot_w*(v_wind**2))/(rho_cgm*(t_arr**2)))**(0.2)

# Plot v_TS versus time, in useful units (e.g., km/s)
fig2, ax2 = plt.subplots(figsize=(8,6))
ax2.plot(t_arr/t0, v_TS/kms, color = 'red', label = r'$v_{TS}$')
ax2.plot(t_arr/t0, v_OS/kms, color = 'blue', label = r'$v_{OS}$')
# ax2.plot(np.log10(t_arr/t0), v_TS/kms, color = 'red', label = r'$v_{TS}$')
# ax2.plot(np.log10(t_arr/t0), v_OS/kms, color = 'blue', label = r'$v_{OS}$')
ax2.set_xlabel(r'$t/t_0$')
ax2.set_ylabel('v (km/s)')
ax2.set_xscale('log')
ax2.grid(True, which='both', linestyle='-', alpha=0.8)
ax2.minorticks_on()
ax2.set_title(f"v vs t for $t_0$ = {t0/s_Myr:.0f} Myr")
plt.legend(loc='best')
plt.tight_layout()

save_filename2 = os.path.join(output_dir, "v_vs_t.png")
plt.savefig(save_filename2)
plt.close()



# Plot v_TS versus time, in useful units (e.g., km/