import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
script_dir = os.path.dirname(os.path.abspath(__file__))

#Constants
G       = 6.67e-8       #cm^3/g/s^2
kB      = 1.38e-16      #erg/K
mp      = 1.67e-24      #g
pc      = 3.086e18      #cm
Msun    = 1.989e33      #g
mu      = 1.0
gamma   = 5.0/3.0       #Polytropic index
gm1     = gamma-1
s_yr    = 3.154e7       #s

#Parameters
v_bh    = 1000*1e5                  #cm/s
M_bh    = 2e7*Msun                  #g
R_bondi = 2*G*M_bh/(v_bh**2)        #cm
print(f"Bondi radius: {R_bondi/pc:.2e} pc")
epsilon = 10*pc                   #cm

#Calculate Virial Temp at 1kpc
r_vir   = 1e3*pc        #cm          
T_vir   = G*M_bh*mp/(r_vir*kB)  #K
print(f"Virial Temperature at 1kpc = {T_vir:.2e}K")

#Calculate Density if pressure equilibrium with CGM at 1kpc
rho_cgm = 1e-3*mp       #g/cm^3
T_cgm   = 1e6           #K
rho_vir = rho_cgm*T_cgm/T_vir
print(f"Virial density at 1kpc = {rho_vir/mp:.2e} mp/cc")

#Calculate Polytropic Constant
K = (kB*T_vir/(mu*mp))*(rho_vir**(1-gamma))
print(f"Polytropic Constant = {K:.2e}")

#Assume: T(r=0) = 1e7K, calculate rho(r=0) assuming polytropic equilibrium
T0      = 1e7       #K
rho0    = rho_vir*((T0/T_vir)**(1/(gm1)))       #g/cm^3
print(f"Density at r=0 = {rho0/mp:.2e} mp/cc")

 # ============================================================================
# COOLING TABLE
# ============================================================================

data = np.loadtxt(os.path.join(script_dir, 'cooltable.dat'))
T_tab = data[:, 0]
Lambda_tab = data[:, 1]
Lambda_interp = interp1d(T_tab, Lambda_tab, kind='linear', 
                        bounds_error=False, fill_value='extrapolate')




#Gravitational Potential
def Phi(r):
    return -1*G*M_bh/np.sqrt(r**2+epsilon**2)
    # return -1*G*M_bh/r
    

def rho(r):
    # term1 = (-1*(gm1/(K*gamma))*Phi(r))
    # term2 = (-1*(gm1/(K*gamma))*Phi(0))
    # if (term1 < 0).any() or (term2 < 0).any():
    #     print("term1 =", term1)
    #     print("term2 =", term2)
    # print((term2**(1/gm1))/mp)
    # return rho0+term1**(1/gm1)-term2**(1/gm1)
    term1 = (-1*(gm1/(K*gamma))*Phi(r))
    term2 = (-1*(gm1/(K*gamma))*Phi(r_vir))
    # print(term2**(1/gm1)/mp)
    print("term1=",(term1**(1/gm1))/mp)
    if (term1 < 0).any() or (term2 < 0).any():
        print("term1 =", term1)
        print("term2 =", term2)
    print((term2**(1/gm1))/mp)
    return rho_vir+term1**(1/gm1)-term2**(1/gm1)
def pressure(r):
    return K*(rho(r)**gamma)

def Temperature(r):
    return (K*mu*mp/kB)*(rho(r)**(gm1))

def entropy(r):
    return pressure(r)/(rho(r)**gamma)

def cool_lambda(T):
    return Lambda_interp(T)

def cooling_time(r):
    n = rho(r)/(mu*mp)
    return (gamma*pressure(r))/(gm1*(n**2)*cool_lambda(Temperature(r)))

def free_fall_time(r):
    return np.sqrt(np.pow((r**2+epsilon**2),1.5)/(2*G*M_bh))

#Plot profiles
r = np.linspace(1e-2*pc,1e3*pc,10000)
print(rho(r))

# Create 3x2 subplots
fig, axes = plt.subplots(2, 3, figsize=(16, 10))
fig.suptitle(f"Adiabatic Equilibrium around Black Hole for $\gamma$={gamma:.3f},  $\epsilon$={epsilon/pc:.2f} pc")

# Top-left: Density
axes[0,0].plot(np.log10(r/pc), np.log10(rho(r)/mp), 'b-', label='Density')
axes[0,0].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[0,0].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length ($\epsilon$)')
axes[0,0].set_xlabel('Radius (pc)')
axes[0,0].set_ylabel('Density ($m_p/cm^3$) [log]')
axes[0,0].set_title('Density Profile')
# axes[0,0].set_xscale('log')
# axes[0,0].set_yscale('log')
axes[0,0].legend()
axes[0,0].grid(True, alpha=0.3)

# Top-Middle: Pressure  
axes[0,1].plot(np.log10(r/pc), np.log10(pressure(r)), 'g-', label='Pressure')
axes[0,1].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[0,1].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length ($\epsilon$)')
axes[0,1].set_xlabel('Radius (pc)')
# axes[0,1].set_ylabel('Pressure ($10^{-13}$ dyne/cm$^2$)')
axes[0,1].set_ylabel('Pressure (dyne/cm$^2$) [log]')
axes[0,1].set_title('Pressure Profile')
# axes[0,1].set_xscale('log')
# axes[0,1].set_yscale('log')
axes[0,1].legend()
axes[0,1].grid(True, alpha=0.3)


# Top-right: Gravitational Potential
axes[0,2].plot(np.log10(r/pc),Phi(r),color='#cc6600', label='Gravitational Potential')
axes[0,2].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[0,2].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length ($\epsilon$)')
axes[0,2].set_xlabel('Radius (pc)')
axes[0,2].set_ylabel('$\Phi$ (erg/g)')
axes[0,2].set_title('Gravitational Potential')
axes[0,2].legend()
axes[0,2].grid(True, alpha=0.3)

# Bottom-left: Temperature
axes[1,0].plot(np.log10(r/pc), np.log10(Temperature(r)), 'r-', label='Temperature')
axes[1,0].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[1,0].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length ($\epsilon$)')
axes[1,0].set_xlabel('Radius (pc)')
axes[1,0].set_ylabel('Temperature (K) [log]')
axes[1,0].set_title('Temperature Profile (Isothermal)')
axes[1,0].legend()
axes[1,0].grid(True, alpha=0.3)



# Bottom-middle: Time scales
axes[1,1].plot(np.log10(r/pc),np.log10(cooling_time(r)/(s_yr)),color='purple', label='Cooling Time')
axes[1,1].plot(np.log10(r/pc),np.log10(free_fall_time(r)/(s_yr)),color='#b81d1d', label='Free Fall Time')
axes[1,1].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[1,1].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length ($\epsilon$)')
axes[1,1].set_xlabel('Radius (pc)')
axes[1,1].set_ylabel(r'Time (yr) [log]')
axes[1,1].set_title('Time Scales')
axes[1,1].legend()
axes[1,1].grid(True, alpha=0.3)

# Bottom-right: Entropy
axes[1,2].plot(np.log10(r/pc),entropy(r),color='orange', label='Entropy')
axes[1,2].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[1,2].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length ($\epsilon$)')
axes[1,2].set_xlabel('Radius (pc)')
axes[1,2].set_ylabel(r'P/$\rho^{\gamma}$')
axes[1,2].set_title('Entropy Profile')
axes[1,2].legend()
axes[1,2].grid(True, alpha=0.3)

plt.tight_layout()

save_path = os.path.join(script_dir, 'adiabatic_profiles.png')
fig.savefig(save_path, dpi=300)
print("Saved plots to:", save_path)
plt.close(fig)