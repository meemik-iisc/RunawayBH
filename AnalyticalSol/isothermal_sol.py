import numpy as np
import matplotlib.pyplot as plt

#Constants
G   = 6.67e-8       #cm^3/g/s^2
kB  = 1.38e-16      #erg/K
mp  = 1.67e-24      #g
pc  = 3.086e18      #cm
Msun= 1.989e33      #g
mu  = 1.0

#Parameters
T_amb   = 1e7                       #K
cs_amb  = np.sqrt(kB*T_amb/(mu*mp)) #cm/s
v_bh    = 1000*1e5                  #cm/s
print(f"Ambient sound speed: {cs_amb/1e5:.2e} km/s")
rho0    = 1e-3*mp                   #g/cm^3
M_bh    = 2e7*Msun                  #g
R_bondi = 2*G*M_bh/(cs_amb**2+v_bh**2)      #cm
print(f"Bondi radius: {R_bondi/pc:.2e} pc")
epsilon = R_bondi                    #cm
# epsilon = R_bondi

#Virial Tenp at r=1kpc
T_vir = G*M_bh*mp/(1e3*pc*kB)
print(f"Virial Temperature at 1kpc={T_vir:.2e}")
def Phi(r):
    return -1*G*M_bh/np.sqrt(r**2+epsilon**2)
def density_profile(r):
    return rho0*np.exp((-Phi(r)+Phi(0))/(cs_amb**2))

def pressure_profile(r):
    return density_profile(r)*cs_amb**2

def temperature_profile(r):
    return (pressure_profile(r)*mu*mp)/(density_profile(r)*kB)

#Plot profiles
r = np.linspace(0.01*pc,1e3*pc,10000)

# Create 2x2 subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('Isothermal Hydrostatic Equilibrium around Black Hole')

# Top-left: Density
axes[0,0].plot(np.log10(r/pc), np.log10(density_profile(r)/mp), 'b-', label='Density')
axes[0,0].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[0,0].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length')
axes[0,0].set_xlabel('Radius (pc)')
axes[0,0].set_ylabel('Density ($m_p/cm^3$) [log]')
axes[0,0].set_title('Density Profile')
# axes[0,0].set_xscale('log')
# axes[0,0].set_yscale('log')
axes[0,0].legend()
axes[0,0].grid(True, alpha=0.3)

# Top-right: Pressure  
axes[0,1].plot(np.log10(r/pc), np.log10(pressure_profile(r)), 'g-', label='Pressure')
axes[0,1].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[0,1].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length')
axes[0,1].set_xlabel('Radius (pc)')
# axes[0,1].set_ylabel('Pressure ($10^{-13}$ dyne/cm$^2$)')
axes[0,1].set_ylabel('Pressure (dyne/cm$^2$) [LOG]')
axes[0,1].set_title('Pressure Profile')
# axes[0,1].set_xscale('log')
# axes[0,1].set_yscale('log')
axes[0,1].legend()
axes[0,1].grid(True, alpha=0.3)

# Bottom-left: Temperature
axes[1,0].plot(np.log10(r/pc), temperature_profile(r), 'r-', label='Temperature')
axes[1,0].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[1,0].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length')
axes[1,0].set_xlabel('Radius (pc)')
axes[1,0].set_ylabel('Temperature (K)')
axes[1,0].set_title('Temperature Profile (Isothermal)')
axes[1,0].legend()
axes[1,0].grid(True, alpha=0.3)

# Bottom-right: Gravitational Potential
axes[1,1].plot(np.log10(r/pc),Phi(r)/1e12,color='orange', label='Gravitational Potential')
axes[1,1].axvline(np.log10(R_bondi/pc), color='k', linestyle='--', label='Bondi Radius')
axes[1,1].axvline(np.log10(epsilon/pc), color='gray', linestyle='--', label='Softening Length')
axes[1,1].set_xlabel('Radius (pc)')
axes[1,1].set_ylabel('$\Phi$ ($10^{12}$erg/g)')
axes[1,1].set_title('Gravitational Potential')
axes[1,1].legend()
axes[1,1].grid(True, alpha=0.3)

plt.tight_layout()
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
save_path = os.path.join(script_dir, 'isothermal_profiles.png')
fig.savefig(save_path, dpi=300)
print("Saved plots to:", save_path)
plt.close(fig)