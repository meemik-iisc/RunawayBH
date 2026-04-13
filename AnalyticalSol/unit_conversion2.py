import numpy as np
import matplotlib.pyplot as plt
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
# K       = 2.91e28       #g^(-2/3)cm^4s^(-2)
gamma   = 5.0/3.0       #Polytropic index
gm1     = gamma-1

#Code Units
L_code  = 1.0*pc        #cm
v_code  = 1e8            #cm/s
rho_code= 1.0*mp         #g/cm^3

M_code  = rho_code*L_code**3
t_code  = L_code/v_code
P_code  = rho_code*v_code**2
T_code  = (mu*mp/kB)*(v_code**2)
#Parameters
v_bh    = 1000*1e5                  #cm/s
M_bh    = 2e7*Msun                  #g
R_bondi = 2*G*M_bh/(v_bh**2)        #cm
epsilon = 0.1*R_bondi               #cm
rho_vir = mp                        #g/cm^3
rho_cgm = 1e-3*mp                   #g/cm^3
T_cgm   = 1e6                       #K
cs_cgm  = np.sqrt(kB*T_cgm/(mu*mp))
P_cgm   = rho_cgm*kB*T_cgm/(mu*mp)
P_ram   = rho_cgm*(v_bh**2)

#Calculate Virial Temp at 1kpc
r_vir   = 0.5*kpc                   #cm          

# #Calculate Polytropic Constant
# K = ((rho_cgm*(cs_cgm**2))/(rho_vir**(gamma)))
# # print(f"Polytropic Constant = {K:.2e}")
# #Calculate Temperature
# T_vir = K*mu*mp*(rho_vir**(gm1))/kB

# #Outflow Parameters
# v_wind  = 1e2*1e5       #cm/s
# R0      = 0.25*kpc
# t_sys   = 70*1e6*s_yr   #s
# eta     = 0.1
# M_dot_w = ((v_bh**2)*4*np.pi*(R0**2)*rho_cgm)/(v_wind)
# Rinj    = 50*pc

# Bubble Parameters
n_b         = 1.0
rho_b       = n_b*mu*mp
# K_b         = P_cgm/(rho_b**gamma)
K_b         = P_ram/(rho_b**gamma)
r_b         = 0.1*pc
#Inner bubble parameters
n_c         = 1e2
rho_c       = n_c*mu*mp

# Calculate bubble density
def phi(r):
    return (-1*G*M_bh)/np.sqrt(r**2+epsilon**2)

def rho_poly(r):
    phi_b, phi_r = phi(r_b), phi(r)
    rhs = np.power(rho_b, gm1) - (gm1 / (gamma * K_b)) * (phi_r - phi_b)
    rho = np.zeros_like(r)
    rho[rhs > 0] = np.power(rhs[rhs > 0], 1/gm1)
    return rho

# =============================================================================
# FIND r_cap WHERE ρ_env(r_cap) = ρ_cap
# =============================================================================
def r_crossing(r):
    rho_r = rho_poly(r)

    # Find crossing
    diff = rho_r - rho_c
    crossings = np.where(np.diff(np.sign(diff)))[0]

    if len(crossings) > 0:
        i_cross = crossings[0]
        # Linear interpolation
        frac = (rho_c - rho_r[i_cross]) / (rho_r[i_cross+1] - rho_r[i_cross])
        r_cap = r[i_cross] + frac * (r[i_cross+1] - r[i_cross])
        print(f"r_cap               = {r_cap:.2e} cm \t\t = {r_cap/L_code:.2e} L_code")
    else:
        print("No crossing! Check rho_env max vs rho_cap")
        r_cap = 0
    return r_cap










print("="*10,"Code Units","="*90)
print(f"Code Length         = {L_code:.2e} cm \t\t = {L_code/kpc:.2e} kpc")
print(f"Code Mass           = {M_code:.2e} g \t\t = {M_code/Msun:.2e} Msun")
print(f"Code Time           = {t_code:.2e} s \t\t = {t_code/(1e6*s_yr):.2e} Myr")
print(f"Code Density        = {rho_code:.2e} g/cm^3 \t\t = {rho_code/mp:.2e} mp/cm^3")
print(f"Code Velocity       = {v_code:.2e} cm/s \t\t = {v_code/1e5:.2e} km/s")
print(f"Code Temperature    = {T_code:.2e} K")
print(f"Code Pressure       = {P_code:.2e} dyne/cm^2")
# print(f"Code Pressure1    = {P_code1:.2e} dyne/cm^2")
#Calculate constants in code units
G_code  = G*(rho_code*t_code**2)
Kb_code  = K_b*((rho_code**gamma)/P_code)
# Kc_code  = K_c*((rho_code**gamma)/P_code)


print("="*10,"Constants in Code Units","="*77)
print(f"G                   = {G:.2e}cm^3/g/s^2 \t = {G_code:.2e} G_code")
# print(f"K                   = {K:.2e} K_cgs \t\t = {K_code:.2e} K_code")

#Calculate Parameters in code units
print("="*10,"BH Parameters in Code Units","="*76)
# print(f"Virial Radius       = {r_vir:.2e} cm \t\t = {r_vir/L_code:.2e} L_code")
# print(f"Virial Density      = {rho_vir:.2e} g/cm^3 \t\t = {rho_vir/rho_code:.2e} rho_code")
# print(f"Virial Temp         = {T_vir:.2e} K \t\t = {T_vir/T_code:.2e} T_code")
print(f"Black Hole Mass     = {M_bh:.2e} g \t\t = {M_bh/M_code:.2e} M_code")
print(f"Black Hole Vel      = {v_bh:.2e} cm/s \t\t = {v_bh/v_code:.2e} v_code")
print(f"Bondi Radius        = {R_bondi:.2e} cm \t\t = {R_bondi/L_code:.2e} L_code")
print(f"Epsilon             = {epsilon:.2e} cm \t\t = {epsilon/L_code:.2e} L_code")
print(f"CGM Density         = {rho_cgm:.2e} g/cm^3 \t\t = {rho_cgm/rho_code:.2e} rho_code")
print(f"CGM Temp            = {T_cgm:.2e} K \t\t = {T_cgm/T_code:.2e} T_code")
print(f"CGM Sound Speed     = {cs_cgm:.2e} cm/s \t\t = {cs_cgm/v_code:.2e} v_code")
print(f"CGM Pressure        = {P_cgm:.2e} dyne/cm^2 \t = {P_cgm/P_code:.2e} P_code")

print("="*10,"Bubble parameters in Code Units","="*77)
print(f"r_b                 = {r_b:.2e} cm \t\t = {r_b/L_code:.2e} L_code")
print(f"rho_b               = {rho_b:.2e} g/cm^3 \t\t = {rho_b/rho_code:.2e} rho_code")
print(f"K_b                 = {K_b:.2e} dyne/cm^2 \t = {Kb_code:.2e} K_code")

# print(f"r_cap               = {r_cap:.2e} cm \t\t = {r_cap/L_code:.2e} L_code")
print(f"rho_c               = {rho_c:.2e} g/cm^3 \t\t = {rho_c/rho_code:.2e} rho_code")


# print("="*10,"Outflow Parameters in Code Units","="*76)
# print(f"Outflow Velocity    = {v_wind:.2e} cm/s \t\t = {v_wind/v_code:.2e} v_code")
# print(f"Outflow M_dot       = {M_dot_w*s_yr/Msun:.2e} Msun/yr \t\t = {(M_dot_w*t_code)/M_code:.2e} M_code/t_code")
# print(f"Injection radius    = {Rinj:.2e} cm \t\t = {Rinj/L_code:.2e} L_code")


# =============================================================================
# PIECEWISE DENSITY
# =============================================================================
r_plot = np.logspace(np.log10(1e-3*pc), np.log10(pc), 2000)
r_cap = r_crossing(r_plot)

def dens(r):
    """rho(r): polytrope until r_cap, then flat core, then CGM"""
    rho_out = np.zeros_like(r)
    mask = r<r_cap
    mask2 = (r >= r_cap) & (r < r_b)
    # Core: ρ = rho_c for r<r_cap
    rho_out[mask] = rho_c 
    #envelope: polytrope for r_cap ≤ r < r_b
    rho_out[mask2] = rho_poly(r[mask2])   
    # CGM: r ≥ r_b
    rho_out[r >= r_b] = rho_cgm
    return rho_out

# rho_total = dens(r_plot)

def bubble_mass(r):
    rho_vals = dens(r)
    M_b = 4*np.pi*np.trapezoid(r**2*rho_vals, r)
    return M_b

def pressure(r):
    pres_out = np.zeros_like(r)
    mask = r < r_cap
    mask2 = (r >= r_cap) & (r < r_b)
    #Core: P = P_c for r<r_cap
    pres_out[mask] = K_b*(rho_c**gamma)
    # Envelope: polytrope for r_cap ≤ r < r_b  
    pres_out[mask2] = K_b*np.power(rho_poly(r[mask2]),gamma)
    # CGM: r ≥ r_b
    pres_out[r >= r_b] = P_cgm
    
    return pres_out
    # if r < r_cap:
    #     return K_b * np.power(rho_c, gamma)
    # elif r_cap <= r < r_b:
    #     return K_b * np.power(rho_poly(r), gamma)
    # else:
    #     return P_cgm

def temperature(r):
    return (pressure(r)*mu*mp)/(dens(r)*kB)



M_b = bubble_mass(r_plot)
print("="*10,"Inner bubble parameters","="*77)

# print(f"(r_cap) = {r_cap:.2e} cm \t\t = {r_cap/L_code:.2e} L_code")
# print(f"rho_c = {rho_c:.2e} g/cm^3 \t\t = {rho_c/rho_code:.2e} rho_code")
print(f"M_b  = {M_b/Msun:.2e} Msun \t\t = {M_b/M_code:.2e} M_code")



print("="*112)
# =============================================================================
# PLOT
# =============================================================================
fig, axes = plt.subplots(3, 1, figsize=(12, 10))

# Density
axes[0].loglog(r_plot/pc, dens(r_plot)/mp, 'b-', lw=2)
axes[0].axvline(r_b/pc, ls='--', c='k', lw=2)
axes[0].axhline(rho_b/mp, ls=':', c='r', lw=2)
axes[0].set_title('Density') 
axes[0].set_ylabel('rho [mp/cm³]')
axes[0].grid(True, alpha=0.3)

# Pressure  
axes[1].loglog(r_plot/pc, pressure(r_plot), 'g-', lw=2)
axes[1].axvline(r_b/pc, ls='--', c='k', lw=2)
axes[1].axhline(P_cgm, ls=':', c='r', lw=2)
axes[1].set_title('Pressure')
axes[1].set_ylabel('P [dyne/cm²]')
axes[1].grid(True, alpha=0.3)
# Mass
# axes[1,1].semilogx(r_int[:-1]/pc, M_enc/Msun, 'm-', lw=2)
# axes[1,1].axvline(r_b/pc, ls='--', c='k', lw=2)
# axes[1,1].axhline(M_b/Msun, ls=':', c='r', lw=2)
# axes[1,1].set_xlabel('r [pc]'); axes[1,1].set_ylabel('M_enc [Msun]')
# axes[1,1].set_title('Enclosed Mass'); axes[1,1].grid(True, alpha=0.3)
# Temperature
axes[2].loglog(r_plot/pc, temperature(r_plot), 'r-', lw=2)
axes[2].axvline(r_b/pc, ls='--', c='k', lw=2)
axes[2].axhline(T_cgm, ls=':', c='orange', label='T_cgm')
axes[2].set_xlabel('r [pc]')
axes[2].set_ylabel('T [K]')
axes[2].set_title('Temperature')
axes[2].grid(True, alpha=0.3)
axes[2].legend()

# # Mass
# axes[1,1].semilogx(r_int[:-1]/pc, M_enc/Msun, 'm-', lw=2)
# axes[1,1].axvline(r_b/pc, ls='--', c='k', lw=2)
# axes[1,1].axhline(M_b/Msun, ls=':', c='r', lw=2)
# axes[1,1].set_xlabel('r [pc]'); axes[1,1].set_ylabel('M_enc [Msun]')
# axes[1,1].set_title('Enclosed Mass'); axes[1,1].grid(True, alpha=0.3)

# plt.suptitle(f'Polytropic Bubble: ρ(r_b)={rho_b/mp:.0f} mp/cm³, M_b={M_b_Msun:.0f} Msun', fontsize=14)
plt.tight_layout()
script_dir = os.path.dirname(os.path.abspath(__file__))
# filename = os.path.join(script_dir, 'bubble_full_cgm.png')
filename = os.path.join(script_dir, 'bubble_full_ram.png')
plt.savefig(filename, dpi=150, bbox_inches='tight')
plt.show()

