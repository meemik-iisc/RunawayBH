#!/usr/bin/env python3
"""
Simple bubble mass calculator for your AthenaK setup.
Just edit the PARAMETERS section and run.
"""

import numpy as np
import matplotlib.pyplot as plt

#Constants
G           = 6.67e-8       #cm^3/g/s^2
kB          = 1.38e-16      #erg/K
mp          = 1.67e-24      #g
kpc         = 3.086e21      #cm
pc          = 3.086e18      #cm
Msun        = 1.989e33      #g
mu          = 1.0
s_yr        = 3.154e7       #s
gamma       = 5.0/3.0       #Polytropic index
gm1         = gamma-1


# =============================================================================
# PARAMETERS (CGS UNITS)
# =============================================================================

#BH parameters
M_bh        = 2e7*Msun   
v_bh        = 1000*1e5             
R_bondi     = 2*G*M_bh/(v_bh**2)        
epsilon     = 0.1*R_bondi               
# CGM parameters
n_cgm       = 1.0e-3
rho_cgm     = n_cgm*mu*mp
T_cgm       = 1e6
cs_cgm      = np.sqrt(kB*T_cgm/(mu*mp))
p_cgm       = rho_cgm * cs_cgm**2
#bubble parameters
r_b         = 0.1*pc   
n_b         = 1.0e3
rho_b       = n_b*mu*mp
K_b         = p_cgm/np.power(rho_b, gamma)   

# Target mass check
M_target    = 1e6*Msun  # desired bubble mass

# =============================================================================

def Phi(r):
    """Bondi potential"""
    return -1*G*M_bh/np.sqrt(r**2 + epsilon**2)

def rho_poly(r):
    """Polytropic hydrostatic density"""
    phi_b = Phi(r_b)
    phi_r = Phi(r)
    
    rhs = np.power(rho_b, gm1) - (gm1/(gamma*K_b))*(phi_b-phi_r)
    rho = np.zeros_like(r)
    rho[rhs > 0] = np.power(rhs[rhs > 0], 1/gm1)
    return rho

# =============================================================================
# COMPUTE
# =============================================================================
print(f"Target      = {M_target/Msun:,.1e} Msun")
print("================ BUBBLE PROFILE =================")
print(f"r_b         = {r_b/pc:.2f} pc")
print(f"rho_b       = {rho_b/mp:.1e}mp/cc")
print(f"p_cgm       = {p_cgm:>8.2e} dyne/cm^2")
print(f"K_bubble    = {K_b:>8.2e}")

# Radial integration
N = 1000   
r = np.logspace(np.log10(1.0e-4*pc), np.log10(r_b), N)
# print(f"r_min = {r[1]/pc:.1e} pc")
rho = rho_poly(r)
# print(f"rho(r_b)    = {rho[-1]/mp:.1e}mp/cc")

# Mass
dr = np.diff(r)
dM = 4 * np.pi * r[:-1]**2 * rho[:-1] * dr
M_enc = np.cumsum(dM)
M_b = M_enc[-1]

print("================ BUBBLE MASS =================")
print(f"M_bubble = {M_b:>9.2e}g \t = {M_b/Msun:>9.2e} Msun")
# print(f"vs target = {M_target_Msun:>9.2e} Msun (factor {M_b_Msun/M_target_Msun:.2f})")

# # =============================================================================
# # PLOT
# # =============================================================================
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

# ax1.loglog(r, rho, 'b-', lw=2)
# ax1.axvline(r_b, color='k', ls='--', lw=2, label=f'r_b={r_b}')
# ax1.axhline(rho_edge, color='r', ls=':', lw=2, label=f'ρ(r_b)={rho_edge:.1e}')
# ax1.set_xlabel('r [code]'); ax1.set_ylabel('ρ [code]')
# ax1.legend(); ax1.grid(True, alpha=0.3)

# ax2.semilogx(r, M_enc, 'g-', lw=2)
# ax2.axvline(r_b, color='k', ls='--', lw=2)
# ax2.axhline(M_b, color='r', ls=':', lw=2, label=f'M_b={M_b:.1e}')
# ax2.set_xlabel('r [code]'); ax2.set_ylabel('Enclosed mass [code]')
# ax2.legend(); ax2.grid(True, alpha=0.3)

# plt.suptitle(f'Bubble: ρ_edge={rho_edge:.1e}, M_b={M_b_Msun:.0f} Msun')
# plt.tight_layout()
# plt.savefig('bubble_mass.png', dpi=150, bbox_inches='tight')
# plt.show()

# print("\n=== TUNING GUIDE ===")
# factor = M_target_Msun / M_b_Msun
# print(f"To get {M_target_Msun:,.0f} Msun:")
# print(f"  Multiply rho_edge by {factor:.1f}")
# print(f"  New rho_edge = {rho_edge*factor:.1e}")
