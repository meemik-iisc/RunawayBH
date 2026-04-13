#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt

# Constants (CGS)
G, kB, mp, pc, Msun = 6.67e-8, 1.38e-16, 1.67e-24, 3.086e18, 1.989e33
mu, gamma, gm1 = 1.0, 5/3, 2/3

# BH/CGM (CGS)
M_bh = 2e7 * Msun
v_bh = 1000 * 1e5
R_bondi = 2 * G * M_bh / v_bh**2
epsilon = 0.1 * R_bondi

n_cgm, T_cgm = 1e-3, 1e6
rho_cgm = n_cgm * mu * mp
cs_cgm = np.sqrt(kB * T_cgm / (mu * mp))
p_cgm = rho_cgm * cs_cgm**2

# Bubble
r_b = 0.1 * pc
n_b = 0.1
rho_b = n_b * mu * mp
K_b = p_cgm / np.power(rho_b, gamma)

M_target = 1e6 * Msun

def Phi(r): return -G * M_bh / np.sqrt(r**2 + epsilon**2)

def rho_poly(r):
    phi_b, phi_r = Phi(r_b), Phi(r)
    rhs = np.power(rho_b, gm1) - (gm1 / (gamma * K_b)) * (phi_r - phi_b)
    rho = np.zeros_like(r)
    rho[rhs > 0] = np.power(rhs[rhs > 0], 1/gm1)
    return rho

# =============================================================================
print("=== SETUP ===")
print(f"r_b = {r_b/pc:.2f} pc")
print(f"rho_b(r_b) = {rho_b/mp:.1e} mp/cm³")
print(f"K_b = {K_b:.2e}")
print()

# Radial grid
r_min = 1e-5 * pc
N = 2000
r = np.logspace(np.log10(epsilon), np.log10(r_b), N)

rho = rho_poly(r)
p = K_b * np.power(rho, gamma)
T = p / rho * (mu * mp / kB)  # T = p / (ρ k_B / μ m_p)
frac_positive = np.sum(rho > 0) / len(r)
print(f"ρ>0 fraction: {frac_positive:.1%} ({np.sum(rho>0)}/{N} cells)")

# Better mass integration (trapezoid + check where rho>0)
mask = rho > 1e-10 * rho_b
if np.sum(mask) > 0:
    r_int = r[mask]
    rho_int = rho[mask]
    dr = np.diff(r_int)
    dM = 4 * np.pi * r_int[:-1]**2 * rho_int[:-1] * dr
    M_enc = np.cumsum(dM)
    M_b = M_enc[-1]
else:
    print("ERROR: rho everywhere <= 0!")
    M_b = 0

M_b_Msun = M_b / Msun
factor = M_target / M_b

print("\n=== RESULTS ===")
print(f"M_bubble = {M_b:.2e} g = {M_b_Msun:.1e} Msun")
print(f"vs target = {factor:.1f}x")

print("\n=== TUNING ===")
print(f"For 1e6 Msun → ρ_b(r_b) = {rho_b/mp * factor:.1e} mp/cm³")
print(f"Edit n_b = {n_b * factor:.1e}")

# =============================================================================
# PLOTS: ρ, p, T, M
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Density
axes[0,0].loglog(r/pc, rho/mp, 'b-', lw=2)
axes[0,0].axvline(r_b/pc, ls='--', c='k', lw=2)
axes[0,0].axhline(rho_b/mp, ls=':', c='r', lw=2)
axes[0,0].set_title('Density'); axes[0,0].set_ylabel('rho [mp/cm³]')
axes[0,0].grid(True, alpha=0.3)

# Pressure  
axes[0,1].loglog(r/pc, p, 'g-', lw=2)
axes[0,1].axvline(r_b/pc, ls='--', c='k', lw=2)
axes[0,1].axhline(p_cgm, ls=':', c='r', lw=2)
axes[0,1].set_title('Pressure'); axes[0,1].set_ylabel('P [dyne/cm²]')
axes[0,1].grid(True, alpha=0.3)

# Temperature
axes[1,0].loglog(r/pc, T, 'r-', lw=2)
axes[1,0].axvline(r_b/pc, ls='--', c='k', lw=2)
axes[1,0].axhline(T_cgm, ls=':', c='orange', label='T_cgm')
axes[1,0].set_xlabel('r [pc]'); axes[1,0].set_ylabel('T [K]')
axes[1,0].set_title('Temperature'); axes[1,0].grid(True, alpha=0.3); axes[1,0].legend()

# Mass
axes[1,1].semilogx(r_int[:-1]/pc, M_enc/Msun, 'm-', lw=2)
axes[1,1].axvline(r_b/pc, ls='--', c='k', lw=2)
axes[1,1].axhline(M_b/Msun, ls=':', c='r', lw=2)
axes[1,1].set_xlabel('r [pc]'); axes[1,1].set_ylabel('M_enc [Msun]')
axes[1,1].set_title('Enclosed Mass'); axes[1,1].grid(True, alpha=0.3)

plt.suptitle(f'Polytropic Bubble: ρ(r_b)={rho_b/mp:.0f} mp/cm³, M_b={M_b_Msun:.0f} Msun', fontsize=14)
plt.tight_layout()
plt.savefig('bubble_full.png', dpi=150, bbox_inches='tight')
plt.show()