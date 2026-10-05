import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os


# Constants
G       = 6.67e-8       # cm^3/g/s^2
kB      = 1.38e-16      # erg/K
mp      = 1.67e-24      # g
kpc     = 3.086e21      # cm
pc      = 3.086e18      # cm
Msun    = 1.989e33      # g
mu      = 1.0
s_yr    = 3.154e7       # s
s_Myr   = 1e6*s_yr
gamma   = 5.0/3.0       # Polytropic index
gm1     = gamma-1
kms     = 1e5


# ============================================================================
# COOLING TABLE
# ============================================================================
cool_dir = "/data/data_meemik/RunawayBH/AnalyticalSol/cooltable.dat"
data = np.loadtxt(cool_dir)
T_tab = data[:, 0]
Lambda_tab = data[:, 1]
Lambda_interp = interp1d(T_tab, Lambda_tab, kind='linear', 
                         bounds_error=False, fill_value='extrapolate')


def cool_lambda(T):
    """Safe cooling function with clamping"""
    T_safe = np.clip(T, T_tab.min() * 0.95, T_tab.max() * 1.05)
    return Lambda_interp(T_safe)


def t_cool(n, T):
    return (gamma * kB * T) / (gm1 * n * cool_lambda(T)) / s_Myr


# Create arrays
n_arr = np.logspace(-4, 2, 100)
T_arr = np.logspace(4, 8, 100)


# Compute cooling time array
t_cool_arr = np.zeros((len(T_arr), len(n_arr)))
for i in range(len(n_arr)):
    for j in range(len(T_arr)):
        t_cool_arr[j, i] = t_cool(n_arr[i], T_arr[j])


# For LogNorm, get valid (non-nan) range
t_cool_valid = t_cool_arr[~np.isnan(t_cool_arr)]
vmin = np.min(t_cool_valid)
vmax = np.max(t_cool_valid)


# Create the plot with improved aesthetics
fig, ax = plt.subplots(figsize=(10, 8))


# Set bad values color on the colormap
cmap = "coolwarm"


# Create pcolormesh with LogNorm for logarithmic colorbar
t_cool_masked = np.ma.masked_invalid(t_cool_arr)
pcm = ax.pcolormesh(n_arr, T_arr, t_cool_masked, 
                    cmap=cmap, shading='auto',
                    norm=plt.matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax))


# Set scales
ax.set_xscale('log')
ax.set_yscale('log')


# Labels with better formatting
ax.set_xlabel(r'$n$ (cm$^{-3}$)', fontsize=14, fontweight='bold')
ax.set_ylabel(r'$T$ (K)', fontsize=14, fontweight='bold')
ax.set_title(r'Cooling Time $t_{\rm cool}$ (Myr)', fontsize=16, fontweight='bold', pad=15)


# Add grid
ax.grid(True, which='major', linestyle='--', alpha=0.8, linewidth=0.5, color='k')
ax.grid(True, which='minor', linestyle='--', alpha=0.3, linewidth=0.5)


# ============================================================================
# PLOT MULTIPLE CONTOURS WITH CUSTOM LABELS
# ============================================================================
contour_levels = [0.1, 1.0, 10.0, 100.0]

# Plot contours
cs = ax.contour(n_arr, T_arr, t_cool_arr, levels=contour_levels, 
                colors='black', linewidths=2.0, linestyles="dashed", alpha=0.9)

# Add custom labels to contours
ax.clabel(cs, levels=contour_levels, inline=True, 
          fmt={lvl: f'{lvl:g} Myr' for lvl in contour_levels}, 
          fontsize=11)


# Create colorbar with log scale
cbar = plt.colorbar(pcm, ax=ax, pad=0.02)
cbar.set_label(r'$t_{\rm cool}$ (Myr)', fontsize=13, fontweight='bold')
cbar.ax.tick_params(labelsize=11)



# Add minor ticks
ax.minorticks_on()


# Tight layout
plt.tight_layout()


# Get script directory and save
script_dir = os.path.dirname(os.path.abspath(__file__))
output_path = os.path.join(script_dir, 'cooling_map_with_contour.png')


plt.savefig(output_path, dpi=300, bbox_inches='tight')
print(f"Figure saved to: {output_path}")
print(f"Cooling time range: {vmin:.3e} to {vmax:.3e} Myr")
print(f"t_cool = 1 Myr contour plotted")