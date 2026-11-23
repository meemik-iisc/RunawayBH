import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import os
output_dir = os.path.dirname(os.path.abspath(__file__))

#Constants
kB      = 1.38e-16      #erg/K
mp      = 1.67e-24      #g
kpc     = 3.086e21      #cm
Msun    = 1.989e33      #g
mu      = 1.0
s_yr    = 3.154e7       #s
kms     = 1e5

#Flags
fill_regions = True

#Parameters
v_bh    = 1000*1e5                  #cm/s
R0      = 1*kpc                     #cm, minimum injection radius
Mdot_max= 0.3*Msun/s_yr             #g/s, mass budget
n_cgm   = 1                         #cm^-3
rho_cgm = n_cgm*mp*mu
T_cgm   = 1e6                       #K

#Outflow Parameters
v_wind  = 1e2*kms       #cm/s

#Code Units
v_code  = 1e3*kms            #cm/s
rho_code= 1.0*mp             #g/cm^3
T_code  = (mu*mp/kB)*(v_code**2)

# ============================================================================
# COOLING TABLE
# ============================================================================
cool_dir = os.path.join(output_dir,"cooltable.dat")
data = np.loadtxt(cool_dir)
T_tab = data[:, 0]
Lambda_tab = data[:, 1]
Lambda_interp = interp1d(T_tab, Lambda_tab, kind='linear',
                         bounds_error=False, fill_value='extrapolate')

def cool_lambda(T):
    """Safe cooling function with clamping"""
    T_safe = np.clip(T, T_tab.min() * 0.95, T_tab.max() * 1.05)
    return Lambda_interp(T_safe)

# ============================================================================
# COOLING TIME CALCULATIONS
# ============================================================================

def t_cool_scgm(ncgm):
    T_scgm          = (3*mu*mp*(v_bh**2))/(16*kB)
    return(9*mu*mp*(v_bh**2))/(128*ncgm*cool_lambda(T_scgm))

def t_cool_sw(ncgm, vw):
    T_sw            = (3*mu*mp*(vw**2))/(16*kB)
    tcool_sw        = (9*mu*mp*(vw**4))/(128*ncgm*(v_bh**2)*cool_lambda(T_sw))
    return tcool_sw

def v_w_R0(ncgm, R0=R0, Mdot=Mdot_max):
    """Wind speed (cm/s) at which R0 is reached for mass-loss rate Mdot"""
    return 4*np.pi*mu*mp*ncgm*(v_bh**2)*(R0**2)/Mdot

print("="*30,"Parameters","="*40)
print(f"CGM Density             = {rho_cgm/mp:.2e} mp/cm^3 \t\t = {rho_cgm/rho_code:.2e} rho_code")
print(f"CGM Temp                = {T_cgm:.2e} K \t\t = {T_cgm/T_code:.2e} T_code")
print(f"Outflow Velocity        = {v_wind:.2e} cm/s \t\t = {v_wind/v_code:.2e} v_code")

n_cgm   = np.logspace(-4, 3, 100)
vw      = np.array([1e2, 1e3, 1e4])      #km/s
vw_cgs  = vw*kms

# ============================================================================
# PLOT 1
# ============================================================================

fig1,ax1 = plt.subplots(figsize=(8,6))
ax1.plot(np.log10(n_cgm), np.log10(t_cool_scgm(n_cgm)/s_yr),  linewidth=2, markersize=2, color='red', label="Shocked CGM")
ax1.plot(np.log10(n_cgm), np.log10(t_cool_sw(n_cgm, vw_cgs[0])/s_yr),  linewidth=2, linestyle='solid', markersize=2, color='blue', label="Shocked Wind (1e2 km/s)")
ax1.plot(np.log10(n_cgm), np.log10(t_cool_sw(n_cgm, vw_cgs[1])/s_yr),  linewidth=2, linestyle='dashed', markersize=2, color='blue', label="Shocked Wind (1e3 km/s)")
ax1.plot(np.log10(n_cgm), np.log10(t_cool_sw(n_cgm, vw_cgs[2])/s_yr),  linewidth=2, linestyle='dotted', markersize=2, color='blue', label="Shocked Wind (1e4 km/s)")
ax1.set_xlabel(r'$\rho_{CGM}$ (mp/cc) [log]', fontsize = 14)
ax1.set_ylabel(r'$t_{cool}$ (yr) [log]', fontsize = 14)
ax1.tick_params(axis='both', which='major', labelsize=14)
# Customize log grid
ax1.minorticks_on()
ax1.grid(True, which='both', ls='-', alpha=0.3)
ax1.grid(which='major', alpha=0.7, linewidth=1.2)
ax1.grid(which='minor', alpha=0.5, linewidth=1.0)
ax1.legend(loc='upper right', fontsize=12)
save_filename1 = os.path.join(output_dir,f"t_cool_vs_ncgm.png")
plt.savefig(save_filename1)
plt.close(fig1)


v_w     = np.logspace(1, 4, 100)
v_w_cgs = v_w*kms
ncgm    = [1e-3,1, 1e3]


# ============================================================================
# PLOT 2
# ============================================================================
fig2,ax2 = plt.subplots(figsize=(8,6))
ax2.plot(np.log10(v_w), np.log10(np.ones(len(v_w_cgs))*t_cool_scgm(ncgm[0])/s_yr),  linewidth=2, linestyle='solid', markersize=2, color='red')
ax2.plot(np.log10(v_w), np.log10(t_cool_sw(ncgm[0], v_w_cgs)/s_yr),  linewidth=2, linestyle='dashed', markersize=2, color='red')
ax2.plot(np.log10(v_w), np.log10(np.ones(len(v_w_cgs))*t_cool_scgm(ncgm[1])/s_yr),  linewidth=2, linestyle='solid', markersize=2, color='blue')
ax2.plot(np.log10(v_w), np.log10(t_cool_sw(ncgm[1], v_w_cgs)/s_yr),  linewidth=2, linestyle='dashed', markersize=2, color='blue')
ax2.plot(np.log10(v_w), np.log10(np.ones(len(v_w_cgs))*t_cool_scgm(ncgm[2])/s_yr),  linewidth=2, linestyle='solid', markersize=2, color='green')
ax2.plot(np.log10(v_w), np.log10(t_cool_sw(ncgm[2], v_w_cgs)/s_yr),  linewidth=2, linestyle='dashed', markersize=2, color='green')
ax2.set_xlabel(r'$v_{wind}$ (km/s) [log]', fontsize = 14)
ax2.set_ylabel(r'$t_{cool}$ (yr) [log]', fontsize = 14)
ax2.tick_params(axis='both', which='major', labelsize=14)
legend1 = ax2.legend(
    handles=[
        Line2D([0], [0], color='red',    lw=2, label=r'$n_{CGM} = 10^{-3} cm^{-3}$'),
        Line2D([0], [0], color='blue',  lw=2, label=r'$n_{CGM} = 1 cm^{-3}$'), 
        Line2D([0], [0], color='green',   lw=2, label=r'$n_{CGM} = 10^3 cm^{-3}$'),
    ],
    loc='upper left', fontsize=12, frameon=True, fancybox=True, shadow=True
)

# Legend 2: Linestyles (Cases) - upper-right  
legend2 = ax2.legend(
    handles=[
        Line2D([0], [0], color='k', lw=2, linestyle='solid',    label="Shocked CGM"),
        Line2D([0], [0], color='k', lw=2, linestyle='dotted', label="Shocked Wind"),
    ],
    loc='upper right', fontsize=12, frameon=True, fancybox=True, shadow=True
)
ax2.add_artist(legend1)
# Customize log grid
ax2.minorticks_on()
ax2.grid(True, which='both', ls='-', alpha=0.3)
ax2.grid(which='major', alpha=0.7, linewidth=1.2)
ax2.grid(which='minor', alpha=0.5, linewidth=1.0)
# ax2.legend(loc='upper right', fontsize=12)
save_filename2 = os.path.join(output_dir,f"t_cool_vs_vw.png")
plt.savefig(save_filename2)
plt.close(fig2)






# ============================================================
# Plot 3: cooling-time contours in (n_cgm, v_w) plane + runs
# ============================================================
plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman', 'Times', 'STIXGeneral', 'DejaVu Serif'],
    'mathtext.fontset': 'stix',
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'xtick.top': True, 'ytick.right': True,
    'xtick.major.size': 6, 'ytick.major.size': 6,
    'xtick.minor.size': 3, 'ytick.minor.size': 3,
    'axes.linewidth': 1.0,
})

v_w_kms_2d = np.logspace(1, 4, 400)
n_cgm_2d   = np.logspace(-4, 2, 400)
NCGM, VW   = np.meshgrid(n_cgm_2d, v_w_kms_2d * kms)

TCOOL_SCGM = t_cool_scgm(NCGM) / s_yr
TCOOL_SW   = t_cool_sw(NCGM, VW) / s_yr

time_Myr  = np.array([1, 10, 60])
levels_yr = time_Myr * 1e6
fmt_time  = {t*1e6: rf'$\mathrm{{{t}\ Myr}}$' for t in time_Myr}

# Sequential blues: shortest cooling time = darkest band
band_colors = ['#9ecae1', '#c6dbef', '#eaf3fb']   # t<1, 1-10, 10-60 Myr
both_colors = ['#7a68b3', '#a59cd0', '#d3cfea']   # both phases cool: t<1, 1-10, 10-60 Myr
c_sw   = '#0b2e59'
c_scgm = '#6b6b6b'
c_ucgm = '#e8590c'
c_tcgm = '#111111'
c_R0   = '#8b1a1a'

fig3, ax3 = plt.subplots(figsize=(7.5, 5.5))
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlim(1e-4, 1e2)
ax3.set_ylim(1e1, 1e4)

# Region where the shocked wind cools within t (below each SW line)
if fill_regions:
    ax3.contourf(n_cgm_2d, v_w_kms_2d, TCOOL_SW,
                 levels=[0, *levels_yr], colors=band_colors, zorder=0)
    # Both shocked wind and shocked CGM cool within t: t_cool,max <= t
    TCOOL_BOTH = np.maximum(TCOOL_SW, TCOOL_SCGM)
    ax3.contourf(n_cgm_2d, v_w_kms_2d, TCOOL_BOTH,
                 levels=[0, *levels_yr], colors=both_colors, zorder=1)

cs_sw = ax3.contour(n_cgm_2d, v_w_kms_2d, TCOOL_SW, levels=levels_yr,
                    colors=c_sw, linewidths=1.6, linestyles='solid', zorder=2)
cs_cg = ax3.contour(n_cgm_2d, v_w_kms_2d, TCOOL_SCGM, levels=levels_yr,
                    colors=c_scgm, linewidths=1.4, linestyles=[(0, (5, 3))], zorder=2)

# Time labels placed directly on the curves (horizontal, white backing)
lab_kw = dict(fontsize=12, ha='center', va='center', zorder=7,
              bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.9))
n_lab_sw = 2e1
for t in time_Myr:
    v_line = np.interp(np.log10(t*1e6), np.log10(t_cool_sw(n_lab_sw, v_w_kms_2d*kms)/s_yr),
                       v_w_kms_2d)
    ax3.text(n_lab_sw, v_line, fmt_time[t*1e6], color=c_sw, **lab_kw)
for t in time_Myr:
    n_line = np.interp(np.log10(t*1e6), np.log10(t_cool_scgm(n_cgm_2d)/s_yr)[::-1],
                       n_cgm_2d[::-1])
    ax3.text(n_line, 7.0e3, fmt_time[t*1e6], color=c_scgm, **lab_kw)

# Mass budget: R0 >= 1 kpc above the dotted line
v_R0_kms = v_w_R0(n_cgm_2d)/kms
ax3.fill_between(n_cgm_2d, v_R0_kms, 1e4, where=v_R0_kms < 1e4, facecolor='none',
                 hatch='///', edgecolor=c_R0, linewidth=0, alpha=0.6, zorder=3)
ax3.plot(n_cgm_2d, v_R0_kms, color=c_R0, lw=2, ls=':', zorder=4)
# Label along the line, rotated to match its on-screen slope
n_lab = 3.0e-3
p1 = ax3.transData.transform((n_lab, v_w_R0(n_lab)/kms))
p2 = ax3.transData.transform((2*n_lab, v_w_R0(2*n_lab)/kms))
ang = np.degrees(np.arctan2(p2[1]-p1[1], p2[0]-p1[0]))
ax3.text(n_lab, v_w_R0(n_lab)/kms, r'$\mathrm{R_0 \geq 1\ kpc}$', color=c_R0,
         fontsize=13, rotation=ang, rotation_mode='anchor', ha='center', va='center', zorder=8,
         bbox=dict(boxstyle='round,pad=0.15', fc='white', ec='none', alpha=0.9))

# Direction in which each phase cools (t_cool < t)
dir_kw = dict(fontsize=13, ha='center', va='center', zorder=7)
arr_kw = dict(arrowstyle='-|>', lw=1.6, mutation_scale=16)
# Shocked wind: cools below the solid lines
ax3.annotate('', xy=(2.5e-3, 1.5e2), xytext=(2.5e-3, 6.5e2),
             arrowprops=dict(color=c_sw, **arr_kw), zorder=7)
ax3.text(2.5e-3, 9.0e2, 'wind cools', color=c_sw, **dir_kw)
# Shocked CGM: cools to the right of the dashed lines
ax3.annotate('', xy=(1.4e-1, 3.2e1), xytext=(3.0e-2, 3.2e1),
             arrowprops=dict(color=c_scgm, **arr_kw), zorder=7)
ax3.text(2.6e-2, 3.2e1, 'CGM cools', color=c_scgm, **{**dir_kw, 'ha': 'right'})

# Simulation runs
# tCGM: n_CGM decreases from 1 to 1e-3 at v_w = 100 km/s
ax3.plot([1e0, 1e-3], [1e2, 1e2], color=c_tcgm, lw=2.2, zorder=4,
         solid_capstyle='butt')
ax3.annotate('', xy=(10**-1.75, 1e2), xytext=(10**-1.25, 1e2),
             arrowprops=dict(arrowstyle='-|>', color=c_tcgm, lw=0,
                             mutation_scale=22), zorder=5)
ax3.plot(1e0, 1e2, marker='o', ms=7, mfc='white', mec=c_tcgm, mew=1.8,
         ls='None', zorder=5)
# uCGM_vw2: (log n, log v_w) = (-3, 2)
ax3.plot(1e-3, 1e2, marker='D', ms=9, mfc=c_ucgm, mec='white', mew=1.5,
         ls='None', zorder=6)

ax3.set_xlabel(r'$\mathrm{n_{CGM}\ (cm^{-3})}$', fontsize=18)
ax3.set_ylabel(r'$\mathrm{v_w\ (km\,s^{-1})}$', fontsize=18)
ax3.tick_params(axis='both', which='major', labelsize=15, pad=7)

ax3.legend(handles=[
    Line2D([0], [0], color=c_sw,   lw=1.6, label='Shocked wind'),
    Line2D([0], [0], color=c_scgm, lw=1.4, ls=(0, (5, 3)), label='Shocked CGM'),
    Line2D([0], [0], marker='D', ls='None', ms=8, mfc=c_ucgm, mec='white',
           label=r'$\mathtt{uCGM\_vw2}$'),
    Line2D([0], [0], color=c_tcgm, lw=2.2, marker='<', ms=7,
           label=r'$\mathtt{tCGM}$'),
    Patch(facecolor=both_colors[1], edgecolor='none', label='Wind + CGM cool'),
], loc='lower right', fontsize=12, frameon=True, framealpha=0.95,
   edgecolor='0.7', fancybox=False)

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "tcool_sw_scgm_runs.png"), dpi=300)
plt.close(fig3)
