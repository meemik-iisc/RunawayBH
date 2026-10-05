import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch
import os
output_dir = os.path.dirname(os.path.abspath(__file__))
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

#Flags
plot_models = False
fill_regions = True

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
R0      = 1*kpc
R01     = 0.5*kpc
t_sys   = 70*1e6*s_yr   #s
eta     = 0.1
# M_dot_w = ((v_bh**2)*4*np.pi*(R0**2)*rho_cgm)/(v_wind)
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
#Calculate Temperature of Shocked Wind

def t_cool_sw(ncgm, vw):
    T_sw            = (3*mu*mp*(vw**2))/(16*kB)
    tcool_sw        = (9*mu*mp*(vw**4))/(128*ncgm*(v_bh**2)*cool_lambda(T_sw))
    return tcool_sw
def T_sw(vw):
    return (3*mu*mp*(vw**2))/(16*kB)

def M_dot_w(ncgm, vw):
    rho_cgm = ncgm*mu*mp
    M_dot_w = ((v_bh**2)*4*np.pi*(R0**2)*rho_cgm)/(vw)
    return M_dot_w

def M_dot_w1(ncgm, vw):
    rho_cgm = ncgm*mu*mp
    M_dot_w = ((v_bh**2)*4*np.pi*(R01**2)*rho_cgm)/(vw)
    return M_dot_w

def t_cool_mix(ncgm, vw):
    T_mix           = (3*mu*mp*(vw*v_bh))/(16*kB)
    P_mix           = (3*ncgm*mu*mp*(v_bh**2))/4
    tcool_mix       = (3*(kB*T_mix)**2)/(2*P_mix*cool_lambda(T_mix))
    return tcool_mix

print("="*30,"Parameters","="*40)
print(f"CGM Density             = {rho_cgm/mp:.2e} mp/cm^3 \t\t = {rho_cgm/rho_code:.2e} rho_code")
print(f"CGM Temp                = {T_cgm:.2e} K \t\t = {T_cgm/T_code:.2e} T_code")
print(f"Outflow Velocity        = {v_wind:.2e} cm/s \t\t = {v_wind/v_code:.2e} v_code")

n_cgm   = np.logspace(-4, 3, 100)
vw      = np.array([1e2, 1e3, 1e4])      #km/s
vw_cgs  = vw*kms

# ============================================================================
# MODEL DICTIONARY
# ============================================================================
models = {
    'uCGM3_vw2': {
        'x': [1e-3], 'y': [1e2],
        'color': 'orange', 'marker': 'D', 'markersize': 10,
        'markerfacecolor': 'orange', 'markeredgecolor': 'black', 'markeredgewidth': 1.5,
        'linestyle': None, 'linewidth': None
    },
    'uCGM0_vw2': {
        'x': [1], 'y': [1e2],
        'color': 'green', 'marker': 'D', 'markersize': 10,
        'markerfacecolor': 'green', 'markeredgecolor': 'black', 'markeredgewidth': 1.5,
        'linestyle': None, 'linewidth': None
    },
    'tCGM_vw2': {
        'x': [1e-3, 1], 'y': [1e2, 1e2],
        'color': 'green', 'linewidth': 4, 'linestyle': '-',
        'marker': '<', 'markersize': 10,
        'markerfacecolor': 'green', 'markeredgecolor': 'black', 'markeredgewidth': 3,
        'extra_marker_x': 3e-2, 'extra_marker_y': 1e2
    },
    'tCGM01_tvw50': {
        'x': [1e-1, 1e-3], 'y': [5e1, 5e2],
        'color': 'red', 'linewidth': 4, 'linestyle': '-',
        'marker': (3, 0, -63), 'markersize': 12,
        'markerfacecolor': 'red', 'markeredgecolor': 'black', 'markeredgewidth': 2,
        'extra_marker_x': np.sqrt(1e-1*1e-3), 'extra_marker_y': np.sqrt(5e1*5e2)
    },
    'tCGM001_tvw10': {
        'x': [1e-2, 1e-4], 'y': [1e1, 1e2],
        'color': 'blue', 'linewidth': 4, 'linestyle': '-',
        'marker': (3, 0, -63), 'markersize': 12,
        'markerfacecolor': 'blue', 'markeredgecolor': 'black', 'markeredgewidth': 2,
        'extra_marker_x': np.sqrt(1e-2*1e-4), 'extra_marker_y': np.sqrt(1e1*1e2)
    },
    'tCGM10_tvw50': {
        'x': [1e1, 1e-3], 'y': [5e1, 5e2],
        'color': 'brown', 'linewidth': 4, 'linestyle': '-',
        'marker': (3, 0, -45), 'markersize': 12,
        'markerfacecolor': 'brown', 'markeredgecolor': 'black', 'markeredgewidth': 2,
        'extra_marker_x': np.sqrt(1e1*1e-3), 'extra_marker_y': np.sqrt(5e1*5e2)
    },
    'tCGM10_tvw100': {
        'x': [1e1, 1e-1], 'y': [1e2, 1e3],
        'color': 'purple', 'linewidth': 4, 'linestyle': '-',
        'marker': (3, 0, -63), 'markersize': 12,
        'markerfacecolor': 'purple', 'markeredgecolor': 'black', 'markeredgewidth': 2,
        'extra_marker_x': np.sqrt(1e1*1e-1), 'extra_marker_y': np.sqrt(1e2*1e3)
    }
}



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
# ax1.plot(np.log10(n_cgm), np.log10(t_cool_sw(n_cgm, vw[0])/s_yr),  linewidth=2, linestyle='solid', markersize=2, color='blue', label="Shocked Wind (1e2 km/s)")
# ax1.plot(np.log10(n_cgm), np.log10(t_cool_sw(n_cgm, vw[1])/s_yr),  linewidth=2, linestyle='dashed', markersize=2, color='blue', label="Shocked Wind (1e3 km/s)")
# ax1.plot(np.log10(n_cgm), np.log10(t_cool_sw(n_cgm, vw[2])/s_yr),  linewidth=2, linestyle='dotted', markersize=2, color='blue', label="Shocked Wind (1e4 km/s)")
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
# Plot 3: contours in (n_cgm, v_w) plane + shaded regions
# ============================================================
v_w_kms_2d = np.logspace(1, 4, 300)
v_w_2d = v_w_kms_2d * kms
n_cgm_2d = np.logspace(-4, 2, 300)

NCGM, VW = np.meshgrid(n_cgm_2d, v_w_2d)

TCOOL_SCGM = t_cool_scgm(NCGM) / s_yr
TCOOL_SW   = t_cool_sw(NCGM, VW) / s_yr
TCOOL_MIX  = t_cool_mix(NCGM, VW) / s_yr

Mdot_Msun_yr = M_dot_w(NCGM, VW)*s_yr/Msun
Mdot1_Msun_yr = M_dot_w1(NCGM, VW)*s_yr/Msun

time_Myr  = np.array([1, 10, 70])
levels_yr = time_Myr * 1e6
level_colors = ['red', 'green', 'blue']

Mdot_levels = np.array([1e0, 1e1, 1e2])
Mdot_levels_cgs = Mdot_levels * Msun / s_yr
Mdot_colors = ['brown', 'purple', 'magenta']


fig3, ax3 = plt.subplots(figsize=(10, 6))

# contours
cs1 = ax3.contour(
    # np.log10(n_cgm_2d),
    # np.log10(v_w_kms_2d),
    n_cgm_2d,
    v_w_kms_2d,
    TCOOL_SCGM,
    levels=levels_yr,
    colors=level_colors,
    linewidths=2,
    linestyles='dashed'
)

cs2 = ax3.contour(
    # np.log10(n_cgm_2d),
    # np.log10(v_w_kms_2d),
    n_cgm_2d,
    v_w_kms_2d,
    TCOOL_SW,
    levels=levels_yr,
    colors=level_colors,
    linewidths=2,
    linestyles='solid'
)
# cs3 = ax3.contour(
#     # np.log10(n_cgm_2d),
#     # np.log10(v_w_kms_2d),
#     n_cgm_2d,
#     v_w_kms_2d,
#     TCOOL_MIX,
#     levels=levels_yr,
#     colors=level_colors,
#     linewidths=2,
#     linestyles='dashed'
# )

cs4 = ax3.contour(
    # np.log10(n_cgm_2d),
    # np.log10(v_w_kms_2d),
    n_cgm_2d,
    v_w_kms_2d,
    Mdot_Msun_yr,
    levels=Mdot_levels,
    colors=Mdot_colors,
    linewidths=2,
    linestyles='solid'
)
# cs2 = ax3.contour(
#     # np.log10(n_cgm_2d),
#     # np.log10(v_w_kms_2d),
#     n_cgm_2d,
#     v_w_kms_2d,
#     Mdot1_Msun_yr,
#     levels=Mdot_levels,
#     colors=Mdot_colors,
#     linewidths=2,
#     linestyles='dashed'
# )

if plot_models:
    
    # Loop through models to plot
    legend_handles = []

    for label, params in models.items():
        # Plot the line (if x has 2 points)
        if len(params['x']) == 2:
            ax3.plot(params['x'], params['y'], 
                    color=params['color'], 
                    linewidth=params['linewidth'], 
                    linestyle=params['linestyle'])
        
        # Plot the marker
        if len(params['x']) == 1:
            # Single point marker
            ax3.plot(params['x'][0], params['y'][0],
                    marker=params['marker'],
                    markersize=params['markersize'],
                    markerfacecolor=params['markerfacecolor'],
                    markeredgecolor=params['markeredgecolor'],
                    markeredgewidth=params['markeredgewidth'],
                    linestyle=params['linestyle'],
                    color=params['color'])
        else:
            # Line with arrowhead marker at end
            ax3.plot(params['x'][1], params['y'][1],
                    linestyle='None')
        
        # Plot extra marker if exists
        if 'extra_marker_x' in params:
            ax3.plot(params['extra_marker_x'], params['extra_marker_y'],
                    marker=params['marker'],
                    markersize=params['markersize'],
                    markerfacecolor=params['markerfacecolor'],
                    markeredgecolor=params['markeredgecolor'],
                    markeredgewidth=params['markeredgewidth'],
                    linestyle='None')
        
        # Create legend handle
        if len(params['x']) == 1 or params['marker'] == 'D':
            # Point marker legend
            handle = Line2D([0], [0], 
                        color=params['color'], 
                        marker=params['marker'], 
                        linestyle='None', 
                        markeredgecolor=params['markeredgecolor'], 
                        markerfacecolor=params['markerfacecolor'], 
                        markersize=8, 
                        label=label)
        else:
            # Line with arrow legend
            handle = Line2D([0], [0], 
                        color=params['color'], 
                        marker="<", 
                        linestyle='solid', 
                        lw=params['linewidth'], 
                        markeredgecolor=params['markeredgecolor'], 
                        markerfacecolor=params['markerfacecolor'], 
                        markersize=10, 
                        label=label)
        
        legend_handles.append(handle)
        
        legend_right = ax3.legend(
            handles=legend_handles,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            frameon=True,
            title='Models',
            fontsize=11
        )
    ax3.add_artist(legend_right)



#Fill Regions
if fill_regions:
    yvals = np.log10(v_w_kms_2d)
    for t_yr, col in zip(levels_yr, level_colors):
        region =(TCOOL_SW <= t_yr) & (TCOOL_SCGM <= t_yr)
        # region =(TCOOL_SW <= t_yr)
        

        Zfill = np.where(region, 1.0, 0.0)

        ax3.contourf(
            # np.log10(n_cgm_2d),
            # np.log10(v_w_kms_2d),
            n_cgm_2d,
            v_w_kms_2d,
            Zfill,
            levels=[0.5, 1.5],
            colors=[col],
            alpha=0.12,
            zorder=1
        )
    for mdot, mdot_col in zip(Mdot_levels, Mdot_colors):
        region1 = Mdot_Msun_yr < mdot

        Zfill = np.where(region1, 1.0, 0.0)

        ax3.contourf(
            # np.log10(n_cgm_2d),
            # np.log10(v_w_kms_2d),
            n_cgm_2d,
            v_w_kms_2d,
            Zfill,
            levels=[0.5, 1.5],
            colors=[mdot_col],
            alpha=0.12,
            zorder=1
        )
    
ax3.set_xlabel(r'$\mathbf{n_{\rm CGM}\,(cm^{-3})}$', fontsize=20, fontweight='bold')
# ax3.set_xlim([-3,3])
ax3.set_ylabel(r'$\mathbf{v_{\rm w}\,(km/s)}$', fontsize=20, fontweight='bold')
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.tick_params(axis='both', which='major', labelsize=16)
ax3.minorticks_on()
ax3.grid(True, which='both', ls='-', alpha=0.3)
ax3.grid(which='major', alpha=0.7, linewidth=1.0)
ax3.grid(which='minor', alpha=0.5, linewidth=0.8)    

legend_time = ax3.legend(
    handles=[
        Line2D([0], [0], color=level_colors[0], lw=2, linestyle='solid', label=r'$t_{cool} = 1 Myr$'),
        Line2D([0], [0], color=level_colors[1], lw=2, linestyle='solid', label=r'$t_{cool} = 10 Myr$'),
        Line2D([0], [0], color=level_colors[2], lw=2, linestyle='solid', label=r'$t_{cool} = 70 Myr$'),
        Line2D([0], [0], color='k', lw=2, linestyle='solid', label='Shocked Wind'),
        Line2D([0], [0], color='k', lw=2, linestyle='dashed', label='Shocked CGM'),
        # Line2D([0], [0], color='k', lw=2, linestyle='dashed', label='mixed gas'),
    ],
    loc='lower right',
    fontsize=14,
    frameon=True
)

legend_mdot = ax3.legend(
    handles=[
        Line2D([0], [0], color=Mdot_colors[0],    lw=2, linestyle='solid', label=r'$\dot{M}_{out} = 1 M_{\odot}/yr$'),
        Line2D([0], [0], color=Mdot_colors[1],    lw=2, linestyle='solid', label=r'$\dot{M}_{out} = 10 M_{\odot}/yr$'),
        Line2D([0], [0], color=Mdot_colors[2],   lw=2, linestyle='solid', label=r'$\dot{M}_{out} = 100 M_{\odot}/yr$'),
        # Line2D([0], [0], color='k', lw=2, linestyle='solid', label=r"$R_0$" f" = {R0/kpc:.1f}kpc"),
        # Line2D([0], [0], color='k', lw=2, linestyle='dashed', label=f"R0 = {R01/kpc:.1f} kpc"),
    ],
    loc='upper right',
    fontsize=14,
    frameon=True
)


# ax3.add_artist(legend_style)
ax3.add_artist(legend_time)
ax3.add_artist(legend_mdot)



plt.tight_layout()
if plot_models:
    plt.subplots_adjust(right=0.80)
    plt.savefig(os.path.join(output_dir, "tcool_Mdot_vw_vs_ncgm_model.png"),dpi=300)
else:
    plt.savefig(os.path.join(output_dir, "tcool_sw_scgm_mdot.png"),dpi=300)
plt.close(fig3)