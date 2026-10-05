import numpy as np
import os
import math
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid
# Save figure in same directory as script
script_dir = os.path.dirname(os.path.abspath(__file__))
def n_cgm(t):
    # return n_cgm0*np.exp(-t/tau)
    """
    Two-phase ISM+CGM density as a function of time, in mp/cc (~cm^-3).

    ISM: n ~ n_ISM * exp(-L/H_ISM) for L << L_tr
    CGM: n ~ n_CGM * exp(-(L-L_tr)/H_CGM) for L >> L_tr
    Smooth transition around L_tr with width dL_kpc.
    Height is L(t) = v_bh * t, using constant v_bh0_kms.
    """
    # Parameters (tune these as you like)
    n_ISM_mp_cc   = 10.0      # ISM density at midplane [mp/cc]
    n_CGM_mp_cc   = 1e-3     # CGM density far out [mp/cc]
    H_ISM_kpc     = 10.0     # ISM scale height [kpc]
    H_CGM_kpc     = 100.0    # CGM scale height [kpc]
    L_tr_kpc      = 10.0     # transition height [kpc]
    dL_kpc        = 1.0      # transition width [kpc]

    # Ensure array-like
    t_arr = np.asarray(t, dtype=float)
    # BH velocity as function of time in cm/s
    v_bh_cgs = v_bh(t_arr)              # cm/s
    # Convert velocity to kpc/Myr
    v_bh_kpc_per_Myr = v_bh_cgs * s_Myr / kpc

    # Height L(t) via cumulative integral: L_i = ∫_0^{t_i} v(t') dt'
    # t_arr is in Myr, so dt_Myr; v_bh_kpc_per_Myr * dt_Myr -> kpc
    if t_arr.ndim == 0:
        # scalar: do a tiny integral analytically (approximate)
        L_kpc = v_bh_kpc_per_Myr * t_arr
    else:
        # array: use cumulative trapezoidal integration
        # ensure t_arr sorted
        order = np.argsort(t_arr)
        t_sorted = t_arr[order]
        v_sorted = v_bh_kpc_per_Myr[order]
        L_sorted = np.concatenate([[0.0],cumulative_trapezoid(v_sorted, t_sorted)])
        # map back to original order
        L_kpc = np.empty_like(t_arr)
        L_kpc[order] = L_sorted

    # ISM and CGM components as functions of height
    n_ISM_L = n_ISM_mp_cc * np.exp(-L_kpc / H_ISM_kpc)
    n_CGM_L = n_CGM_mp_cc * np.exp(-(L_kpc - L_tr_kpc) / H_CGM_kpc)

    # Smooth transition weight (0 -> ISM, 1 -> CGM)
    w = 0.5 * (1.0 + np.tanh((L_kpc - L_tr_kpc) / dL_kpc))

    # Blend ISM and CGM
    n_total = (1.0 - w) * n_ISM_L + w * n_CGM_L

    return n_total
def M_dot(t):
    rho_cgs = n_cgm(t)*mu*mp
    M_dot_cgs = (4*np.pi*R_0(t)**2*v_bh(t)**2*rho_cgs)/(v_w(t))
    return M_dot_cgs/(Msun/s_yr)
def M_total(t_array):
    """
    Total outflow mass between t_array[0] and t_array[-1], in Msun.
    M_dot(t) must return Msun/yr, t_array in Myr.
    """
    mdot_vals = M_dot(t_array)  # Msun/yr
    M_cum = np.zeros_like(t_array)
    for i in range(1, len(t_array)):
        # trapezoid between t[i-1] and t[i]
        dt_yr = (t_array[i] - t_array[i-1]) * 1e6
        M_seg = 0.5 * (mdot_vals[i] + mdot_vals[i-1]) * dt_yr
        M_cum[i] = M_cum[i-1] + M_seg
    return M_cum

def v_w(t):
    # return 100*kms*np.ones_like(t)
    return vw0*np.exp(-t/tau)

def v_bh(t):
    # return v_bh0*np.exp(-t/tau_vbh)
    return 1000*kms*np.ones_like(t)

def R_0(t):
    return R0*np.exp(t/tau_r)
    # return 1*kpc*np.ones_like(t)

def tcool_scgm(t):
    T_scgm          = (3*mu*mp*(v_bh(t)**2))/(16*kB)   
    return(9*mu*mp*(v_bh(t)**2))/(128*n_cgm(t)*cool_lambda(T_scgm))

def tcool_sw(t):
    T_sw            = (3*mu*mp*(v_w(t)**2))/(16*kB)
    return (9*mu*mp*(v_w(t)**4))/(128*n_cgm(t)*(v_bh(t)**2)*cool_lambda(T_sw))

# ============================================================================
# COOLING TABLE
# ============================================================================
cool_dir = os.path.join(script_dir,"cooltable.dat")
data = np.loadtxt(cool_dir)
T_tab = data[:, 0]
Lambda_tab = data[:, 1]
Lambda_interp = interp1d(T_tab, Lambda_tab, kind='linear', 
                         bounds_error=False, fill_value='extrapolate')

def cool_lambda(T):
    """Safe cooling function with clamping"""
    T_safe = np.clip(T, T_tab.min() * 0.95, T_tab.max() * 1.05)
    return Lambda_interp(T_safe)
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
n_cgm0 = 1
M_inf  = 1e7    #Msun
M_0    = 1e6    #Msun
tau    = 10    #Myr
R0     = 0.1*kpc  #cm
tau_r  = 20    #Myr
v_bh0  = 2000*kms
tau_vbh = 50    #Myr
vw0    = 1e4*kms
t = np.linspace(0,60,1000) #Myr
# Mdot_0 = ((4*np.pi*R0**2*v_bh(0)**2*n_cgm(0)*mu*mp)/(v_w(0)))/(Msun/s_yr)
# print(f"{Mdot_0:.2e}")
# print(M_dot(100))

panels = {
    "Mdot": {
        "yfunc": lambda t: M_dot(t),
        "ylabel": r"$\dot{M}_{\rm w}\ [\mathrm{M_\odot\ yr^{-1}}]$",
        "title": "Mass Outflow Rate",
        "color": "red",
        "yscale": "log",
    },
    "Mtotal": {
        "yfunc": lambda t: M_total(t),
        "ylabel": r"$M_{\rm out}\ [\mathrm{M_\odot}]$",
        "title": "Total Outflow Mass",
        "color": "blue",
        "yscale": "log",
    },
    "vw": {
        "yfunc": lambda t: v_w(t)/kms,
        "ylabel": r"$v_{\rm w}\ [\mathrm{km\ s^{-1}}]$",
        "title": "Wind Velocity",
        "color": "green",
        "yscale": "log",
    },
    "ncgm": {
        "yfunc": lambda t: n_cgm(t),
        "ylabel": r"$n_{\rm cgm}\ [\mathrm{cm^{-3}}]$",
        "title": "CGM Density",
        "color": "purple",
        "yscale": "log",
    },
    "vbh": {
        "yfunc": lambda t: v_bh(t)/kms,
        "ylabel": r"$v_{\rm bh}\ [\mathrm{km\ s^{-1}}]$",
        "title": "Black Hole Velocity",
        "color": "orange",
        "yscale": "log",
    },
    "R0": {
        "yfunc": lambda t: R_0(t)/kpc,
        "ylabel": r"$R_0\ [\mathrm{kpc}]$",
        "title": "Standoff Radius",
        "color": "brown",
        "yscale": "log",
    },
    "tcool_scgm": {
        "yfunc": lambda t: tcool_scgm(t)/s_Myr,
        "ylabel": r"$t_{\rm cool}\ [\mathrm{Myr}]$",
        "title": "Shocked CGM Cooling Time",
        "color": "magenta",
        "yscale": "log",
    },
    "tcool_sw": {
        "yfunc": lambda t: tcool_sw(t)/s_Myr,
        "ylabel": r"$t_{\rm cool}\ [\mathrm{Myr}]$",
        "title": "Shocked Wind Cooling Time",
        "color": "cyan",
        "yscale": "log",
    },
}

# Order in which to plot panels
plot_order = ["ncgm","vw","R0","vbh","Mdot", "Mtotal", "tcool_scgm", "tcool_sw"]

n_panels = len(plot_order)
n_rows = 4
n_cols = math.ceil(n_panels / n_rows)

fig, ax = plt.subplots(n_rows, n_cols, figsize=(6*n_cols, 4*n_rows))

# If n_rows or n_cols is 1, ax may not be 2D; normalize
ax = np.atleast_2d(ax)

for idx, key in enumerate(plot_order):
    row = idx // n_cols
    col = idx % n_cols

    panel = panels[key]
    y = panel["yfunc"](t)
    color = panel.get("color", f"C{idx}")
    yscale = panel.get("yscale", "linear")

    ax_rc = ax[row, col]
    ax_rc.plot(t, y, color=color, lw=2)
    ax_rc.set_yscale(yscale)
    ax_rc.set_xlabel("Time (Myr)")
    ax_rc.set_ylabel(panel["ylabel"])
    ax_rc.set_title(panel["title"])
    ax_rc.grid(True, alpha=0.3)

# If there are empty slots (e.g. last row not full), hide them
for idx in range(n_panels, n_rows*n_cols):
    row = idx // n_cols
    col = idx % n_cols
    ax[row, col].axis("off")

fig.suptitle("Wind–CGM Evolution", fontsize=16, y=0.98)
plt.tight_layout(rect=[0, 0, 1, 0.96])


out_path = os.path.join(script_dir, "wind_cgm_evolution.png")
fig.savefig(out_path, dpi=200)
print(f"Figure saved to {out_path}")
# plt.plot(t, v_w(t)/kms, label = "v_w")
# # plt.yscale('log')
# # plt.plot(t, M_total(t))
# plt.show()