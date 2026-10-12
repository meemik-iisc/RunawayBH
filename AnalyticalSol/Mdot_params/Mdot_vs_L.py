import numpy as np
import os
import math
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.interpolate import interp1d
from scipy.integrate import cumulative_trapezoid

script_dir = os.path.dirname(os.path.abspath(__file__))

# ── Global style ─────────────────────────────────────────────────────────────
font_size   = 18          # axis labels
tick_size   = 16          # tick labels
title_size  = font_size + 6
suptitle_size = font_size + 10

mpl.rcParams.update({
    "font.family":        "sans-serif",
    "axes.titlesize":     title_size,
    "axes.titleweight":   "bold",
    "axes.labelsize":     font_size,
    "axes.labelweight":   "bold",
    "xtick.labelsize":    tick_size,
    "ytick.labelsize":    tick_size,
    "xtick.major.width":  1.5,
    "ytick.major.width":  1.5,
    "xtick.minor.width":  1.0,
    "ytick.minor.width":  1.0,
    "xtick.major.size":   6,
    "ytick.major.size":   6,
    "xtick.minor.size":   3,
    "ytick.minor.size":   3,
    "xtick.direction":    "in",
    "ytick.direction":    "in",
    "axes.linewidth":     1.5,
    "lines.linewidth":    2.5,
    "grid.alpha":         0.25,
    "grid.linewidth":     0.8,
    "legend.fontsize":    tick_size - 2,
})

def bold_ticks(ax):
    """Make all tick labels bold."""
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontweight("bold")

# ── Physics / model functions ────────────────────────────────────────────────
def n_cgm(L_kpc):
    n_ISM_L = n_ISM_mp_cc * np.exp(-L_kpc / H_ISM_kpc)
    n_CGM_L = n_CGM_mp_cc * np.exp(-(L_kpc - L_tr_kpc) / H_CGM_kpc)
    w = 0.5 * (1.0 + np.tanh((L_kpc - L_tr_kpc) / dL_kpc))
    # return (1.0 - w) * n_ISM_L + w * n_CGM_L
    return n_ISM_L*np.exp(-L_kpc/10)

def v_bh(L_kpc):
    return vbh_f * np.ones_like(L_kpc)


def v_w(L_kpc):
    dL = 25
    L_tr = 10
    vw_0_L = vw_0 * np.exp(-L_kpc / (10*H_ISM_kpc))
    vw_f_L = vw_f* np.exp(-(L_kpc - L_tr) / (10*H_CGM_kpc))
    # w = 0.5 * (1.0 + np.tanh((L_kpc - L_tr) / dL))
    # return (1.0 - w) * vw_0_L + w * vw_f_L
    # return vw_f* np.exp(np.exp(-1*(L_kpc / dL)))
    # return vw_0* np.exp(-1*(L_kpc / dL))
    return vw_f * np.ones_like(L_kpc)

def M_dot(L_kpc):
    dL_rise = 10    # scale of rise
    dL_decay = 20  # scale of decay
    M_max  = 0.1
    M_min  = 0.01  # floor it decays TO
    
    # rise  = 1 - np.exp(-L_kpc / dL_rise)
    # decay = (M_max - M_min) * np.exp(-L_kpc / dL_decay) + M_min
    # Mdot  = rise * decay
    # return Mdot
    return 1e-3*np.exp(L_kpc/10)

def R_0(L_kpc):
    Mdot_cgs = M_dot(L_kpc) * Msun / s_yr
    R0 = np.sqrt((Mdot_cgs*v_w(L_kpc))/(4*np.pi*n_cgm(L_kpc)*mu*mp*v_bh(L_kpc)**2))
    # return R0 + ((Rf - R0) / 2) * (1 + np.tanh((L_kpc - L_tr) / dL))
    return R0
def M_total(L_array):
    mdot_vals = M_dot(L_array)
    M_cum = np.zeros_like(L_array)
    for i in range(1, len(L_array)):
        dt_yr = (L_array[i] - L_array[i - 1]) * 1e6
        M_cum[i] = M_cum[i - 1] + 0.5 * (mdot_vals[i] + mdot_vals[i - 1]) * dt_yr
    return M_cum

def tcool_scgm(L):
    T_scgm = (3 * mu * mp * v_bh(L)**2) / (16 * kB)
    return (9 * mu * mp * v_bh(L)**2) / (128 * n_cgm(L) * cool_lambda(T_scgm))
def tcool_scgm_og(L):
    ncgm_og = n_ISM_mp_cc*np.exp(-L/10)
    T_scgm = (3 * mu * mp * v_bh(L)**2) / (16 * kB)
    return (9 * mu * mp * v_bh(L)**2) / (128 * ncgm_og * cool_lambda(T_scgm))
def tcool_sw(L):
    T_sw = (3 * mu * mp * v_w(L)**2) / (16 * kB)
    return (9 * mu * mp * v_w(L)**4) / (128 * n_cgm(L) * v_bh(L)**2 * cool_lambda(T_sw))
def tcool_sw_og(L):
    vw_og = vw_f*np.ones_like(L)
    ncgm_og = n_ISM_mp_cc*np.exp(-L/10)
    T_sw = (3 * mu * mp * vw_og**2) / (16 * kB)
    return (9 * mu * mp * vw_og**4) / (128 * ncgm_og * v_bh(L)**2 * cool_lambda(T_sw))
# ── Cooling table ─────────────────────────────────────────────────────────────
cool_dir = os.path.join(script_dir, "cooltable.dat")
data = np.loadtxt(cool_dir)
T_tab, Lambda_tab = data[:, 0], data[:, 1]
Lambda_interp = interp1d(T_tab, Lambda_tab, kind='linear',
                         bounds_error=False, fill_value='extrapolate')

def cool_lambda(T):
    return Lambda_interp(np.clip(T, T_tab.min() * 0.95, T_tab.max() * 1.05))

# ── Constants ─────────────────────────────────────────────────────────────────
G, kB, mp   = 6.67e-8, 1.38e-16, 1.67e-24
kpc, pc     = 3.086e21, 3.086e18
Msun        = 1.989e33
mu          = 1.0
s_yr, s_Myr = 3.154e7, 3.154e13
gamma, gm1  = 5.0 / 3.0, 2.0 / 3.0
kms         = 1e5

# ── Parameters ────────────────────────────────────────────────────────────────
n_ISM_mp_cc = 10.0
n_CGM_mp_cc = 1e-3
vbh_0, vbh_f = 3000 * kms, 1000 * kms
vw_0,  vw_f  = 1e3  * kms, 1e2  * kms
H_ISM_kpc    = 10.0
H_CGM_kpc    = 100.0
L_tr_kpc     = 10.0
dL_kpc       = 1.0

t = L_kpc = np.linspace(0, 60, 1000)

# ── Panel definitions ─────────────────────────────────────────────────────────
panels = {
    "ncgm":      {"yfunc": lambda L: n_cgm(L),          "ylabel": r"$n_{\rm cgm}\ [\mathrm{cm^{-3}}]$",          "title": "CGM Density",               "color": "mediumpurple"},
    "vw":        {"yfunc": lambda L: v_w(L) / kms,       "ylabel": r"$v_{\rm w}\ [\mathrm{km\ s^{-1}}]$",         "title": "Wind Velocity",             "color": "seagreen"},
    "R0":        {"yfunc": lambda L: R_0(L) / kpc,       "ylabel": r"$R_0\ [\mathrm{kpc}]$",                      "title": "Standoff Radius",           "color": "saddlebrown"},
    "vbh":       {"yfunc": lambda L: v_bh(L) / kms,      "ylabel": r"$v_{\rm bh}\ [\mathrm{km\ s^{-1}}]$",        "title": "Black Hole Velocity",       "color": "darkorange"},
    "Mdot":      {"yfunc": lambda L: M_dot(L),            "ylabel": r"$\dot{M}_{\rm w}\ [\mathrm{M_\odot\ yr^{-1}}]$", "title": "Mass Outflow Rate",    "color": "crimson"},
    "Mtotal":    {"yfunc": lambda L: M_total(L),          "ylabel": r"$M_{\rm out}\ [\mathrm{M_\odot}]$",           "title": "Total Outflow Mass",        "color": "steelblue"},
    "tcool_scgm":{"yfunc": lambda L: tcool_scgm(L)/s_Myr,"ylabel": r"$t_{\rm cool}\ [\mathrm{Myr}]$",             "title": "Shocked CGM Cooling Time",  "color": "deeppink"},
    "tcool_sw":  {"yfunc": lambda L: tcool_sw(L) /s_Myr, "ylabel": r"$t_{\rm cool}\ [\mathrm{Myr}]$",             "title": "Shocked Wind Cooling Time", "color": "dodgerblue"},
}

OG_panels = {
    "ncgm":      {"yfunc": lambda L: 1*np.exp(-L/10),          "ylabel": r"$n_{\rm cgm}\ [\mathrm{cm^{-3}}]$",          "title": "CGM Density",               "color": "mediumpurple"},
    "vw":        {"yfunc": lambda L: vw_f*np.ones_like(L) / kms,       "ylabel": r"$v_{\rm w}\ [\mathrm{km\ s^{-1}}]$",         "title": "Wind Velocity",             "color": "seagreen"},
    "R0":        {"yfunc": lambda L: 1*np.ones_like(L),       "ylabel": r"$R_0\ [\mathrm{kpc}]$",                      "title": "Standoff Radius",           "color": "saddlebrown"},
    "vbh":       {"yfunc": lambda L: v_bh(L) / kms,      "ylabel": r"$v_{\rm bh}\ [\mathrm{km\ s^{-1}}]$",        "title": "Black Hole Velocity",       "color": "darkorange"},
    "Mdot":      {"yfunc": lambda L: 1*np.ones_like(L),            "ylabel": r"$\dot{M}_{\rm w}\ [\mathrm{M_\odot\ yr^{-1}}]$", "title": "Mass Outflow Rate",    "color": "crimson"},
    "Mtotal":    {"yfunc": lambda L: 1e6*L,          "ylabel": r"$M_{\rm out}\ [\mathrm{M_\odot}]$",           "title": "Total Outflow Mass",        "color": "steelblue"},
    "tcool_scgm":{"yfunc": lambda L: tcool_scgm_og(L)/s_Myr,"ylabel": r"$t_{\rm cool}\ [\mathrm{Myr}]$",             "title": "Shocked CGM Cooling Time",  "color": "deeppink"},
    "tcool_sw":  {"yfunc": lambda L: tcool_sw_og(L) /s_Myr, "ylabel": r"$t_{\rm cool}\ [\mathrm{Myr}]$",             "title": "Shocked Wind Cooling Time", "color": "dodgerblue"},
}
plot_order = ["ncgm", "Mdot","vw", "vbh", "R0", "Mtotal", "tcool_scgm", "tcool_sw"]
n_panels   = len(plot_order)
n_rows, n_cols = 2, math.ceil(n_panels / 2)

fig, ax = plt.subplots(n_rows, n_cols, figsize=(7 * n_cols, 6 * n_rows),
                       constrained_layout=True)
ax = np.atleast_2d(ax)

for idx, key in enumerate(plot_order):
    row, col = divmod(idx, n_cols)
    panel    = panels[key]
    color    = panel["color"]
    ax_rc    = ax[row, col]

    y    = panel["yfunc"](t)
    og_y = OG_panels[key]["yfunc"](t)

    ax_rc.plot(t, y,    color=color, lw=2.5, label="Model")
    ax_rc.plot(t, og_y, color=color, lw=2.5, ls="--", label="Reference")
    ax_rc.set_yscale("log")
    ax_rc.set_xlabel("Length (kpc)")
    ax_rc.set_ylabel(panel["ylabel"])
    ax_rc.set_title(panel["title"])
    ax_rc.grid(True, which="both")
    ax_rc.xaxis.set_minor_locator(mpl.ticker.AutoMinorLocator())

    # Bold tick labels
    bold_ticks(ax_rc)

    # Legend only on first panel
    if idx == 0:
        ax_rc.legend(framealpha=0.85)

# Hide empty slots
for idx in range(n_panels, n_rows * n_cols):
    ax[divmod(idx, n_cols)].axis("off")

fig.suptitle("Mass Injection vs Length", fontsize=suptitle_size, fontweight="bold", y=1.05)

out_path = os.path.join(script_dir, "Mass_inj_vs_length.png")
fig.savefig(out_path, dpi=200, bbox_inches="tight")
print(f"Figure saved to {out_path}")
print(f"Final M_out = {M_total(L_kpc)[-1]:.2e}")