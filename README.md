# RunawayBH

Analytical models, analysis scripts and figures for studying a **runaway
(recoiling) supermassive black hole** moving supersonically through the
circumgalactic medium (CGM), and the wind/outflow it drives. The hydrodynamic
simulations themselves are run with [AthenaK](https://github.com/IAS-Astrophysics/athenak);
this repository holds the analysis side of the project.

Fiducial setup used throughout the scripts (CGS units internally):

| Quantity | Value |
|---|---|
| Black hole mass `M_bh` | 2 × 10⁷ M☉ |
| Black hole velocity `v_bh` | 1000 km/s |
| Bondi radius | `R_B = 2 G M_bh / (v_bh² + c_s²)` |
| CGM density `n_cgm` | 10⁻³ – 10⁻¹ cm⁻³ |
| CGM temperature `T_cgm` | 10⁶ K |
| Adiabatic index `γ` | 5/3, `μ = 1` |

## Repository layout

```
RunawayBH/
├── AnalyticalSol/          Analytical / semi-analytical models and their plots
│   ├── bubble/             Hydrostatic bubble around the BH (isothermal & polytropic)
│   ├── shock_cooling/      Cooling times of shocked CGM and shocked wind
│   ├── spherical_wind/     Spherical wind expansion r(t), v(t), termination shock
│   ├── Mdot_params/        Wind mass-injection rate vs. time and path length
│   ├── cooltable.dat       Cooling table: T [K], Λ(T) [erg cm³ s⁻¹]
│   └── *.py                Cooling maps, timescale sweeps, unit conversions
├── plots_for_paper/        Final figures for the paper
│   ├── baseline/
│   └── outflows/{tCGM,uCGM}/
├── sim_data/               Simulation runs and outputs (git-ignored, local only)
└── check_git_size.sh       Checks staged file sizes against GitHub limits
```

### `AnalyticalSol/`

| Script | What it does |
|---|---|
| `bubble/isothermal_sol.py` | Isothermal hydrostatic profile around the BH (softened potential) |
| `bubble/polytropic_sol.py` | Adiabatic/polytropic profiles, entropy and dynamical timescales |
| `bubble/compute_bubblemass.py`, `compute_bubblemass2.py` | Bubble mass and density, pressure, temperature and enclosed-mass profiles for an AthenaK initial condition |
| `shock_cooling/shock_cooling_plot.py` | Cooling time of shocked CGM and shocked wind vs. `n_cgm` and `v_wind`, with `Ṁ` contours |
| `spherical_wind/spherical_wind.py` | Time evolution of a spherical wind: radius, velocity and termination-shock radius |
| `Mdot_params/Mdot_vs_time.py` | Wind mass-injection rate, total injected mass, standoff radius and cooling times vs. time along the BH trajectory |
| `Mdot_params/Mdot_vs_L.py` | The same quantities vs. path length `L` through a two-phase ISM→CGM density profile |
| `Mass_inj_vs_time.py` | Combined wind–CGM evolution panel plot |
| `cooling_curve.py` | Cooling-time map in the (n, T) plane |
| `test.py`, `test1.py` | Hydro vs. cooling timestep comparison (subcycling); `test1.py` ports AthenaK's `ISMCoolFn()` |
| `calculate_alpha.py` | Power-law indices α for time-varying CGM density and wind velocity, `∝ (1 + t/τ)^α` |
| `unit_conversion.py`, `unit_conversion2.py` | Convert physical parameters to AthenaK code units |

The ISM→CGM density model in `Mdot_params` blends two exponentials with a
`tanh` transition:

```
n(L) = (1 − w) n_ISM e^(−L/H_ISM) + w n_CGM e^(−(L − L_tr)/H_CGM),
w    = ½ [1 + tanh((L − L_tr)/ΔL)]
```

### `sim_data/` (not tracked)

Simulation outputs are too large for git and are listed in `.gitignore`. Locally
they are organised by setup:

| Folder | Setup |
|---|---|
| `bondiAccr/` | Bondi accretion tests at ambient T = 10⁴, 10⁷, 10⁹ K |
| `adiaBondi/` | Adiabatic moving-BH runs (10 Myr, 100 Myr) |
| `radBondi/` | The same with radiative cooling |
| `wind_bondi/` | Wind + Bondi test run |
| `windOutflowRadBondiTimeCGM/` | Wind outflow into a time-varying CGM |
| `uCGM/` | Outflows into a uniform CGM |
| `tCGM/` | Outflows into a time-varying CGM (with and without cooling) |

## Requirements

Python 3 with:

```bash
pip install numpy scipy matplotlib
```

## Usage

Each script runs on its own: edit the parameter block at the top, then run it.
Figures are saved next to the script.

```bash
python AnalyticalSol/shock_cooling/shock_cooling_plot.py
python AnalyticalSol/Mdot_params/Mdot_vs_L.py
```

`AnalyticalSol/cooling_curve.py` reads the cooling table from an absolute path;
change `cool_dir` if your clone lives elsewhere.

Before committing new outputs, run `./check_git_size.sh` to catch files that
exceed GitHub's size limits.
