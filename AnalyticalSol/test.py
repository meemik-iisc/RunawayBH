import numpy as np
from scipy.optimize import brentq

#Constants
G       = 6.67e-8       #cm^3/g/s^2
kB      = 1.38e-16      #erg/K
mp      = 1.67e-24      #g
kpc     = 3.086e21      #cm
M_sun   = 1.989e33      #g
mu      = 1.0
s_yr    = 3.154e7       #s
K       = 5.87e26       #g^(-2/3)cm^4s^(-2)
gamma   = 5.0/3.0       #Polytropic index
gm1     = gamma-1
c       = 3e10          #cm/s


#Parameters
M_bh    = 2e7*M_sun     #g
v_bh    = 1000*1e5      #cm/s
M_b     = 2e5*M_sun     #g
n_cgm   = 1e-3          #cm^-3
rho_cgm = n_cgm*mp      #g/cm^3
T_cgm   = 1e6           #K
cs_cgm  = np.sqrt(kB*T_cgm/(mu*mp))

# v_wind  = 1e4*1e5       #cm/s
t_sys   = 70*1e6*s_yr   #s
eta     = 0.1

rho_vir = mp            #g/cm^3
r_vir   = kpc           #cm
epsilon = 0.1*r_vir     #cm
#Calculate Polytropic Constant
K       = ((rho_cgm*(cs_cgm**2))/(rho_vir**(gamma)))
R_b     = ((3*G*M_bh*M_b)/(4*np.pi*n_cgm*T_cgm*kB))**(1/4)
rho_b   = 3*M_b/(4*np.pi*R_b**3)
T_b     = (G*M_bh*mu*mp)/(kB*R_b)


# print(R_b/kpc)
# print(rho_b/mp)
# print(T_b)
# print(rho_b*T_b, rho_cgm*T_cgm)

# def Phi(r):
#     return -G*M_bh/np.sqrt(r**2+epsilon**2)

# def Rho(r):
#     term1 = -1*gm1*Phi(r_vir)/(K*gamma)
#     term2 = -1*gm1*Phi(r)/(K*gamma)
#     return rho_vir-term1**(1/gm1)+term2**(1/gm1)

# def Pres(r):
#     return K*(Rho(r)**gamma)

# P_ram = rho_cgm*(v_bh**2)
# R0 = brentq(lambda r: Pres(r)-P_ram-(rho_cgm*(cs_cgm**2)), 1e-3*kpc, kpc)
# R1 = brentq(lambda r: Pres(r)-(rho_cgm*(cs_cgm**2)), 1e-3*kpc, kpc)
# print(f"Standoff distance = {R0/kpc:.2e} kpc")
# print(R1/kpc)

def calculate_stuff(v_wind):
    M_dot_w     = ((v_bh**2)*4*np.pi*(kpc**2)*rho_cgm)/(v_wind)
    # M_dot_w     = ((cs_cgm**2)*4*np.pi*(kpc**2)*rho_cgm)/(v_wind)
    # M_w         = M_dot_w*t_sys
    # E_wind      = 0.5*M_w*(v_wind**2)
    E_ke_bh     = 0.5*M_bh*(v_bh**2)
    M_dot_bh    = (0.5*M_dot_w*(v_wind**2))/(eta*c**2)
    M_cgm       = (rho_cgm*4*np.pi*(kpc**2)*v_bh*t_sys)
    E_cgm       = 0.5*M_cgm*(v_bh**2)
    f_duty      = (2*E_cgm)/(M_dot_w*t_sys*(v_wind**2))
    M_w         = f_duty*M_dot_w*t_sys
    E_wind      = 0.5*M_w*(v_wind**2)
     
    # if E_wind<=E_ke_bh:
    print("Energy Threshold passed")
    # M_dot_w     = ((cs_cgm**2)*4*np.pi*(kpc**2)*rho_cgm)/(v_wind)
    print(f"Wind velocity   = {v_wind:.2e}cm/s      \t = {v_wind/1e5:.2e} km/s")
    print(f"M dot wind      = {M_dot_w:.2e}g/s      \t = {M_dot_w/(M_sun/s_yr):.2e} M_sun/yr")
    print(f"M wind          = {M_w:.2e}g            \t = {M_w/(M_sun):.2e} M_sun")
    print(f"BH Accretion    = {M_dot_bh:.2e}g/s     \t = {M_dot_bh/(M_sun/s_yr):.2e} M_sun/yr")
    print(f"CGM mass swept  = {M_cgm:.2e}g          \t = {M_cgm/(M_sun):.2e} M_sun")
    print(f"Outflow KE      = {E_wind:.2e}erg")
    print(f"Black hole KE   = {E_ke_bh:.2e}erg")
    print(f"CGM KE          = {E_cgm:.2e}erg")
    print(f"Duty cycle      = {f_duty:.2e}")
# print(v_bh)
# for v_w in np.logspace(7.0,9.0,10):
#     print(v_w)
#     calculate_stuff(v_w)
calculate_stuff(1e8)