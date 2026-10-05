import numpy as np

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

#Parameters
n_amb   = 0.1           #cm^-3
d_amb   = n_amb*mu*mp
T_amb   = 1e6           #k
P_amb   = n_amb*kB*T_amb#dyne/cm^2
cs_amb  = np.sqrt(kB*T_amb/(mu*mp))
#Code Units
L_code  = 1.0*pc        #cm
v_code  = 1e2*kms       #cm/s
rho_code= 1.0*mp        #g/cm^3

M_code  = rho_code*L_code**3
t_code  = L_code/v_code
P_code  = rho_code*v_code**2
T_code  = (mu*mp/kB)*(v_code**2)



print("="*10,"Code Units","="*90)
print(f"Code Length         = {L_code:.2e} cm \t\t = {L_code/kpc:.2e} kpc")
print(f"Code Mass           = {M_code:.2e} g \t\t = {M_code/Msun:.2e} Msun")
print(f"Code Time           = {t_code:.2e} s \t\t = {t_code/(1e6*s_yr):.2e} Myr")
print(f"Code Density        = {rho_code:.2e} g/cm^3 \t\t = {rho_code/mp:.2e} mp/cm^3")
print(f"Code Velocity       = {v_code:.2e} cm/s \t\t = {v_code/1e5:.2e} km/s")
print(f"Code Temperature    = {T_code:.2e} K")
print(f"Code Pressure       = {P_code:.2e} dyne/cm^2")
# print(f"Code Pressure1    = {P_code1:.2e} dyne/cm^2")
print("="*10,"Parameters","="*90)
print(f"Ambient density     = {d_amb:.2e} g/cm^3\t\t = {d_amb/rho_code:.2e} rho_code")
print(f"Ambient Temperature = {T_amb:.2e} K \t\t = {T_amb/T_code:.2e} T_code")
print(f"Ambient Pressure    = {P_amb:.2e} dyne/cm^2 \t = {P_amb/P_code:.2e} P_code")
print(f"Ambient sound speed = {cs_amb:.2e} cm/s \t\t = {cs_amb/v_code:.2e} v_code")



print("="*112)

# P_vir   = K_code*((rho_vir/rho_code)**gamma)*P_code
# print(f"Virial Pressure     = {P_vir:.2e} dyne/cm^2 \t = {P_vir/P_code:.2e} P_code")