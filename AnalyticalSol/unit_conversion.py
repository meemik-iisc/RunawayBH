import numpy as np

#Constants
G       = 6.67e-8       #cm^3/g/s^2
kB      = 1.38e-16      #erg/K
mp      = 1.67e-24      #g
kpc     = 3.086e21      #cm
Msun    = 1.989e33      #g
mu      = 1.0
s_yr    = 3.154e7       #s
# K       = 2.91e28       #g^(-2/3)cm^4s^(-2)
gamma   = 5.0/3.0       #Polytropic index
gm1     = gamma-1

#Parameters
#Parameters
v_bh    = 1000*1e5                  #cm/s
M_bh    = 2e7*Msun                  #g
R_bondi = 2*G*M_bh/(v_bh**2)        #cm
epsilon = 0.1*R_bondi               #cm

#Calculate Virial Temp at 1kpc
r_vir   = kpc        #cm          
T_vir   = G*M_bh*mp/(r_vir*kB)  #K

#Calculate Density if pressure equilibrium with CGM at 1kpc
rho_cgm = 1e-3*mp       #g/cm^3
T_cgm   = 1e6           #K
cs_cgm  = np.sqrt(kB*T_cgm/(mu*mp))
rho_vir = rho_cgm*T_cgm/T_vir

#Calculate Polytropic Constant
K = (kB*T_vir/(mu*mp))*(rho_vir**(1-gamma))
# print(f"Polytropic Constant = {K:.2e}")


#Code Units
L_code  = 1.0*kpc        #cm
v_code  = 1e6            #cm/s
rho_code= 1.0*mp         #g/cm^3

M_code  = rho_code*L_code**3
t_code  = L_code/v_code
P_code  = rho_code*v_code**2
T_code  = (mu*mp/kB)*(v_code**2)



print("="*10,"Code Units","="*90)
print(f"Code Length      = {L_code:.2e} cm \t = {L_code/kpc:.2e} kpc")
print(f"Code Mass        = {M_code:.2e} g \t = {M_code/Msun:.2e} Msun")
print(f"Code Time        = {t_code:.2e} s \t = {t_code/(1e6*s_yr):.2e} Myr")
print(f"Code Density     = {rho_code:.2e} g/cm^3 \t = {rho_code/mp:.2e} mp/cm^3")
print(f"Code Velocity    = {v_code:.2e} cm/s \t = {v_code/1e5:.2e} km/s")
print(f"Code Temperature = {T_code:.2e} K")
print(f"Code Pressure    = {P_code:.2e} dyne/cm^2")
# print(f"Code Pressure1    = {P_code1:.2e} dyne/cm^2")
#Calculate constants in code units
G_code  = G*(rho_code*t_code**2)
K_code  = K*((rho_code**gamma)/P_code)


print("="*10,"Constants in Code Units","="*77)
print(f"G_code = {G_code:.2e}")
print(f"K_code = {K_code:.2e}")

#Calculate Parameters in code units
print("="*10,"Parameters in Code Units","="*76)
print(f"Virial Radius   = {r_vir:.2e} cm \t\t = {r_vir/L_code:.2e} L_code")
print(f"Virial Density  = {rho_vir:.2e} g/cm^3 \t = {rho_vir/rho_code:.2e} rho_code")
print(f"Virial Temp     = {T_vir:.2e} K \t\t = {T_vir/T_code:.2e} T_code")
print(f"Black Hole Mass = {M_bh:.2e} g \t\t = {M_bh/M_code:.2e} M_code")
print(f"Black Hole Vel  = {v_bh:.2e} cm/s \t = {v_bh/v_code:.2e} v_code")
print(f"Bondi Radius    = {R_bondi:.2e} cm \t\t = {R_bondi/L_code:.2e} L_code")
print(f"Epsilon         = {epsilon:.2e} cm \t\t = {epsilon/L_code:.2e} L_code")
print(f"CGM Density     = {rho_cgm:.2e} g/cm^3 \t = {rho_cgm/rho_code:.2e} rho_code")
print(f"CGM Temp        = {T_cgm:.2e} K \t\t = {T_cgm/T_code:.2e} T_code")
print(f"CGM Sound Speed = {cs_cgm:.2e} cm/s \t = {cs_cgm/v_code:.2e} v_code")
print("="*112)

P_cgm   = rho_cgm*kB*T_cgm/(mu*mp)
P_vir   = K_code*((rho_vir/rho_code)**gamma)*P_code
print(f"CGM Pressure    = {P_cgm:.2e} dyne/cm^2 \t = {P_cgm/P_code:.2e} P_code")
print(f"Virial Pressure = {P_vir:.2e} dyne/cm^2 \t = {P_vir/P_code:.2e} P_code")