import numpy as np

#Constants
G           = 6.67e-8       #cm^3/g/s^2
kB          = 1.38e-16      #erg/K
mp          = 1.67e-24      #g
kpc         = 3.086e21      #cm
pc          = 3.086e18      #cm
Msun        = 1.989e33      #g
mu          = 1.0
s_yr        = 3.154e7       #s
s_Myr       = 1e6*s_yr
# K           = 2.91e28       #g^(-2/3)cm^4s^(-2)
gamma       = 5.0/3.0       #Polytropic index
gm1         = gamma-1
kms         = 1e5


rho_cgm0    = 1e1*mp                   #g/cm^3
rho_cgm1    = 1e-1*mp                   #g/cm^3

vw_0        = 100*kms                    #cm/s
vw_1        = 1000*kms                   #cm/s

t1          = 100*s_Myr                  #s
tau         = 1*s_Myr                   #s

alpha_cgm   = np.log(rho_cgm1/rho_cgm0)/np.log(1+(t1/tau))
alpha_vw    = np.log(vw_1/vw_0)/np.log(1+(t1/tau))
print(f"alpha_cgm = {alpha_cgm:.2f}") 
print(f"alpha_vw  = {alpha_vw:.2f}")
