import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

M0 = 0.82
altitude = 35000 * ft_to_m #m
a = speed_of_sound(altitude)
v = M0*a

pi_fp = 1.49
pi_fs = 1.6
m_dot_star = 415    #kg/s
BPR = 5 


P0 = Air(v)     # Ambient air (static, ground cond)
P0.set_isa(altitude)   #isa conditions (ground)
P0.set_total()  #Gives all total(stagnation) values

P1 = P0.RAM()
P0.print_state()
p0_star = 101325 #Pa
T0_star = 288.15 #K
mdot = m_dot_star * P0.p0 / p0_star * np.sqrt(T0_star / P0.T0)

print(f"Mass flow rate : {mdot:.2f} kg/s")
NPR = pi_fp * P0.p0 / P0.p 
NPR_star = ((P0.gamma + 1 )/2)**(P0.gamma/(P0.gamma-1))
print(f"Nozzle pressure ratio : {NPR:.2f}")
print(f"Choking nozzle pressure ratio : {NPR_star:.2f}")

mdot_s = mdot * BPR / (BPR + 1)
mdot_p = mdot / (BPR + 1)
print(f"Mass flow rate through fan : {mdot_s:.2f} kg/s")
print(f"Mass flow rate through core : {mdot_p:.2f} kg/s")

P5 = Air(0)     # Ambient air (static, ground cond)
P5.p0 = P0.p0 * pi_fs
P5.T0 = P0.T0 * (1 + pi_fs**((P0.gamma-1)/P0.gamma) - 1) 
a0 = (P0.v/P0.M) * ((1 + (P0.gamma - 1)/2 * P0.M**2) ** 0.5)

crit_mass_fr = P0.rho0 * a0 * (2/(P0.gamma+1))**((P0.gamma+1)/(2*(P0.gamma-1)))
As = mdot_s / crit_mass_fr
print(f"Fan exit area : {As:.4f} m^2")
#T_s = mdot_s * (vs - v) + (ps - P0.p0) * As