import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

a_ratio = 45.1

Pr = Air()
Pr.p0 = 110e5
mdot = 250
MR = 5.9
MH2 = 2e-3
MO2 = 32e-3
MH2O = 18e-3
gamma = 1.33

Me = 13.8
Re = get_R(Me)
T0 = C_to_K(3300)

#1 ground
Pa = Air(0)
Pa.set_isa(0)

P1 = Air(R = Re)
P1.p0 = 110e5
P1.gamma = 1.33
P1.T0 = C_to_K(3300)-273.15


c_star,A_t,A_e,Me,Pe = space_exhaust(P1,mdot,a_ratio)

print(f"Throat area           : {A_t:.4f} m²")
print(f"Exit area             : {A_e:.4f} m²")
print(f"Exit temperature (T*) : {Pe.T:.0f} K")
print(f"Exit Mach (Me*)       : {Me:.3f}")
print(f"Exit pressure (pe*)   : {Pe.p:.3f} Pa")

pa_vac, T_vac, Isp_vac, c_T_vac = vacuum_operation(Pe,P1.p0,mdot,A_e,A_t,doprint = True)

pa_full, T_full, Isp_full, c_T_full = fully_expanded(Pe, P1.p0, mdot, A_e, A_t, doprint=True)

pa_adapt, T_adapt, Isp_adapt, c_T_adapt = adapted_operation(Pe, P1.p0, mdot, A_t, doprint=True)

Pe_ground, T_ground, Isp_ground, c_T_ground = ground_operation(Pe,P1.p0, mdot,A_e, A_t, doprint=True)




