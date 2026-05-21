import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

A_ratio = 45.1
mdot = 250 #kg/s
MR = 5.9
g = 9.81

MH2 = 1e-3*2 #kg/mole
MO2 = 16e-3*2 #kg/mole
MH2O = (MO2 + 2*MH2)/2

mH2 = mdot/(1+MR)
mO2 = mdot*MR/(1+MR) 

results = combustion_analysis(
    species_in=["H2", "O2"],nu_in=[2, 1],M_in=[MH2, MO2],quantities_in=[mH2, mO2],
    species_out=["H2O"],nu_out=[2],M_out=[MH2O],
)


P1 = Air(R = results["R"])
P1.p0 = 110e5
P1.gamma = 1.33
P1.T0 = C_to_K(3300)-273.15
c_star = charac_vel(P1)
print(c_star)
A_t = A_nozzle(P1,mdot,1)
A_e = A_t * A_ratio
print(f"Throat area : {A_t:.4f} m²")
Me = M_nozzle(P1,mdot,A_e)
Pe = exhaust_cds(P1,Me)
Pe.print_state()

T = space_thrust(mdot,Pe,A_e)
Isp = T/mdot/g
c_T = T/P1.p0/A_t
print(f"Thrust (actual)           : {T:.0f} N")
print(f"Isp (actual)              : {Isp:.1f} s")
print(f"Thrust coefficient (act.) : {c_T:.3f}")

pa_full, T_full, Isp_full, c_T_full = fully_expanded(Pe, P1.p0, mdot, A_e,A_t, doprint=True)
pa_adapt, T_adapt, Isp_adapt, c_T_adapt = adapted_operation(Pe, P1.p0, mdot, A_t, doprint=True)
Pe_ground, T_ground, Isp_ground, c_T_ground = ground_operation(P1, mdot,A_t, doprint=True)

