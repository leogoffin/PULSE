import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *


M = 0.75
v = M*speed_of_sound(0)
mdot = 165 * lb_to_kg #kg/s
pi_c = 15
eta_is_c = 0.88
delta_hf = 17800 * BTUlb_to_kJkg*1000
TIT = R_to_K(2500)
eta_cc = 0.98
pi_cc = 0.95
eta_is_t = 0.915
RR = 0.92
eta_shaft = 0.995


P0 = Air(v)
P0.set_isa(0)
P0.set_total()

P1 = P0.RAM(RR) #RAM
P2,Pc,Cp12 = compr.with_eta(P1,mdot,pi_c,eta_is_c) #compressor
P3,mdotf,f = cc.with_Tout(P2,P0.T0,mdot,delta_hf,TIT,eta_cc = eta_cc, pi_cc = pi_cc) # combustion chamber
P4,Cp34,T4is = turb.with_P(Pc,mdot,f,P3,TIT,eta_is = eta_is_t) # turbine

P5,chocked,NPR,A = adiabatic_convergent_nozzle(P4,mdot,f,P0.p) #nozzle
Tp,Ts = turbofan_thrust(mdot,0,f,P0,P5,AP= A) #thrust
SFC = mdotf/Tp

print(f"Dry thrust conv : {Tp:.0f} N")
print(f"Dry SFC : {SFC*kgN_to_kghDaN:.2f} kg/(hDaN)")
print(f"Nozzle thoughflow section : {A:.2f} m²")

P5b,Ab = adiabatic_convdiv_nozzle(P4,mdot,f,P0.p)
Tpb,Ts = turbofan_thrust(mdot,0,f,P0,P5b,AP= Ab) #thrust
print(f"Dry thrust conv+div : {Tpb:.0f} N")
P5b.print_state()