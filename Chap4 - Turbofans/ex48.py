import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

alt = 35000* ft_to_m
M = 0.8
RR = 0.97
BPR = 5.8
OPR = 30

pi_f = 1.7 
eta_p_fan = 0.9 #polytropic efficiency

Asec = 1.87 #m²

eta_p_c = 0.89 #polytropic
pi_c = OPR/pi_f

eta_p_t = 0.9

pi_cc = 1 - 0.05
eta_cc = 0.99
delta_hf = 43.1e6
TIT = C_to_K(1130)

v = M*speed_of_sound(alt)
P0 = Air(v)
P0.set_isa(alt)
P0.set_total()

P1 = P0.RAM(RR)
P2,Pf,Cp12 = compr.with_eta_poly(P1,1,pi_f,eta_p_fan)
P8b,chokedb,NPRb,Ab = adiabatic_convergent_nozzle(P2,1,0,P0.p)

mdots = Asec/Ab
mdotp = mdots/BPR
mdot = mdots + mdotp

print(f"Primary mass flow rate : {mdotp:.2f} kg/s")
print(f"Secondary mass flow rate : {mdots:.2f} kg/s")

P2,Pf,Cp12 = compr.with_eta_poly(P1,mdot,pi_f,eta_p_fan)
P6b,chokedb,NPRb,Ab = adiabatic_convergent_nozzle(P2,mdots,0,P0.p)
P3,Pc,Cp23 = compr.with_eta_poly(P2,mdotp,pi_c,eta_p_c)
P4,mdotf,f = cc.with_Tout(P3,P0.T0,mdotp,delta_hf,TIT,eta_cc = eta_cc,pi_cc = pi_cc)
Pt = Pc + Pf
P5,Cp45,T5s = turb.with_P_poly(Pt,mdotp,f,P4,eta_p = eta_p_t)

print(f"Fan power : {Pf/1e6:.2f} MW")
print(f"Compressor power : {Pc/1e6:.2f} MW")
print(f"Turbine power : {Pt/1e6:.2f} MW")

P6,choked,NPR,A = adiabatic_convergent_nozzle(P5,mdotp,f,P0.p)
Tp,Ts = turbofan_thrust(mdotp,mdots,f,P0,P6,P6b,AP = A,AS = Ab)
T = Tp + Ts
print(f"Thrust : {T/1e3:.2f} kN")
SFC = mdotf/T
print(f"SFC : {SFC*10*3600:.3f} kg/DaNh")
