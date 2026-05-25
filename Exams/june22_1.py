import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

OPR = 40
BPR = 9
pi_f = 1.9
eta_f_p = 0.95
n_stage_boost = 3
pi_stage_boost = 1.4
eta_c_p = 0.94
eta_t_p = 0.96
TIT = C_to_K(1550)
eta_intercooler = 0.8
mdot = 1
mdotp = mdot/(1+BPR)
mdots = mdotp * BPR
delta_hf = 43.8e6

M = 0.8
alt = 10000
v = M*speed_of_sound(alt)

Pref = Air(0)
Pref.set_isa(0)

Pa = Air(v)
Pa.set_isa(alt)

P0 = Air(v)
P0.set_isa(alt)
P0.set_total()

P1,pi_r = P0.RAM()

P2,Pf,cp12 = compr.with_eta_poly(P1,mdot,pi_f,eta_f_p)
pi_LPC = pi_stage_boost**n_stage_boost
P3,PLPC,cp23 = compr.with_eta_poly(P2,mdotp,pi_LPC,eta_c_p)
pi_HPC = OPR/pi_f/pi_LPC
P5,PHPC,cp45 = compr.with_eta_poly(P3,mdotp,pi_HPC,eta_c_p)
P6,mdotf,fcc = cc.with_Tout(P5,Pref.T,mdotp,delta_hf,TIT)
P7,cp67,T7s,pi_LPT = turb.with_eta_poly(PHPC,mdotp,fcc,P6,eta_p = 0.96)
P8,cp78,T8s,pi_HPT = turb.with_eta_poly(PLPC+Pf,mdotp,fcc,P7,eta_p = 0.96)
#P7,cp67,pi_HPT = turb.with_P_poly(PHPC,mdotp,fcc,P6,eta_p = eta_t_p)
#P8,cp78,pi_LPT = turb.with_P_poly(PLPC+Pf,mdotp,fcc,P7,eta_p = eta_t_p)

P9p,NPR,Ap = adiabatic_convdiv_nozzle(P8,mdotp,fcc,P0.p)
P9s,NPR,As = adiabatic_convdiv_nozzle(P2,mdots,0,P0.p)

Tp = simple_thrust(mdotp,fcc,P0,P9p,Ap)
Ts = simple_thrust(mdots,0,P0,P9s,As)
T = Tp + Ts
SFC = mdotf/T

print("="*50)
print("No-Intercooler")
"""
print(f"FAR : {fcc:.3f}")
print(f"Fan power                       : {Pf:.0f} W")
print(f"Booster power                   : {PLPC:.0f} W")
print(f"HPC power                       : {PHPC:.0f} W")
"""
print(f"Total Thrust                    : {T:.0f} N \n")
print(f"Thrust (primary)                : {Tp:.0f} N ({Tp/T*100:.0f} %)")
print(f"Thrust (secondary)              : {Ts:.0f} N ({Ts/T*100:.0f} %)")
print(f"SFC                             : {SFC*kgN_to_kghDaN:.2f} kg/DaN.h")
print(f"Specific Thrust                 : {T/mdot:.0f} N \n")





P4p,P4s = intercooler(P3,P2,eta_intercooler,mdotp,mdots)
pi_HPC = OPR/pi_f/pi_LPC
P5,PHPC,cp45 = compr.with_eta_poly(P4p,mdotp,pi_HPC,eta_c_p)
P6 = cc.with_f(P5,Pref.T,mdotp,delta_hf,fcc)
P7,cp67,T7s,pi_LPT = turb.with_eta_poly(PHPC,mdotp,fcc,P6,eta_p = 0.96)
P8,cp78,T8s,pi_HPT = turb.with_eta_poly(PLPC+Pf,mdotp,fcc,P7,eta_p = 0.96)
#P7,cp67,pi_LPT = turb.with_P_poly(PHPC,mdotp,fcc,P6,eta_p = eta_t_p)
#P8,cp78,pi_HPT = turb.with_P_poly(PLPC+Pf,mdotp,fcc,P7,eta_p = eta_t_p)
P9p,NPR,Ap = adiabatic_convdiv_nozzle(P8,mdotp,fcc,P0.p)
P9s,NPR,As = adiabatic_convdiv_nozzle(P4s,mdots,0,P0.p)

Tp = simple_thrust(mdotp,fcc,P0,P9p,Ap)
Ts = simple_thrust(mdots,0,P0,P9s,As)
T = Tp + Ts
SFC = mdotf/T
Tsp = T/mdot

print("="*50)
print("Intercooler")
print(f"Total Thrust                    : {T:.0f} N \n")
print(f"Thrust (primary)                : {Tp:.0f} N ({Tp/T*100:.0f} %)")
print(f"Thrust (secondary)              : {Ts:.0f} N ({Ts/T*100:.0f} %)")
print(f"Fuel consumption                : {mdotf:.3f} kg/s")

print(f"SFC                             : {SFC*kgN_to_kghDaN:.2f} kg/DaN.h")

