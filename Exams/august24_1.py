import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

BPR = 0.63
OPR = 25
eta_c_p = 0.95
eta_t_is = 0.98
eta_cc = 0.95
lam_cc = 0.03
pi_cc = 1 - lam_cc

eta_m = 0.95
delta_hf = 43e6
mdotfcc = 1.35
mdotfab = 2.4
mdot = 103
mdotp = mdot/(1+BPR)
mdots = BPR*mdotp

Pref = Air(0)
Pref.set_isa(0)

M = 1.2
alt = 10000
v = M*speed_of_sound(alt)
P0 = Air(v)
P0.set_isa(alt)



P1 = Air(v)
P1.set_isa(alt)

RR = P0.Ram_recovery()
print(f"RAM recovery factor : {RR:.3f}")
P2,pi_r = P1.RAM(RR)

def error(pi_f):

    pi_HPC = OPR/pi_f
    P3,Pf,cp23 = compr.with_eta_poly(P2,mdot,pi_f,eta_c_p)
    P4,PHPC,cp34 = compr.with_eta_poly(P3,mdotp,pi_HPC,eta_c_p) 
    fcc = mdotfcc/mdotp
    P5 = cc.with_f(P4,Pref.T,mdotp,delta_hf,fcc,eta_cc = eta_cc,pi_cc = pi_cc)

    PLPT = Pf/eta_m
    PHPT = PHPC/eta_m
    P6,cp56,T6s = turb.with_P(PHPT,mdotp,fcc,P5,eta_is =eta_t_is)
    P7,cp67,T7s = turb.with_P(PLPT,mdotp,fcc,P6,eta_is =eta_t_is)
    return P7.p0 - P3.p0

pi_f,i = bisection(error,1.1,5)

pi_HPC = OPR/pi_f

P3,Pf,cp23 = compr.with_eta_poly(P2,mdot,pi_f,eta_c_p)
P4,PHPC,cp34 = compr.with_eta_poly(P3,mdotp,pi_HPC,eta_c_p) 
fcc = mdotfcc/mdotp
P5 = cc.with_f(P4,Pref.T,mdotp,delta_hf,fcc,eta_cc = eta_cc,pi_cc = pi_cc)

PLPT = Pf/eta_m
PHPT = PHPC/eta_m
P6,cp56,T6s = turb.with_P(PHPT,mdotp,fcc,P5,eta_is =eta_t_is)
P7,cp67,T7s = turb.with_P(PLPT,mdotp,fcc,P6,eta_is =eta_t_is)

P8,f_8 = mixing(P7,P3,Pref.T,mdotp,mdots,fcc)

P9d = P8

P10d,choked,NPR,Ad = adiabatic_convergent_nozzle(P9d,mdot,f_8,P0.p)

Td = simple_thrust(mdot,f_8,P0,P10d,Ad)
SFCd = mdotfcc/Td
veff_d = effective_velocity(P10d,mdot,P0.p,Ad)
eta_td = thermal_efficiency(v,veff_d,mdot,f_8,delta_hf)
eta_pd = propulsive_efficiency(v,veff_d,mdot,f_8,Td)


print('-'*50)
print(f"Fan pressure ratio          : {pi_f:.2f}")
print(f"Nozzle choked ?             : {choked}")
print(f"Thrust (dry)                : {Td:.0f} N")
print(f"Throat area (dry)           : {Ad:.3f} m²\n")
print(f"SFC (dry)                   : {SFCd*kgN_to_kghDaN:.2f} kg/DaN.h")
print(f"Thermal efficiency (dry)    : {eta_td:.3f}")
print(f"Propulsive efficiency (dry) : {eta_pd:.3f}")
print(f"Overall efficiency (dry)    : {eta_pd*eta_td:.3f}")

fab = mdotfab/mdot
f_final = fab + f_8
P9w = cc.with_f(P8,Pref.T,mdot,delta_hf,fab,f_8,eta_cc = eta_cc,pi_cc = pi_cc)

P10w,choked,NPR,Aw = adiabatic_convergent_nozzle(P9w,mdot,f_final,P0.p)

Tw = simple_thrust(mdot,f_final,P0,P10w,Aw)
SFCw = mdotfcc/Tw
veff_w = effective_velocity(P10w,mdot,P0.p,Aw)
eta_tw = thermal_efficiency(v,veff_w,mdot,f_final,delta_hf)
eta_pw = propulsive_efficiency(v,veff_w,mdot,f_final,Tw)


print('-'*50)
print(f"Fan pressure ratio          : {pi_f:.2f}")
print(f"Nozzle choked ?             : {choked}")
print(f"Thrust (wet)                : {Tw:.0f} N")
print(f"Throat area (wet)           : {Aw:.3f} m²\n")
print(f"SFC (wet)                   : {SFCw*kgN_to_kghDaN:.2f} kg/DaN.h")
print(f"Thermal efficiency (wet)    : {eta_tw:.3f}")
print(f"Propulsive efficiency (wet) : {eta_pw:.3f}")
print(f"Overall efficiency (wet)    : {eta_pw*eta_tw:.3f}")

P1.sprint(1)
P2.sprint(2)
P3.sprint(3)
P4.sprint(4)
P5.sprint(5)
P6.sprint(6)
P6.sprint(6)
P7.sprint(7)
P8.sprint(8)
P9d.sprint(9)
P10d.sprint(10)