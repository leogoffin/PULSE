import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

OPR = 50
RR = 0.995

eta_f_p = 0.95
eta_g = 0.995

pi_LPC = 2.15
eta_c_p = 0.94

eta_t_is = 0.92

eta_cc = 0.999
lam_cc = 0.03
pi_cc = 1 - lam_cc
delta_hf = 43.1e6

mdotf = 2.4
mdotp = 80

Pref = Air()
Pref.set_isa(0)

alt = 10000
M = 0.8
v = speed_of_sound(alt)*M
P0 = Air(v)
P0.set_isa(alt)

P1 = Air(v)
P1.set_isa(alt)


def ex1(BPR,pi_f,verbose = False):
    """Compute ex1."""
    mdot = mdotp*(BPR+1)
    mdots = mdotp * BPR

    # RAM
    P2,pi_r = P1.RAM(RR,True)
    #Fan
    P3,Pf,Cp23 = compr.with_eta_poly(P2,mdot,pi_f,eta_f_p)

    # Booster (LPC) 
    P4,PLPC,Cp34 = compr.with_eta_poly(P3,mdotp,pi_LPC,eta_c_p)

    # HPC
    pi_HPC = OPR/pi_f/pi_LPC
    P5,PHPC,Cp45 = compr.with_eta_poly(P4,mdotp,pi_HPC,eta_c_p)

    #combustion chamber 
    fcc = mdotf/mdotp
    P6 = cc.with_f(P5,Pref.T,mdotp,delta_hf,fcc,eta_cc = eta_cc, pi_cc = pi_cc)

    # HPT
    PHPT = PHPC
    P7,cp67,T7s = turb.with_P(PHPT,mdotp,fcc,P6,eta_is = eta_t_is)

    # LPT
    PLPT = Pf/eta_g + PLPC
    P8,cp78,T8s = turb.with_P(PLPT,mdotp,fcc,P7,eta_is = eta_t_is)

    # Primary exhaust
    P9,chokedp,NPRp,Ap = adiabatic_convergent_nozzle(P8,mdotp,fcc,P0.p)
    Tp = simple_thrust(mdotp,fcc,P0,P9,Ap)
    effv_p = effective_velocity(P9,mdotp*(1+fcc),P0.p,Ap)

    # Secondary exhaust
    P9s,chokeds,NPRs,As = adiabatic_convergent_nozzle(P3,mdots,0,P0.p)
    Ts = simple_thrust(mdots,0,P0,P9s,As)
    effv_s = effective_velocity(P9s,mdots,P0.p,As)

    T = Tp + Ts
    SFC = mdotf/T
    eta_p = propulsive_efficiency_turbofan(T,P0.v,mdotp,mdotf,mdots,effv_p,effv_s,mdot)
    eta_T = thermal_efficiency_turbofan(mdotp,mdotf,mdots,effv_p,effv_s,mdot,P0.v,delta_hf)

    if not verbose :
        return T

    print("="*50)
    print(f"Total Thrust                    : {T:.0f} N \n")
    print(f"Thrust (primary)                : {Tp:.0f} N ({Tp/T*100:.0f} %)")
    print(f"Choked nozzle ? (primary)       : {chokedp}")
    print(f"Thrust (secondary)              : {Ts:.0f} N ({Ts/T*100:.0f} %)")
    print(f"Choked nozzle ? (secondary)     : {chokeds}\n")


    print(f"SFC                    : {SFC*kgN_to_kghDaN:.2f} kg/DaN.h")
    print(f"Thermal efficiency     : {eta_T:.3f}")
    print(f"Propulsive efficiency  : {eta_p:.3f}")
    print(f"Overall efficiency     : {eta_p*eta_T:.3f}")
    print("="*50)
    return T


"""
1. BPR = 10 
"""
BPR = 10
pi_f = 1.5
ex1(BPR,pi_f,True)

"""
2. BPR = 12.5 
"""
BPR = 12.5
pi_f = 1.5

ex1(BPR,pi_f,True)

"""
3. BPR = 12.5 T = 140 kN
"""
BPR = 12.5
Thrust = 140e3 #N

def error(pi_f):
    """Return the error for the given input."""
    return ex1(BPR,pi_f,False) - Thrust

pi_f3,i = bisection(error,1.1,1.5)
ex1(BPR,pi_f3,True)