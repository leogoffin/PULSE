import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

pi_LPC_s = 1.3
ns_LPC = 3
pi_LPC = pi_LPC_s**ns_LPC

pi_HPC_s = 1.5
ns_HPC = 5
pi_HPC = pi_HPC_s**ns_HPC

eta_c_p = 0.91

eta_t_is = 0.89
TOT = C_to_K(1250)
eta_m = 0.98

intercooler_DT = 45 #K 
#constant TIT
delta_hf = 42.6e6 #J/kg

mdot = 1 #we look at specific

Pa = Air(0)
Pa.set_isa(0)

P0 = Pa
P1 = P0

# LPC
P2,PLPC,Cp12 = compr.with_eta_poly(P1,mdot,pi_LPC,eta_c_p)

# Intercooler
#P2b = P2
#P2b.T0 -= intercooler_DT

# HPC
P3,PHPC,Cp23 = compr.with_eta_poly(P2,mdot,pi_HPC,eta_c_p)

Pt = (PLPC + PHPC) / eta_m

def error(fcc):
    P4 = cc.with_f(P3,P0.T,mdot,delta_hf,fcc)
    P5,cp45,T5s = turb.with_P(Pt,mdot,fcc,Pin = P4,eta_is = eta_t_is)
    return P5.T0 - TOT

fcc,i = bisection(error,0,0.05)

P4 = cc.with_f(P3,P0.T,mdot,delta_hf,fcc)
TIT = P4.T0
P5,cp45,T5s = turb.with_P(Pt,mdot,fcc,P4,eta_is = eta_t_is)
P6,choked,NPR,A = adiabatic_convergent_nozzle(P5,mdot,fcc,Pa.p)
T = simple_thrust(mdot,fcc,P0,P6,A)
SFC = fcc/T
print(f"TIT                                    : {TIT:.0f} K")
print(f"Thrust (no_intercooler)                : {T:.0f} N")
print(f"Exhaust Area (no-intercooler)          : {A:.3f} m²")
print(f"SFC (no-intercooler)                   : {SFC*kgN_to_kghDaN:.2f} kg/DaN.h")
print(f"FAR (no-intercooler)                   : {fcc:.3f}\n")


# LPC
P2b,PLPC,Cp12 = compr.with_eta_poly(P1,mdot,pi_LPC,eta_c_p)

# Intercooler
P2b.T0 -= intercooler_DT

# HPC
P3b,PHPC,Cp23 = compr.with_eta_poly(P2b,mdot,pi_HPC,eta_c_p)

Pt = (PLPC + PHPC) / eta_m

P4b,mdotfb,fccb = cc.with_Tout(P3b,P0.T,mdot,delta_hf,TIT)
P5b,cp45,T5s = turb.with_P(Pt,mdot,fccb,P4b,eta_is = eta_t_is)
P6b,choked,NPR,A = adiabatic_convergent_nozzle(P5b,mdot,fccb,Pa.p)
T = simple_thrust(mdot,fccb,P0,P6b,A)
SFC = fccb/T

print(f"Thrust (intercooler)                : {T:.0f} N")
print(f"Exhaust Area (intercooler)          : {A:.3f} m²")
print(f"SFC (intercooler)                   : {SFC*kgN_to_kghDaN:.2f} kg/DaN.h")
print(f"FAR (intercooler)                   : {fccb:.3f}")

P1.sprint(1)
P2.sprint(2)
P3.sprint(3)
P4.sprint(4)
P5.sprint(5)
P6.sprint(6)

P1.sprint(1,"intercooler")
P2b.sprint(2,"intercooler")
P3b.sprint(3,"intercooler")
P4b.sprint(4,"intercooler")
P5b.sprint(5,"intercooler")
P6b.sprint(6,"intercooler")