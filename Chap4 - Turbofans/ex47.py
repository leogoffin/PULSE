import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *


M = 0.8
alt = 350 *fl_to_m
TIT = 1460
OPR = 35.5
RR = 0.995
BPR = 5.2
Treq = 28558 #N

pi_f = 1.7
eta_is_f = 0.895

pi_boost = 1.8
eta_is_boost = 0.875

eta_is_HPC = 0.875
pi_HPC = OPR/pi_f/pi_boost

pi_cc = 1 - 0.06
eta_is_cc = 0.999
delta_hf = 42.9e6

eta_is_HPT = 0.91
eta_is_LPT = 0.92

eta_m = 0.995

mdotp = 1 #per unit mass
mdots = mdotp*BPR
mdot = mdotp*(1+BPR)


v = M*speed_of_sound(alt)
P0 = Air(v)
P0.set_isa(alt)
P0.set_total()

P1 = P0.RAM(RR)
P2,Pfan,Cp12 = compr.with_eta(P1,mdot,pi_f,eta_is_f)
P3,Pboost,Cp23 = compr.with_eta(P2,mdotp,pi_boost,eta_is_boost)
P4,PHPC,Cp34 = compr.with_eta(P3,mdotp,pi_HPC,eta_is_HPC)
P5,mdotf,f = cc.with_Tout(P4,P0.T0,mdotp,delta_hf,TIT,eta_cc = eta_is_cc,pi_cc = pi_cc)
P6,Cp56,T6s = turb.with_P(PHPC/eta_m,mdotp,f,Pin = P5,eta_is = eta_is_HPT)
P7,Cp67,T7s = turb.with_P((Pfan+Pboost)/eta_m,mdotp,f,Pin = P6,eta_is= eta_is_LPT)
P8,choked,NPR,A = adiabatic_convergent_nozzle(P7,mdotp,f,P0.p)

P8b,chokedb,NPRb,Ab = adiabatic_convergent_nozzle(P2,mdots,0,P0.p)
Tp,Ts = turbofan_thrust(mdotp,mdots,f,P0,P8,P8b,AP=A,AS=Ab)
T = Tp + Ts
mdotp_req = mdotp*Treq/T
print(f"Thrust per unit mass : {T:.0f} N/kg")
print(f"Thrust required : {Treq:.0f} N")
print(f"Required primary mass flow rate : {mdotp_req:.2f} kg/s")


"""
Mixed - iterative fan pressure ratio
"""
mdotp = mdotp_req
mdots = mdotp*BPR
mdot = mdotp*(1+BPR)
converged = False
errs = []
while not converged and len(errs)<100 :

    converged,var,errs = Pcontrol(P7.p0,P2.p0,errs,Kp = 0.3)
    pi_f *= (1+var)
    pi_HPC = OPR/pi_f/pi_boost
    if not converged :
        P2,Pfan,Cp12 = compr.with_eta(P1,mdot,pi_f,eta_is_f)
        P3,Pboost,Cp23 = compr.with_eta(P2,mdotp,pi_boost,eta_is_boost)
        P4,PHPC,Cp34 = compr.with_eta(P3,mdotp,pi_HPC,eta_is_HPC)
        P5,mdotf,f = cc.with_Tout(P4,P0.T0,mdotp,delta_hf,TIT,eta_cc = eta_is_cc,pi_cc = pi_cc)
        P6,Cp56,T6s = turb.with_P(PHPC,mdotp,f,Pin = P5,eta_is = eta_is_HPT)
        P7,Cp67,T7s = turb.with_P(Pfan+Pboost,mdotp,f,Pin = P6,eta_is= eta_is_LPT)

print(f"Perfect mixing fan pressure ratio : {pi_f:.3f}")

Pm,global_f = mixing(P7,P2,P0.T0,mdotp,mdots,f)
P8,chokedm,NPRm,Am = adiabatic_convergent_nozzle(Pm,mdot,global_f,P0.p)
Tm = simple_thrust(mdot,global_f,P0,P8,AP=Am)
SFC = mdotp*f/Tm
print(f)
print(f"Mixed exhaust Thurst : {Tm:.0f} N")
print(f"Mixed SFC : {SFC*10*3600:.3f} kg/DaNh")