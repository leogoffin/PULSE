import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

OPR = 24.5
BPR = 0.3
mdot = 65 #kg/s
mdotp = mdot/(1+BPR)
mdots = mdotp*BPR

TIT = 1850 #K
TiT = 1560 #K
SFC_dry = 0.8 / 10 / 3600 #kg/sN
T_dry = 11250*lbf_to_N #N
mdotf_dry = SFC_dry*T_dry
f_dry = mdotf_dry/mdotp
f_m = mdotf_dry/mdot

SFC_wet = 1.7 / 10 / 3600 #kg/sN
T_wet = 16860*lbf_to_N #N
mdotf_wet = SFC_wet*T_wet
f_final = mdotf_wet/mdot
pi_s = 0.99
eta_p_c = 0.875  #polytropic, guess 12.5% loss

pi_cc = 1- 0.03
eta_is_cc = 1
pi_ab = 1 - 0.02
delta_hf = 42.8e6
eta_is_t = 0.9

pi_m = 1 - 0.03

P1 = Air()
P1.set_isa(0)
P1.set_total()

P4 =  Air()
P4.T0 = TiT
P4.p0 = P1.p0*OPR*pi_cc
P3 = cc.rev_with_f(P4,delta_hf,P1.T0,mdotp,f_dry)
eta_p_c = compr.get_eta_poly(P1,P3,OPR)
eta_is_c = eta_is_from_eta_poly(OPR,eta_p_c)

pi_f1 = 3
converged = False
errs = []
while not converged and len(errs)<100 :
     
    P2,Pf,Cp12 = compr.with_eta_poly(P1,mdot,pi_f1,eta_p = eta_p_c) 
    pi_c = OPR/pi_f1
    P3,Pc,Cp23 = compr.with_eta_poly(P2,mdotp,pi_c,eta_p = eta_p_c)
    P5,Cp45,T5s = turb.with_P(Pf+Pc,mdotp,f_dry,Pin = P4,eta_is = eta_is_t)
    converged,var,errs = Pcontrol(P2.p0*pi_s,P5.p0, last_errors=errs, Kp=0.7,tol = 1e-4)
    pi_f1*= (1-var)

P6,mixed_f = mixing(P5,P2,P1.T0,mdotp,mdots,f_dry,pi_m = pi_m) 

P7,choked,NPR,A = adiabatic_convergent_nozzle(P6,mdot,mixed_f,P1.p)
T = simple_thrust(mdot,mixed_f,P1,P7,A = A)
print(f"Dry Thrust {T:.0f} N")

P7b = cc.with_f(P6,P1.T0,mdot,delta_hf,f_final-mixed_f,mixed_f,pi_ab)
P8b,chokedb,NPRb,Ab = adiabatic_convergent_nozzle(P7b,mdot,f_final,P1.p)
Tb = simple_thrust(mdot,f_final,P1,P8b,A = Ab)
print(f"Wet Thrust {Tb:.0f} N")
#P7.print_state()