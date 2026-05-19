import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *
import matplotlib.pyplot as plt

mdot = 164 * lbs_to_kgs # kg
pi_s = 1.3
eta_p = 0.9
stages = 10 

P = Air()
P.p0 = 57452 #Pa
P.T0 = 266.42 #K

total_pi = pi_s**stages
P1,Pc,Cp01,eta_is = compr.with_eta_poly_per_stage(P,mdot,total_pi,eta_p) 
eta_is_stage = eta_is_from_eta_poly(pi_s,eta_p,get_gamma(Cp01))

print(f"Total isentropic efficiency : {eta_is:.3f}")
print(f"Isentropic efficiency per stage : {eta_is_stage:.3f}")
print(f"Total compressor power : {Pc/1e6:.2f} MW")





Pc = 0
Ps = []
p = []
p.append(P.p0)
Ps.append(P)

for i in range(1,stages+1):
    P,dPc,Cp = compr.with_eta_poly(Ps[-1],mdot,pi_s,eta_p)
    Ps.append(P)
    p.append(P.p0)
    Pc+=dPc
print(f"Total compressor power (chain) : {Pc/1e6:.2f} MW")
