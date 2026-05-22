import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

A_ratio = 100
MHe = 4
R = get_R(MHe)

Pr = Air(R = R)
Pr.T0 = 310 #K
Pr.gamma = 5/3

Pa = Air()
Pa.p = 0
gamma = 5/3

Me = M_nozzle2(Pr,A_ratio,Mmax = 20)
Pe,p_ratio = exhaust_cds(Pr,Me)
c_T = get_cT(p_ratio,A_ratio,Pe.gamma,Pe.M)
print(f"Thrust Coeff : {c_T:.3f}")
c_star = get_c_star(Pr)
print(f"C_star : {c_star:.0f} m/s")
Isp = get_Isp_with_cT_star(c_T,c_star)
print(f"Isp : {Isp:.0f} s")

Pr.p0 = 10*1e5 #Pa
T = 0.2 #N
Pe = exhaust_cds(Pr,Me)
A_t = At_from_cT(T,c_T,Pr.p0)
A_e = A_t * A_ratio

mdot = mdot_from_Isp(T,Isp)
print(f"Mass flow rate : {mdot:.3e} kg")
print(f"Throat area : {A_t:.3e} mm²")

tot_impulse = 5000 #N*s
resm = 0.05
dt = tot_impulse/T
mT = mdot*dt
print(f"Useful mass : {mT:.3f} kg") 
m_init = mT/(1-resm)
print(f"Initial mass : {m_init:.3f} kg") 
pfinal = Pr.p0
pinit = pfinal/resm
#V is constant 
V = V_perfect_gas(m_init*resm,Pr)
print(f"Tank Volume : {V:.3e} m³")