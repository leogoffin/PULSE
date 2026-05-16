import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

pi_c = 35
m_dot_e = 159.21
T0_4 = 658.15
Pe = 51.2e6
HR = 2.40611
delta_hf = 47110e3 #J/kg
eta_c_is = 0.88
eta_cc = 1

Qc = HR*Pe
m_dot_f = Qc/eta_cc/delta_hf
m_dot = m_dot_e - m_dot_f
f = m_dot_f/m_dot

P1 = Air(0)     # Ambient air (static, ground cond)
P1.set_isa(0)   #isa conditions (ground)
P1.set_total()  #Gives all total(stagnation) values

P2, Pc, Cp12 = compr.with_eta(P1, m_dot, pi_c, eta_c_is) # compressor
P3 = cc.with_f(P2, Qc, P1.T0, m_dot, m_dot_e, f) # combustion chamber
P4 = Air()
P4.T0 = T0_4
Cp34 = findCp(av(P3.T0, P4.T0), f)
Pt = P_turb(m_dot_e, Cp34, f, P3, P4) # turbine

Pfa = Pt - Pc - Pe
print(f"Compressor Power : {Pc/1e6:.1f} MW")
print(f"Turbine Power : {Pt/1e6:.1f} MW")
print(f"Power Output : {Pe/1e6:.1f} MW")
print(f"Power losses : {Pfa/1e6:.1f} MW")

