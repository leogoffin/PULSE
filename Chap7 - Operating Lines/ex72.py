import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *


BPR = 20
Pref = Air()
Pref.set_isa(0)
Pref.set_total()

M = 0.85
alt = 10000* ft_to_m

r_in = 55/100/2
r_ext = 80/100/2

A_nozzle = np.pi*(r_ext**2-r_in**2)


m_dot_star_values = [74.04,    75.55,   78.27,   79.59,   81.10, 81.86]
pi_fs_values  = [1.3054,  1.3040,  1.2932,  1.2837,  1.2672, 1.2552]
eta_fs_values = [0.9016, 0.91377, 0.92997, 0.93454, 0.93422, 0.92911]

v = M* speed_of_sound(alt)
P1 = Air(v)
P1.set_isa(alt)
P1.set_total()

m_dot_real, m_dot_s_conv, pi_fs_conv, eta_fs_conv, P2s = compr.IterateMatchingFan_and_SecondaryNozzle(P1,Pref,pi_fs_values,eta_fs_values,m_dot_star_values,BPR,A_nozzle)
P2,Pf,Cp12 = compr.with_eta(P1,m_dot_real,pi_fs_conv,eta_fs_conv)
P3s,choked,NPR,A = adiabatic_convergent_nozzle(P2,m_dot_real,0,P1.p)
Ts = simple_thrust(m_dot_real,0,P1,P3s,A_nozzle)
print(f"Mass flow rate : {m_dot_real:.2f} kg/s")
print(f"Corrected mass flow rate : {corrected_massflow(m_dot_real,P1,Pref):.2f} kg/s")
print(f"Thrust {Ts:.0f} N")
print(f"Fan power : {Pf/1e6:.3f} MW")
print(f"Nozzle is choked. - {choked}")