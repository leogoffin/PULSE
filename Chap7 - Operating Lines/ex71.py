import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *



M0 = 0.82
altitude = 35000 * ft_to_m #m
a = speed_of_sound(altitude)
v = M0*a

pi_fp = 1.49
pi_fs = 1.6
eta_is_fs = 0.84 #from table
eta_is_fp = 0.83
m_dot_corr = 415    #kg/s

BPR = 5 
print(275*BPR/(BPR+1))
m_dot_corr_p = m_dot_corr /(BPR+1)
print(f"Corrected core mass flow rate : {m_dot_corr_p:.0f} kg/s")
#ref cds (ground)
Pref = Air()
Pref.set_isa(0)
Pref.set_total()  

P0 = Air(v)     # Ambient air (static, ground cond)
P0.set_isa(altitude)   #isa conditions (ground)
P0.set_total()  #Gives all total(stagnation) values

P1 = P0.RAM()

mdot = true_massflow_from_corrected(m_dot_corr,P1,Pref)
print(f"Mass flow rate : {mdot:.2f} kg/s")

mdot_s = mdot * BPR / (BPR + 1)
mdot_p = mdot / (BPR + 1)
print(f"Mass flow rate through fan : {mdot_s:.2f} kg/s")
print(f"Mass flow rate through core : {mdot_p:.2f} kg/s")

P2s,Pfs,Cp12s = compr.with_eta(P1,mdot_s,pi_fs,eta_is_fs)

P5s,chockeds,NPRs,As = adiabatic_convergent_nozzle(P2s,mdot_s,0,P0.p)
print(f"Nozzle pressure ratio : {NPRs:.2f}")
print(f"Fan exit area : {As:.3f} m^2")

Ts = simple_thrust(mdot_s,0,P0,P5s,A = As)
print(f"Secondary thrust : {Ts:.0f} N")

P2,Pfp,Cp12 =  compr.with_eta(P1,mdot_p,pi_fp,eta_is_fp)

print(f"Primary fan power : {Pfp/1e6:.2f} MW")
print(f"Secondary fan power : {Pfs/1e6:.2f} MW")
print(f"Total fan power : {(Pfs+Pfp)/1e6:.2f} MW")


m_dot_star_values = [275, 300, 325, 339]
eta_fs_values = [0.746, 0.793, 0.839, 0.82]
pi_fs_values = [1.34, 1.37, 1.38, 1.25]

m_dot_real, m_dot_s_conv, pi_fs_conv, eta_fs_conv, P2 = compr.IterateMatchingFan_and_SecondaryNozzle(P1,Pref,pi_fs_values,eta_fs_values,m_dot_star_values,BPR,As)

print(f"{m_dot_real:.0f} kg/s")



