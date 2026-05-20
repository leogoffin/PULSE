import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

pi_c = 4
n = 21600           # RPM
mdot = 7.6          # kg/s
mdot_f = 0.1191     # kg/s
f = mdot_f/mdot
delta_hf = 42.7e6   # J/kg
TIT = C_to_K(800)   # C -> K
eta_is_t = 0.8

P1 = Air()
P1.set_isa(0)
P1.set_total()


P3 = Air()
P3.T0 = TIT
P2 = cc.rev_with_f(P3,delta_hf,P1.T0,mdot,f)
eta_is_c,Pc = compr.get_eta(P1,P2,mdot,pi_c)
print("\n")
print(f"Compressor isentropic efficiency : {eta_is_c:.3f}")
print(f"Compressor work : {Pc/1e6:.2f} MW")
P3.p0 = P2.p0
P4, Cp34,T4_is = turb.with_P(Pc,mdot,f,Pin = P3,eta_is = eta_is_t)
print(f"Turbine isentropic outlet temperature : {T4_is:.0f} K")

P5, choked,NPR,A = adiabatic_convergent_nozzle(P4,mdot,f,P1.p)
Tp,Ts = turbofan_thrust(mdot,0,f,P1,P5)
print(f"Thrust : {Tp:.0f} N")
SFC = mdot_f/Tp
print(f"SFC : {SFC*1000*3600:.3f} kg/kNh")
print(f"SFC : {SFC*10*3600:.3f} kg/DaNh")
eta_t = thermal_efficiency(P1.v,P5.v,mdot,f,delta_hf)
eta_p = propulsive_efficiency(P1.v,P5.v,mdot,f,Tp)
print(f"Thermal efficiency : {eta_t:.3f}")
print(f"Propulsive efficiency : {eta_p:.3f}")
print("\n \n \n")




M = 0.65
alt = 9000
v = M*speed_of_sound(alt)
P0 = Air(v)
P0.set_isa(alt)
P0.set_total()

P1_corr = P0.RAM()

mdot_corr = true_massflow_from_corrected(mdot,P1_corr,P1)
n_corr = true_rpm_from_corrected(n,P1_corr,P1)

print(f"Corrected mass flow : {mdot_corr:.3f} kg/s")
print(f"Corrected rpm : {n_corr:.0f} RPM")

P2_corr,Pc_corr,Cp12_corr = compr.corrected(P1,P2,P1_corr,mdot_corr)
print(f"Corrected Compressor power : {Pc_corr/1e3:.0f} kW")

P3_corr,mdot_f_corr,f_corr = cc.with_Tout(P2_corr,P0.T0,mdot_corr,delta_hf,TIT)
print(f"Corrected fuel mass flow rate : {mdot_f_corr:.3f}")

P4_corr,Cp34_corr,T4s_corr = turb.with_P(Pc_corr,mdot_corr,f_corr,P3_corr,eta_is = eta_is_t)
P5_corr, choked_corr, NPR_corr, A_corr = adiabatic_convergent_nozzle(P4_corr,mdot_corr,f,P0.p)

Tp_corr,Ts_corr = turbofan_thrust(mdot_corr,0,f_corr,P0,P5_corr,AP = A_corr) 
print(f"Corrected Thrust : {Tp_corr:.0f} N")
SFC_corr = mdot_f_corr/Tp_corr
print(f"Corrected SFC : {SFC_corr*1000*3600:.3f} kg/kNh")
print(f"Corrected SFC : {SFC_corr*10*3600:.3f} kg/DaNh")

eta_t_corr = thermal_efficiency(P0.v,P5_corr.v,mdot_corr,f_corr,delta_hf)
eta_p_corr = propulsive_efficiency(P0.v,P5_corr.v,mdot_corr,f_corr,Tp_corr)
print(f"Corrected Thermal efficiency : {eta_t_corr:.3f}")
print(f"Corrected Propulsive efficiency : {eta_p_corr:.3f}")