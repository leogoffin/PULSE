import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *


M = 0.3
alt = 121.92
m_dot_p = 90
alpha = 6
m_dot_s = alpha*m_dot_p
m_dot = m_dot_p+m_dot_s

pi_f = 2.2
PI = 36
pi_C = PI/pi_f
T4_tot = 1680
pi_cc = 0.95
eta_cc = 0.99
Delta_hf = 41.4*1e6 #(J/kg)
eta_f_s = 0.92
eta_c_s = 0.9
eta_hpt_s = 0.91
eta_lpt_s = 0.91
RR = 0.97

v = M* speed_of_sound(alt)
P0 = Air(v)
P0.set_isa(alt)
P0.set_total()
T0_ref = P0.T0

P1 = P0.RAM(RR) # RAM effect, adaiabatic inlet
P2,P_FAN,Cp12 = compr.with_eta(P1,m_dot,pi_f,eta_f_s) # Homogeneous Fan == compressor
P3,P_C,Cp23 = compr.with_eta(P2,m_dot_p,pi_C,eta_c_s) # Compressor
P4,m_dot_f, f = cc.with_Tout(P3,T0_ref,m_dot_p,Delta_hf, T4_tot, eta_cc = eta_cc, pi_cc = pi_cc) # Combustion chamber : know T_out, f unknown

P_T = P_FAN + P_C
P5,Cp45,T5s = turb.with_P(P_T,m_dot_p,f,Pin = P4,eta_is = eta_hpt_s) # Low and high pressure turbine together eta_hpt_s = eta_lpt_s => single turb
P8,chocked,NPR,A = adiabatic_convergent_nozzle(P5,m_dot_p+m_dot_f,f,P0.p) # Adiabatic convergent nozzle (core)
P8_b,chocked_b,NPR_b,A_b = adiabatic_convergent_nozzle(P2,m_dot_s,0,P0.p) # Adiabatic convergent nozzle (bypass)

Tp,Ts = turbofan_thrust(m_dot_p,m_dot_s,f,P0,P8,P8_b,A,A_b)
T = Tp + Ts
SFC = m_dot_f/T

print(f"Fan power : {(P_FAN/1e6):.2f} MW")
print(f"Compressor work : {P_C/1e6:.2f} MW\n")
print(f"Fuel to air ratio : {f:.4f}")
print(f"Turbine Inlet Temperature :  {P4.T0:.0f} K")
print(f"Turbine work : {P_T/1e6:.2f} MW")
print(f"Exhaust velocity : {P8.v:.0f}  m/s\n")
print(f"Primary flow Thrust : {Tp/1000:.3f} kN")
print(f"Secondary flow Thrust : {Ts/1000:.3f} kN")
print(f"Total Thrust : {T/1000:.3f} kN\n")
print(f"SFC : {SFC*1000*3600:.3f} kg/kNh")