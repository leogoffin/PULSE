import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

alt = 35000*ft_to_m #m
v = 590.3 * mph_to_ms
mdot = 77 #kg/s

RR = 0.99
OPR = 13.5 
eta_p_c = 0.98

eta_m = 0.99
TIT = 1210 #K
eta_t_is = 0.98

delta_hf = 42.8e6 #J/kg
eta_cc = 0.97
ploss_cc = 0.03
pi_cc = 1 - ploss_cc

mdotfab = 3.2 #kg/s
eta_ab = 0.97
ploss_ab = 0.05
pi_ab = 1 - ploss_ab

Pref = Air() # sets air  with default R (287.058 kg/JK)
Pref.set_isa(0) # sets the properties from ISA tables at ground

#A. Ambiant air
P0 = Air(v = v) # sets air at velocity v with default R (287.058 kg/JK)
P0.set_isa(alt) # sets the properties from ISA tables at the right altitude

#B. Intake
P1,pi_r = P0.RAM(RR) #cds using RAM effect
print(f"RAM pressure ratio : {pi_r:.3f}")

#C. Compressor
pi_c = OPR/pi_r
print(f"Compressor pressure ratio : {pi_c:.3f}")
P3,Pc,Cp23 = compr.with_eta_poly(P1,mdot,pi_c,eta_p_c)
#print(f"Heat capacity (compr) : {Cp23:.0f} J/kg.K")
print(f"Compressor Work : {Pc/1e6:.2f} MW")

#D. Combustion chamber
P4,mdotfcc,fcc =cc.with_Tout(P3,P0.T,mdot,delta_hf,TIT,eta_cc = eta_cc,pi_cc = pi_cc)

#E. Turbine
Pt = Pc/eta_m
P5,Cp45,T5s = turb.with_P(Pt,mdot,fcc,P4,eta_is = eta_t_is)

#F. Afterburner
fab = mdotfab/mdot
f_tot = fab + fcc

P6w = cc.with_f(P5,Pref.T,mdot,delta_hf,fab,fin = fcc,eta_cc = eta_ab, pi_cc = pi_ab) 

#G Exhaust nozzle

P8d, NPR,Atd = adiabatic_convdiv_nozzle(P5,mdot,f_tot,P0.p)
Td = simple_thrust(mdot,f_tot,P0,P8d,Atd)

P8w, NPR,Atw = adiabatic_convdiv_nozzle(P6w,mdot,f_tot,P0.p)
Tw = simple_thrust(mdot,f_tot,P0,P8w,Atw)

Qd = delta_hf*mdotfcc
Qw = delta_hf*(mdotfcc+mdotfab)

SFCd = Qd/Td
SFCw = Qw/Tw

eta_td = thermal_efficiency(P0.v,P8d.v,mdot,fcc,delta_hf)
eta_tw = thermal_efficiency(P0.v,P8w.v,mdot,f_tot,delta_hf)

eta_pd = propulsive_efficiency(P0.v,P8d.v,mdot,fcc,Td)
eta_pw = propulsive_efficiency(P0.v,P8w.v,mdot,fcc,Tw)

P1.sprint(1)
P3.sprint(3)
P4.sprint(4)
P5.sprint(5)
P8d.sprint(8, "dry")

P6w.sprint(6,"wet")
P8d.sprint(8, "wet")

print('-'*50)

print(f"Thrust (dry)                : {Td:.0f} N")
print(f"Throat area (dry)           : {Atd:.3f} m²\n")
print(f"SFC (dry)                   : {SFCd*kgN_to_kghDaN:.2f} kg/DaN.h")
print(f"Thermal efficiency (dry)    : {eta_td:.3f}")
print(f"Propulsive efficiency (dry) : {eta_pd:.3f}")
print(f"Overall efficiency (dry)    : {eta_pd*eta_td:.3f}")

print('-'*50)

print(f"Thrust (wet)                : {Tw:.0f} N")
print(f"Throat area (wet)           : {Atw:.3f} m²\n")
print(f"SFC (wet)                   : {SFCw*kgN_to_kghDaN:.2f} kg/DaN.h")
print(f"Thermal efficiency (wet)    : {eta_tw:.3f}")
print(f"Propulsive efficiency (wet) : {eta_pw:.3f}")
print(f"Overall efficiency (wet)    : {eta_pw*eta_tw:.3f}")
print('-'*50)