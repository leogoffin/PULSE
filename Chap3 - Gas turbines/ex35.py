import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *


eta_is_c = 0.85
delta_hf = 42.979e6 # J/kg
pi_c = 2.86
P_t2 = 173e3 #W
mdot = 1.75 #kg
TIT = C_to_K(983) # K
SFC = 878/3600/1000/1000 #kg/(sw)
mdot_f = SFC * P_t2
f = mdot_f/mdot

P1 = Air(0)   # no velocity
P1.set_isa(0) # ground stdatm
P1.set_total()
P2,P_c,Cp12 = compr.with_eta(P1,mdot,pi_c,eta_is_c) # compressor

P4,Cp34,T4s = turb.with_P(P_c,mdot*(1+f),f,Pin = P2) # turbine
eta_cc = cc.get_eta(P2,P1.T0,mdot,f,delta_hf,TIT)
print(eta_cc)
print(P4.T0)



