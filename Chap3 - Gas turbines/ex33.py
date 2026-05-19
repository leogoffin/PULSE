import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *
import matplotlib.pyplot as plt

Pt = 30e6 #W
Tout = R_to_K(1472) #K
mdot = 164 * lb_to_kg
f = 0.015
eta_is_t = 0.93


P1,Cp,T1is = turb.rev_with_P(Pt,mdot,f,Tout = Tout,eta_is = eta_is_t)
gamma = get_gamma(Cp)
eta_p = turb.get_eta_poly(4.683,eta_is_t,gamma)
print(f"Polytropic efficiency : {eta_p:.3f}")