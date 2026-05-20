import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

alt = 30000 * ft_to_m
M = 0.5
n = 37780 #RPM

Pref = Air()
Pref.set_isa(0)
Pref.set_total()

v = M*speed_of_sound(alt)
P0 = Air(v)
P0.set_isa(alt)
P0.set_total()

P1 = P0.RAM()

n_corr = corrected_rpm(n,P1,Pref)
print(f"Corrected rotation speed : {n_corr:.0f} RPM")