import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

Pa10 = Air(0)
Pa10.set_isa(10000)
Pa10.set_total()

Pa20 = Air(0)
Pa20.set_isa(20000)
Pa20.set_total()

Pe10 = Air()
Pe10.gamma = 1.25
Pe10.p0 = 60e5
pi10 = Pe10.p0/Pa10.p
Pe10.M = Mach_from_pratio(pi10,Pe10.gamma)
a_ratio10 = A_ratio_nozzle(Pe10)
CT = get_cT(pi10,a_ratio10,Pe10.gamma,Pe10.M,Pa10.p)
print(Pe10.M)
print(a_ratio10)
print(CT)

Pe20 = Air()
Pe20.gamma = 1.25
Pe20.p0 = 60e5
pi20 = Pe20.p0/Pa20.p
Pe20.M = Mach_from_pratio(pi20,Pe20.gamma)
a_ratio20 = A_ratio_nozzle(Pe20)
CT = get_cT(pi20,a_ratio20,Pe20.gamma,Pe20.M,Pa20.p)
print(Pe20.M)
print(a_ratio20)
print(CT)
print(pi10)
#switch when CT are equal
print(pi_switch_nozzle_(pi10,pi20,a_ratio10,a_ratio20,Pe10.M,Pe20.M,Pe10.gamma))