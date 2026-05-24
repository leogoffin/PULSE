import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

gamma = 1.25

Pe10 = Air()
Pe10.p0 = 60e5
Pe10.gamma = 1.25

Pa10 = Air()
Pa10.set_isa(10000)
Pa10.gamma = 1.25

pi10 = Pe10.p0/Pa10.p
Pe10.M = Mach_from_pratio(pi10,Pa10.gamma)
a_ratio10 = A_ratio_nozzle(Pe10)
CT = get_cT(pi10,a_ratio10,Pe10.gamma,Pe10.M,Pa10.p)

print("-"*50)
print(f"10 km Altitude")
print(f"Exhaust Mach : {Pe10.M:.2f}")
print(f"Expansion ratio : {a_ratio10:.2f}")
print(f"Thrust Coefficient : {CT:.3f}")

Pe20 = Air()
Pe20.p0 = 60e5
Pe20.gamma = 1.25

Pa20 = Air()
Pa20.set_isa(20000)
Pa20.gamma = 1.25

pi20 = Pe20.p0/Pa20.p
Pe20.M = Mach_from_pratio(pi20,Pa20.gamma)
a_ratio20 = A_ratio_nozzle(Pe20)
CT = get_cT(pi20,a_ratio20,Pe20.gamma,Pe20.M,Pa20.p)

print("-"*50)
print(f"20 km Altitude")
print(f"Exhaust Mach : {Pe20.M:.2f}")
print(f"Expansion ratio : {a_ratio20:.2f}")
print(f"Thrust Coefficient : {CT:.3f}")

pi_switch = pi_switch_nozzle_(pi10,pi20,a_ratio10,a_ratio20,Pe10.M,Pe20.M,gamma)
p_switch = Pe20.p0/pi_switch
def error(alt):
    return pressure(alt) - p_switch
alt_switch,i = bisection(error,10000,20000)
print(f"Switching altitude : {alt_switch/1e3:.3f} km")
p = pressure(alt_switch)
crit = p/Pe20.p0 * pi20
print(f"Summerfield criterion : {crit:.0f} < K")