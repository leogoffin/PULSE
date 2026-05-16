import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

mdot = 68 #kg/s
mdotf = 1.2 #kg/s
f = mdotf/mdot 
mdotf2 = 2.2 #kg/s
f2 = mdotf2/mdot + f 

P0 = Air()
P0.p = 1.033e5 #Pa
P0.T0 = C_to_K(15)

P3 = Air(150)
P3.T0 = C_to_K(890) #K
P3.p0 = 5.47e5 #Pa

P4 = Air(346)
P4.T0 = C_to_K(685) #K
P4.p0 = 2.36e5 #Pa
eta_is_t,Pt,T02_is = turb.get_eta(P3,P4,mdot*(1+f),f)
print(f"Isentropic efficiency : {eta_is_t:.3f}")

P5a = Air()
P5a.p0 = 2.16e5
P5a.T0 = C_to_K(1600)
P5a.set_static()
Pout = 115e3


P6a,chokeda,NPRa,Aa = adiabatic_convergent_nozzle(P5a,mdot,f2,Pout)
print(f"Nozzle choked ? -{chokeda}")
print(f"Exhaust velocity : {P6a.v:.0f} m/s\n")


P6,choked,NPR,A = adiabatic_convergent_nozzle(P4,mdot,f,P0.p)
Tp,Ts = turbofan_thrust(mdot,0,f,P0,P6,AP=A)

Ta = 58680
print(f"Thurst (dry) : {Tp/1e3:.2f} kN")
print(f"Thurst (wet) : {Ta/1e3:.2f} kN")
print(f"Thrust ratio : {Ta/Tp:.2f} kN\n")

SFCa = mdot*f2/Ta
SFC = mdot*f/Tp
print(f"SFC (dry) : {SFC*10*3600:.3f} kg/DaNh")
print(f"SFC (wet) : {SFCa*10*3600:.3f} kg/DaNh")
print(f"SFC ratio : {SFCa/SFC:.2f}")