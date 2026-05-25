import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *



MCH4 = 16e-3 #kg/mole
MO2 = 32e-3 #kg/mole
MCO2 = 44e-3 #kg/mole
MH2O = 18e-3 #kg/mole

mdot_CH4 = 140 #kg/s
mdot_O2 = 650 #kg/s
mdot = mdot_CH4+mdot_O2
MR = mdot_O2/mdot_CH4
MR_star = stoich_mixture_ratio(MCH4,1,MO2,2)
print(MR_star)
delta_hf = 47.1e6
Cp = 1000
p0 = 380e5

res = combustion_analysis(["CH4","O2"],[1,2],[MCH4,MO2],[mdot_CH4,mdot_O2],
                          ["CO2","H2O"],[1,2],[MCO2,MH2O])

res2 = combustion_enthalpy(mdot,delta_hf,MR,MR_star)

dT = res2["specific_enthalpy"] / Cp
Tout = dT + 288.15
print(f"Reservoir temperature : {Tout:.0f} K")

Pr = Air(R = res["R"])
Pr.cp = Cp
Pr.p0 = p0
Pr.T0 = Tout
Pr.get_gamma_from_cp()
At = A_nozzle(Pr,mdot,M=1)
print(f"Throat area : {At:.4f} m²")