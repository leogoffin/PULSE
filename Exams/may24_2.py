import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

MR = 4
Cp = 2610
MH2 = 1e-3*2 #kg/mole
MO2 = 16e-3*2 #kg/mole
MH2O = (MO2 + 2*MH2)/2

H2_mf = 1/(1+MR)
O2_mf = 1*MR/(1+MR) 

moles_H2 = H2_mf / MH2
moles_O2 = O2_mf / MO2

H2_MF = moles_H2/(moles_H2+moles_O2)
O2_MF = moles_O2/(moles_H2+moles_O2)


results = combustion_analysis(
    species_in=["H2", "O2"],nu_in=[2, 1],M_in=[MH2, MO2],quantities_in=[H2_mf, O2_mf],
    species_out=["H2O"],nu_out=[2],M_out=[MH2O],
)

MR_star = stoich_mixture_ratio(2,2,32,1)
print(1/2*MO2/MH2)



mdot = 1
res = combustion_enthalpy(mdot,120e6,MR,MR_star,print_results=True)
dT = res["energy_released"]/Cp/mdot
print(dT+ 288.15)