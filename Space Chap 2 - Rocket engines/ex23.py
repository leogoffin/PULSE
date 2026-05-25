import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

a_ratio = 34.34
a_ratios = 80
delta_hf = 50e6

MR = 3.6
MO2 = 32e-3
MCH4 = 16e-3
MCO2 = 44e-3
MH2O = 18e-3

ptot = 300e5
gamma = 1.28
K = 2.5

MR_star = stoich_mixture_ratio(MCH4,1,MO2,2)
res = combustion_enthalpy(1,delta_hf,MR,MR_star)
delta_h = res["specific_enthalpy"]

res = combustion_analysis(["CH4","O2"],[1,2],[MCH4,MO2],[1,MR],
                          ["CO2","H2O"],[1,2],[MCO2,MH2O])

Pref = Air()
Pref.set_isa(0)


Pr = Air()
Pr.p0 = ptot
Pr.R = res["R"]
Pr.gamma = gamma
cp = get_cp(gamma,Pr.R)
Pr.cp = cp

dT = delta_h/cp
Pr.T0 = 288.15 + dT


c_star = charac_vel(Pr)

Me = M_nozzle2(Pr,a_ratio,Mmax = 20)
Pe = exhaust_cds(Pr,Me)
p_ratio = ptot/Pe.p
cTg = get_cT(p_ratio,a_ratio,gamma,Me,Pref.p/Pe.p)
Ispg = get_Isp_with_cT_star(cTg,c_star)

Me = M_nozzle2(Pr,a_ratios,Mmax = 20)
Pe = exhaust_cds(Pr,Me)
p_ratio = ptot/Pe.p
cTs = get_cT(p_ratio,a_ratios,gamma,Pe.M,0)
Isps = get_Isp_with_cT_star(cTs,c_star)

print(f"R                         : {Pr.R:.0f} J/kgK")
print(f"Cp                        : {cp:.0f} J/kg.K")
print(f"Char vel                  : {c_star:.0f} m/s")
print(f"Thrust coeff (ground)     : {cTg:.2f}")
print(f"Specific impulse (ground) : {Ispg:.0f} s")
print(f"Thrust coeff (space)      : {cTs:.2f} ")
print(f"Specific impulse (space)  : {Isps:.0f} s")
print(f"")

# 3. Thrust min = > a la plus basse altitude possible (sea level)
# la pression du resérvoir diminue jusqu à sa limite max : cas critique r = 0.85 re => pn prend l'aire qui correspond
# nouveau A_ratio 
# Nouvelle inconnue est pression totale, utiliser pressure ratio = F(M) àpd du mach obtenu

new_A_ratio = a_ratio * (0.85**2)
p_s =  Pref.p / K # On a la pression statique avec le critère
Me = M_nozzle2(0,new_A_ratio, gamma)

p_tot = TotalPressure(p_s, gamma, Me)

__, C_T_min, pe = IterateMeTotalCds_FindThrust(gamma, Pr.R, p_tot,Pr.T0, new_A_ratio, 0, Pref.p)

print(C_T_min/cTg)