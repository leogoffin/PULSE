import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

alt = 35000 * ft_to_m
M = 0.8
#n/n* = 0.9

BPR = 6
mdot_ref = 485


mdot_ratios = [0.80354, 0.80524, 0.84259, 0.84829, 0.85411, 0.88808, 0.89208, 0.91418, 0.92323, 0.92789, 0.93022, 0.93255, 0.93959, 0.94443]
pi_f = [1.62579, 1.62569, 1.62243, 1.61954, 1.61529, 1.6, 1.59749, 1.57742, 1.5652, 1.55497, 1.54359, 1.53439, 1.50394, 1.44982]
eta_p = [0.712973, 0.715436, 0.76635, 0.772995, 0.779305, 0.8125, 0.816895, 0.843755, 0.850232, 0.850359, 0.848718, 0.843809, 0.836905, 0.837882]

pi_f = interpolate(mdot_ratios,pi_f)
eta_p= interpolate(mdot_ratios,eta_p)
mdot_ratios = np.linspace(mdot_ratios[0],mdot_ratios[-1],100)


Pref = Air()
Pref.set_isa(0)

v = speed_of_sound(alt)*M
P1 = Air(v)
P1.set_isa(alt)

Tmax = 0
m_ratio_opt = 0
mdot_opt = 0
A_opt = 0
P_opt = 0
for m_ratio in mdot_ratios : 
    mdot = true_massflow_from_corrected(m_ratio*mdot_ref,P1,Pref)*BPR/(1+BPR)
    eta_is = eta_is_from_eta_poly(pi_f(m_ratio),eta_p(m_ratio))
    P2,Pc,cp = compr.with_eta(P1,mdot,pi_f(m_ratio),eta_is)
    P3,choked,NPR,A = adiabatic_convergent_nozzle(P2,mdot,0,P1.p)
    T = simple_thrust(mdot,0,P1,P3,A)
    if T > Tmax : 
        Tmax = T
        m_ratio_opt = m_ratio
        mdot_opt = mdot * (BPR+1)/BPR
        A_opt = A
        P_opt = Pc

print(f"Fan Power                      : {P_opt/1e6:.2f} MW")
print(f"Thrust                         : {Tmax:.0f} N")
print(f"Exhaust area                   : {A_opt:.3f} m²")   
print(f"Best relative design flow rate : {m_ratio_opt:.3f}")
print(f"Optimal True mas flow rate     : {mdot_opt:.3f} kg/s")
print('-'*50)

P1 = Pref
def error(m_ratio):
    mdot = m_ratio*mdot_ref *BPR/(1+BPR)
    eta_is = eta_is_from_eta_poly(pi_f(m_ratio),eta_p(m_ratio))
    P2,Pc,cp = compr.with_eta(P1,mdot,pi_f(m_ratio),eta_is)
    P3,choked,NPR,A = adiabatic_convergent_nozzle(P2,mdot,0,P1.p)
    return A - A_opt

m_ratio,i = bisection(error,mdot_ratios[0],mdot_ratios[-1])


mdot = m_ratio*mdot_ref *BPR/(1+BPR)
eta_is = eta_is_from_eta_poly(pi_f(m_ratio),eta_p(m_ratio))
P2,Pc,cp = compr.with_eta(P1,mdot,pi_f(m_ratio),eta_is)
P3,choked,NPR,A = adiabatic_convergent_nozzle(P2,mdot,0,P1.p)
T = simple_thrust(mdot,0,P1,P3,A)

print(f"Fan Power                 : {Pc/1e6:.2f} MW")
print(f"Thrust                    : {T:.0f} N")
print(f"Exhaust area              : {A:.3f} m²")   
print(f"Relative design flow rate : {m_ratio:.3f}")
print(f"True mas  flow rate       : {mdot:.3f} kg/s")