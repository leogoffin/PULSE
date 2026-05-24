import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *

pi_fan_star = 1.54
eta_f_is_star = 0.92
alt = 30000*ft_to_m
M = 0.85
Ga_star = 192

G_aOverG_a_star_points = [ 0.894791 ,0.959345 ,1.007327 ,1.032927, 1.050198]
PiOverPi_star_points = [1.035155, 1.028333, 1.018101, 1.00062, 0.969922]

PiOverPi_star_points_2 = [ 0.971654, 0.992714 ,1.005178, 1.017584, 1.02435, 1.030144,1.034616]
EtaOverEtastar_points_2 = [0.905892, 1, 1.008299, 1.002905, 0.981743 ,0.93361 ,0.877178]


Pref = Air(0)
Pref.set_isa(0)
Pref.set_total()

v = M * speed_of_sound(alt)
P0 = Air(v)
P0.set_isa(alt)
P0.set_total()

P1 = P0.RAM()


# We sweep over the table  G_aOverG_a_star_points and compute the propulsive power for ecah possible case
x = np.linspace(np.min(G_aOverG_a_star_points), np.max(G_aOverG_a_star_points), 100) 

f_pi = interp1d(G_aOverG_a_star_points, PiOverPi_star_points, kind='linear', bounds_error=False, fill_value='extrapolate')
f_eta  = interp1d(PiOverPi_star_points_2, EtaOverEtastar_points_2,  kind='linear', bounds_error=False, fill_value='extrapolate')
 
eta_prop = np.zeros(len(x))
eta_th = np.zeros(len(x))
A_ratio = np.zeros(len(x))
A1 = np.zeros(len(x))
T = np.zeros(len(x))

for i in range (len(x)):
    pi_fs  = float(f_pi(x[i]))
    eta_fs = float(f_eta(pi_fs))
    pi_fs_real = pi_fs * pi_fan_star 

    mdot = true_massflow_from_corrected(x[i]*Ga_star,P1,Pref)
    P2,Pf,Cp12 = compr.with_eta(P1,mdot,pi_fs_real,eta_fs)

    P3,choked,NPR,A = adiabatic_convergent_nozzle(P2,mdot,0,P0.p)
   
    T[i] = simple_thrust(mdot,0,P0,P3,A)
    eta_prop [i] = PropulsiveEff(mdot, 0, P3.v,P0.v , T[i])
    eta_th [i] = ThermalEff_from_Power(mdot, P3.v, P0.v, Pf)
    A1 [i] = mdot / (P0.rho * P0.v) 
    A_ratio[i] = A/A1[i]

index = np.argmax(eta_prop)
A1_opt = A1[index]
eta_overall = overalleff(eta_th[index], eta_prop[index])
m_star_opt = x[index] * Ga_star * A1_opt         
pi_opt = f_pi(x[index]) * pi_fan_star 
eta_is_opt = f_eta(f_pi(x[index])) * eta_f_is_star  


print(f"Area ratio {A_ratio[index]:.3f}")
print(f'Fan operating point:\n  m_dot_star: {m_star_opt:.0f}\n  pi_f_star: {pi_opt:.3f}\n eta_s_star: {eta_is_opt:.3f}')
print(f'In-flight values:\n T/A1: {T[index]/A1_opt:.0f} \n m/A1: {x[index] * Ga_star:.0f} \n Overall eff: {eta_overall:.3f}')

