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

mdot_corrs = [2.5,2.7,2.8,2.9,3]
pi_fs = [1.19,1.19,1.18,1.175,1.15]
eta_ss = [0.62,0.75,0.82,0.83,0.75]

interp_pi = interp1d(mdot_corrs,pi_fs, kind='linear', bounds_error=False, fill_value='extrapolate')
interp_eta = interp1d(mdot_corrs,eta_ss, kind='linear', bounds_error=False, fill_value='extrapolate')
mdot_trials = np.linspace(min(mdot_corrs),max(mdot_corrs),100)

T = np.zeros(len(mdot_trials))
pi_s = np.zeros(len(mdot_trials))
eta_iss = np.zeros(len(mdot_trials))
As = np.zeros(len(mdot_trials))
for i,mdot_corr in enumerate(mdot_trials) : 
    mdot = true_massflow_from_corrected(mdot_corr,P1,Pref)
    pi_f = interp_pi(mdot_corr)
    eta_is = interp_eta(mdot_corr)

    P2s,Pfs,Cp12s = compr.with_eta(P1,mdot,pi_f,eta_is)
    P3,choked,NPR,A = adiabatic_convergent_nozzle(P2s,mdot,0,P0.p)
    
    T[i] = simple_thrust(mdot,0,P0,P3,A=A)
    eta_iss[i] = eta_is
    pi_s[i] = pi_f 
    As[i] = A


index = np.argmin(abs(T-71))
print(f"Thrust : {T[index]:.1f} N")
print(f"Corrected mass flow rate : {mdot_trials[index]:.3f}")
print(f"Pressure ratio : {pi_s[index]:.2f}")
print(f"Nozzle outlet area : {As[index]*10000:.4f} cm²")