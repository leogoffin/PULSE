import numpy as np
import matplotlib.pyplot as plt
from findCp import findCp

P_generated = 30000000
T_out_C = 1472
m_dot_lb = 164
f = 0.015
eta_is = 0.93

T_out_0 = 817.7777
m_dot_air = 74.389
m_dot = m_dot_air *(1+f)

error = 1
T_in = 800

while(error > 0.00001):
    Cp = findCp((T_in+T_out_0)/2,f)
    T_in_new = P_generated/(m_dot*Cp) + T_out_0
    error = np.abs(T_in-T_in_new)/T_in
    T_in = T_in_new

gamma = Cp/(Cp-287.05)
T_out_is = T_in - (T_in-T_out_0)/eta_is
pi_tot = (T_in/T_out_is)**(gamma/(gamma-1))
eta_p = eta_p = np.log(1 - eta_is*(1 - pi_tot**(-(gamma-1)/gamma))) / np.log(pi_tot**(-(gamma-1)/gamma))
print(eta_p)
print(pi_tot)

