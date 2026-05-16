import os
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))
from Importer import *


pi_c_a = 4.5
eta_is_c = 0.85

p_dis_b = 90 * psi_to_Pa
T0_dis_b = F_to_K(480)
v_a = 500*ft_to_m
T_a = F_to_K(60)
p_a = 14.7 *psi_to_Pa

def polytropic_eff(pi_c,eta,gamma = 1.4,Ta = None,T0b = None):
    num = np.log10(pi_c**((gamma-1)/gamma))
    den = np.log10((pi_c**((gamma-1)/gamma)-1)/eta + 1)
    return num/den

def polytropic_eff_discharge(P_a,P_dis,gamma = 1.4):
    num = np.log10((P_dis.p0/P_a.p)**((gamma-1)/gamma))
    den = np.log10(P_dis.T0/P_a.T)
    return num/den

P_dis = Air(v_a)
P_dis.p = p_dis_b
P_dis.T0 = T0_dis_b
T= T_a
P_dis.gamma = 1.4
diff = np.inf
tol = 1e-5
while diff > tol :
    M = v_a/np.sqrt(P_dis.gamma*P_dis.R*T)
    new_T = P_dis.T0 / (1 + (P_dis.gamma-1)/2 * M**2) 
    diff = abs(new_T-T)/T
    T = new_T

P_dis.M = M
P_dis.T = T
P_dis.set_total()


Pa = Air(v_a)
Pa.p = p_a
Pa.T = T_a

print(f"Compressor A polytropic efficiency : {polytropic_eff(pi_c_a,eta_is_c):.3f}")
print(f"Compressor B polytropic efficiency : {polytropic_eff_discharge(Pa,P_dis):.3f}")


