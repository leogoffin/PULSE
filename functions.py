import numpy as np 
from scipy.interpolate import interp1d
from findCp import findCp
from Convert_unit_and_ISA import *


''' Total/Static Values ''' 

def TotalTemp(T, gamma, M):
    T_tot = T * (1 + M**2 * (gamma-1)/2)
    return T_tot

def TotalPressure(P, gamma, M):
    P_tot = P * (1 + M**2 * (gamma-1)/2)**(gamma/(gamma-1))
    return P_tot

def Static_Temp(T_tot, gamma, M):
    T = T_tot / (1 + M**2 * (gamma-1)/2)
    return T

def StaticPressure(P_tot, gamma, M):
    P = P_tot / (1 + M**2 * (gamma-1)/2)**(gamma/(gamma-1))
    return P

def TotalTemp_gamma(T, f, M):
    R = 287.058
    Cp = findCp(T, f)
    gamma = Cp / (Cp-R)
    T_tot = T * (1 + M**2 * (gamma-1)/2)
    return T_tot, gamma

def findRho(p, T):
    R = 287.058
    rho = p / (R * T)
    return rho


''' Efficiencies '''

def IsoEff_Comp(gamma, pi, T1_tot, T2_tot):
    eta = (pi**((gamma-1)/gamma) - 1) / (T2_tot/T1_tot - 1)
    return eta

def PolyEff_Comp(gamma, pi_c, T_1_tot, T_2_tot):
    eta_p_c = (gamma - 1) * np.log(pi_c) / (gamma * np.log(T_2_tot/T_1_tot))
    return eta_p_c

def PolyEff_from_IsoEff_Comp(isoEff, gamma, pi):
    PolyEff = (gamma - 1) / (gamma * np.log((pi**((gamma-1)/gamma) - 1) / isoEff + 1) / np.log(pi))
    return PolyEff

def PolyEff_from_IsoEff_Turb(isoEff, gamma, pi):
    k = (gamma - 1) / gamma
    PolyEff = -np.log(1 - isoEff * (1 - pi**(-k))) / (k * np.log(pi))    
    return PolyEff

def IsoEffStageComp_from_PolyEff(poly_eff_stage, gamma, pi_stage):
    eta_is_stage = (pi_stage**((gamma-1)/gamma) - 1  ) / (pi_stage**((gamma-1)/(gamma*poly_eff_stage)) - 1)
    return eta_is_stage

def ThermalEff(m_dot_a, m_dot_f, v_j, v_1, delta_hf):
    """Computes thermal efficiency from mass flow rates, exit velocity vj and inlet velocity v_1"""
    eta_th = ((m_dot_a + m_dot_f) * v_j**2 - m_dot_a * v_1**2 ) / (2*m_dot_f*delta_hf)
    return eta_th

def ThermalEff_from_PowerConsumed(m_dot_a, v_j, v_1, P):
    eta_th = 0.5 * m_dot_a * (v_j**2 - v_1**2) / P
    return eta_th

def PropulsiveEff(m_dot_a, m_dot_f, v_j, v_1, T):
    """Computes propulsive efficiency from mass flow rates, exit velocity vj, inlet velocity v_1 and thrust T"""
    eta_p = (2*T*v_1) / ((m_dot_a + m_dot_f) * v_j**2 - m_dot_a * v_1**2 )
    return eta_p

def OverallEff(eta_th, eta_prop):
    return eta_th * eta_prop

''' Similiarity'''
# Pressure and temperature ratios are maintained : p2tot/p1tot = p2tot_star/p1tot_star

def m_dot_a_FromCorrectedCds(m_dot_a_star, p_1_tot, p_1_tot_star, T_1_tot, T_1_tot_star):
    """Stars denote standars cds"""
    m_dot_a = m_dot_a_star * p_1_tot * np.sqrt(T_1_tot_star/T_1_tot) / p_1_tot_star
    return m_dot_a

def n_FromCorrectedCds(n_star, T_1_tot, T_1_tot_star):
    n = n_star * np.sqrt(T_1_tot/T_1_tot_star)
    return n

def Pc_fromCorrectedCds(Pc_star, p_1_tot, p_1_tot_star, T_1_tot, T_1_tot_star):
    Pc = Pc_star * p_1_tot * np.sqrt(T_1_tot/T_1_tot_star) / p_1_tot_star
    return Pc 

def Uncorrected_m_dot(m_dot_star, p_1_tot, T_1_tot):
    p_r_tot = 101325
    T_r_tot =  288.15
    m_dot = m_dot_star * p_1_tot / p_r_tot * np.sqrt(T_r_tot/T_1_tot)
    return m_dot

def nStar(n, T_1_tot):
    T_r = 288.15 
    n_star = n * np.sqrt(T_r/T_1_tot)
    return n_star

'''Mass flow rates'''

def Primary_Secondary_M_dot(alpha, m_dot_air):
    m_dot_p =  m_dot_air / (1 + alpha)
    m_dot_s = alpha* m_dot_air / (1 + alpha)
    return m_dot_p, m_dot_s

''' Sound Speed, Mach, Gamma '''

def SpeedSound (gamma, T):
    R = 287.058 
    return (gamma * R * T)**0.5

def findGamma(Cp): 
    R = 287.058
    gamma = Cp / (Cp-R)
    return gamma

def Cp_from_Gamma(gamma):
    R = 287.058 
    Cp = gamma / (gamma - 1) * R
    return Cp

def Vel_from_Mach(a, M):
    return a*M


''' Compressor/fan'''

def IterateM_static_T(T_guess_init, gamma, T_tot, V):
    """Used to compute M from gamma, velocity and total_temp"""
    error = 1 
    R = 287.058 
    T_guess_new = T_guess_init
    while error > 1e-5:
        T_guess = T_guess_new
        a = np.sqrt(gamma*R*T_guess)
        M = V / a 
        T_guess_new = Static_Temp(T_tot, gamma, M)
        error = np.abs((T_guess_new- T_guess)/T_guess)
    return M

def IterateGamma_Comp(gamma_init, T1_tot, pi, eta_is_init, eta_prop):
    """Set eta_is or eta_prop to 0 depending on the given value"""
    T1_tot = T1_tot
    R = 287.058 
    error = 1
    T_2s_new = T1_tot * pi**((gamma_init - 1)/ gamma_init)
    if (eta_is_init == 0):
            eta_is_new = (pi**((gamma_init-1)/gamma_init) - 1) / (pi**((gamma_init-1)/(gamma_init*eta_prop)) - 1)
    else:
        eta_is_new = eta_is_init

    while error > 1e-5:
        T_2s = T_2s_new
        eta_is = eta_is_new
        T2_tot = (T_2s - T1_tot)/ eta_is + T1_tot
        Cp = findCp((T2_tot+T1_tot)/2, 0)
        gamma  = Cp / (Cp - R)
        T_2s_new = T1_tot * pi **((gamma - 1)/ gamma)
        error = np.abs((T_2s_new - T_2s)/T_2s)
        if(eta_is_init == 0):
            eta_is_new = (pi**((gamma-1)/gamma) - 1) / (pi**((gamma-1)/(gamma*eta_prop)) - 1)
        else:
            eta_is_new = eta_is_init

    if(eta_is_init == 0):
        return gamma, Cp, T2_tot, eta_is
    else:
        return gamma, Cp, T2_tot 
    
def Iterate_Tout_Comp_from_TIT(gamma2r_init, m_dot_a, m_dot_f, delta_hf, TIT):
    R = 287.058 
    error = 1
    T_r_tot = 288.15
    Cp3r = findCp((TIT+T_r_tot)/2, m_dot_f/m_dot_a)
    m_dot_e = m_dot_a + m_dot_f

    # Initialisation
    gamma2r_new = gamma2r_init
    Cp2r = gamma2r_new / (gamma2r_new - 1) * R
    
    while error > 1e-5:
        gamma2r = gamma2r_new
        T_2_tot = T_r_tot + (m_dot_e * Cp3r * (TIT - T_r_tot) - m_dot_f * delta_hf) / (m_dot_a * Cp2r)
        Cp2r = findCp((T_2_tot+T_r_tot)/2, 0)
        gamma2r_new  = Cp2r / (Cp2r - R)
        error = np.abs((gamma2r_new - gamma2r)/gamma2r)

    return gamma2r, Cp2r, T_2_tot

def PowerComp(m_dot_a, Cp, T_1_tot, T_2_tot):
    Pc = m_dot_a * Cp * (T_2_tot - T_1_tot)
    return Pc

def Iterate_lpc(pi_tot, T_1_tot, p_1_tot, TIT, p_4_tot, m_dot_p, m_dot_a, m_dot_f, Cp_air, gamma_air,
                Cp13, gamma13, eta_p_c, eta_s_t, lambda_sec_duct):
    """
    Exercise 4.1.9
    Finds Pi_LPC such that the energy balance is satisfied AND the mixer total pressure equality (4.88) holds simultaneously.

    The constraint merged into a single residual:
        R(Pi_LPC) = p5_primary - p5_secondary = 0

    where:
        p5_secondary = 0.99 * Pi_LPC * p1          
        p5_primary   = p4 * (Wt / (eta_s_t*Cp_45*T4) + 1)^(gamma45/(gamma45-1))
        Wt from energy balance = [(1+B)*W_LPC + W_HPC] / (1+f)
    """
    tol = 1e-8
    max_iter = 100
    f = m_dot_f / m_dot_p

    def residual(Pi_LPC):
        # Secondary duct total pressure at mixer inlet
        p5s_tot = (1 - lambda_sec_duct) * Pi_LPC * p_1_tot

        # T2: cold air properties (secondary flow assumption, gamma_air constant)
        T_2_tot = T_1_tot * Pi_LPC ** ((gamma_air - 1) / (gamma_air * eta_p_c))

        # T3: primary flow, single mean Cp13/gamma13 over LPC+HPC
        Pi_HPC = pi_tot / Pi_LPC
        T_3_tot = T_2_tot * Pi_HPC ** ((gamma13 - 1) * eta_p_c / gamma13)

        # Compressor specific works for energy balance (Cp13 for primary flow)
        W_LPC = Cp13 * (T_2_tot - T_1_tot)
        W_HPC = Cp13 * (T_3_tot - T_2_tot)

        # Required turbine specific work from energy balance (4.90)
        W_t = (m_dot_a * W_LPC + m_dot_p * W_HPC) / ((1 + f) * m_dot_p)

        # Hot gas properties at mean turbine temperature
        # First estimate of T5 using Cp at TIT
        Cp45_est = findCp(TIT, f)
        T_5_est = TIT - W_t / Cp45_est
        # Corrected Cp45 at mean turbine temperature
        Cp45 = findCp((TIT + T_5_est) / 2, f)
        gamma_45 = findGamma(Cp45)

        # Invert eq. 4.93: Wt = eta_s_t * Cp45 * T4 * [1 - (p5/p4)^((gamma-1)/gamma)]
        bracket = 1.0 - W_t / (eta_s_t * Cp45 * TIT)
        if bracket <= 0:
            return np.nan
        p_5_tot = p_4_tot * bracket ** (gamma_45 / (gamma_45 - 1))

        return (p_5_tot - p5s_tot) / p_1_tot

    # Secant method, starting from the two extrema
    Pi_a = 1.0 + 1e-6
    Pi_b = pi_tot

    R_a = residual(Pi_a)
    R_b = residual(Pi_b)

    if R_a * R_b > 0:
        raise ValueError(
            f"Residual does not change sign: R(1)={R_a:.4f}, R(Pi_c={pi_tot:.2f})={R_b:.4f}. "
            "Check input parameters."
        )

    for i in range(max_iter):
        if abs(R_b - R_a) < 1e-15:
            raise RuntimeError("Secant method stalled: denominator ≈ 0.")
        Pi_new = Pi_b - R_b * (Pi_b - Pi_a) / (R_b - R_a)
        Pi_new = np.clip(Pi_new, 1.0 + 1e-6, pi_tot)
        R_new = residual(Pi_new)
        Pi_a, R_a = Pi_b, R_b
        Pi_b, R_b = Pi_new, R_new
        if abs(R_b) < tol:
            break
    else:
        raise RuntimeError(f"Secant method did not converge in {max_iter} iterations.")

    # Recompute all quantities at converged Pi_LPC
    Pi_LPC = Pi_b
    Pi_HPC = pi_tot / Pi_LPC

    T_2_tot = T_1_tot * Pi_LPC ** ((gamma_air - 1) / (gamma_air * eta_p_c))
    T_3_tot = T_2_tot * Pi_HPC ** ((gamma13 - 1) * eta_p_c / gamma13)

    W_LPC = Cp13 * (T_2_tot - T_1_tot)
    W_HPC = Cp13 * (T_3_tot - T_2_tot)
    W_t   = (m_dot_a * W_LPC + m_dot_p * W_HPC) / ((1 + f) * m_dot_p)

    Cp45_est = findCp(TIT, f)
    T_5_est  = TIT - W_t / Cp45_est
    Cp45 = findCp((TIT + T_5_est) / 2, f)
    gamma_45 = findGamma(Cp45)
    T_5_tot  = TIT - W_t / Cp45

    p_5s_tot = (1 - lambda_sec_duct) * Pi_LPC * p_1_tot

    return Pi_LPC, Pi_HPC, T_2_tot, T_3_tot, T_5_tot, p_5s_tot

''' Combustion chamber '''

def IterateGamma_CC (gamma_3r_init, m_dot_a, m_dot_f, delta_Hf, eta_cc, T_2_tot):
    R = 287.058 
    error = 1
    T_r_tot = 288.15
    gamma3r_new = gamma_3r_init
    Cp3r = gamma3r_new / (gamma3r_new - 1) * R

    f = m_dot_f/m_dot_a
    m_dot_e = m_dot_a + m_dot_f
    Cp2r = findCp((T_2_tot+ T_r_tot)/2,0)  #ATTENTION: f = 0 pour la référence
    
    while error > 1e-5:
        gamma3r = gamma3r_new
        #print(gamma3r_new)
        T_3_tot =(m_dot_a * Cp2r * (T_2_tot-T_r_tot) + eta_cc * (delta_Hf*1000) * m_dot_f) / (m_dot_e * Cp3r) + T_r_tot #fuel heating value in J
        Cp3r = findCp((T_3_tot + T_r_tot)/2, f) #We are after the combustion chamber so we use f
        gamma3r_new  = Cp3r / (Cp3r - R)
        error = np.abs((gamma3r_new - gamma3r)/gamma3r)
    
    return gamma3r, Cp3r, T_3_tot

def IterateMF_CC (m_dot_a, m_dot_f_init, delta_Hf, eta_cc, T_2_tot, T_3_tot):
    error = 1
    T_r_tot = 288.15
    m_dot_f_new = m_dot_f_init
    Cp2r = findCp((T_2_tot+ T_r_tot)/2,0)  #ATTENTION: f = 0 pour la référence
    while error > 1e-5:
        m_dot_f = m_dot_f_new 
        Cp3r = findCp((T_3_tot + T_r_tot)/2, m_dot_f / m_dot_a) #We are after the combustion chamber so we use f
        m_dot_f_new = m_dot_a*(Cp3r*(T_3_tot - T_r_tot) - Cp2r*(T_2_tot - T_r_tot))/ (eta_cc*delta_Hf - Cp3r*(T_3_tot-T_r_tot))
        error = np.abs((m_dot_f_new - m_dot_f)/m_dot_f)
    return m_dot_f, Cp3r

def PressureAfterLossesCC(p_3_tot, lambda_cc):
    p_4_tot = (1-lambda_cc)*p_3_tot
    return p_4_tot


''' Turbine '''

def IterateTinlet_Turb(gamma_init, P_t, T_4_tot, eta_m, m_dot_a, m_dot_f):

    R = 287.058 
    error = 1
    gamma34_new = gamma_init
    Cp34 = gamma34_new / (gamma34_new - 1) * R
    f = m_dot_f / m_dot_a

    while error > 1e-5:
        gamma34 = gamma34_new
        T_3_tot = (P_t) / (eta_m * (m_dot_a + m_dot_f) * Cp34) + T_4_tot
        Cp34 = findCp((T_4_tot + T_3_tot)/2, f) #We are after the combustion chamber so we use f
        gamma34_new  = Cp34 / (Cp34 - R)
        error = np.abs((gamma34_new - gamma34)/gamma34)
    return T_3_tot, gamma34 

def IterateGamma_Turb (gamma_init, T_4_tot, P_c, P_f,  eta_m, m_dot_a, m_dot_f):
    """Set Pf to 0 if no losses/fan"""
    R = 287.058 
    error = 1
    gamma45_new = gamma_init
    Cp45 = gamma45_new / (gamma45_new - 1) * R
    f = m_dot_f / m_dot_a

    while error > 1e-5:
        gamma45 = gamma45_new
        T_5_tot = - (P_c + P_f) / (eta_m * (m_dot_a + m_dot_f) * Cp45) + T_4_tot
        Cp45 = findCp((T_4_tot + T_5_tot)/2, f) #We are after the combustion chamber so we use f
        gamma45_new  = Cp45 / (Cp45 - R)
        error = np.abs((gamma45_new - gamma45)/gamma45)
    
    return gamma45, Cp45, T_5_tot 

def P_Out_Turb_PolyEff(T_5_tot, T_4_tot, p_4_tot, gamma45, eta_poly):
    p_5_tot = p_4_tot * (T_5_tot/T_4_tot)**(gamma45/((gamma45-1)*eta_poly))
    return p_5_tot

''' Afterburner '''

def IterateMF_AB (m_dot_a, m_dot_f, m_dot_ab_init, delta_Hf, eta_ab, T_4_tot, T_5_tot):
    error = 1
    T_r_tot = 288.15
    m_dot_ab_new = m_dot_ab_init
    Cp4r = findCp((T_4_tot+ T_r_tot)/2, m_dot_f / m_dot_a)  #ATTENTION: f not 0 at the reference since we are after the CC
    while error > 1e-5:
        m_dot_ab = m_dot_ab_new 
        Cp5r = findCp((T_5_tot + T_r_tot)/2, (m_dot_ab + m_dot_f )/ m_dot_a) #We are after the AB so we use m_dot_f + m_dot_ab
        m_dot_ab_new = (m_dot_a + m_dot_f)*(Cp5r*(T_5_tot - T_r_tot) - Cp4r*(T_4_tot - T_r_tot))/ (eta_ab*delta_Hf - Cp5r*(T_5_tot-T_r_tot))
        error = np.abs((m_dot_ab_new - m_dot_ab)/m_dot_ab)
    return m_dot_ab, Cp5r


''' Nozzle '''

def NozzlePressureRatio(gamma):
    """Return NPR_star, to determine whether the nozzle is choked or not"""
    pi_n = (1 + (gamma - 1)/2)**(gamma/(gamma-1))
    return pi_n

def ExhaustArea(m_dot_e, v_e, T_5, p_5):
    R = 287.058
    rho = p_5 / (R * T_5)
    A5 = m_dot_e / (rho * v_e)
    return A5

def m_dot_from_ExhaustAera(A, v, rho):
    m_dot = rho * v * A
    return m_dot

def M_ConvDivNozzle(gamma, p5_tot, p5):
     M = np.sqrt(2 * ((p5_tot/p5)**((gamma-1)/gamma) - 1 ) / (gamma - 1))
     return M

def MachAdaptedNozzle(p_8_tot, p_8, gamma):
    M8 = np.sqrt( 2 * ((p_8_tot/p_8)**((gamma - 1)/gamma) - 1) / (gamma - 1))
    return M8

def IterateGamma_ConvNozzle(gamma_init, T_4_tot, T_5_tot, m_dot_a, m_dot_f):
    R = 287.058 
    f = m_dot_f / m_dot_a
    
    error = 1
    gamma_new = gamma_init
    while error > 1e-5:
        gamma = gamma_new
        T_5 = Static_Temp(T_5_tot, gamma, 1)
        Cp = findCp((T_4_tot + T_5)/2, f) # Heat capacity is computed over the expansion
        gamma_new  = Cp/ (Cp - R)
        error = np.abs((gamma_new - gamma)/gamma)
    return gamma, T_5

def IterateGammaM_ConvDivNozzle(gamma_init, T_4_tot, T_5_tot, m_dot_a, m_dot_f, p5, p5_tot):
    R = 287.058 
    f = m_dot_f / m_dot_a
    error = 1
    gamma_new = gamma_init
    while error > 1e-5:
        gamma = gamma_new
        M_5 = M_ConvDivNozzle(gamma, p5_tot, p5)
        T_5 = Static_Temp(T_5_tot, gamma, M_5)
        Cp = findCp((T_4_tot + T_5)/2, f) # Heat capacity is computed over the expansion
        gamma_new  = Cp/ (Cp - R)
        error = np.abs((gamma_new - gamma)/gamma)
    return gamma, T_5, M_5

def IterateGammaFromT_ConvNozzle(gamma_init, T_4_tot, T_5_tot, m_dot_a, m_dot_f):
    """For CHOKED NOZZLES
    If m_dot_a is not known and f is zero (secondary flow), set m_dot_a to 1
    """
    R = 287.058 
    f = m_dot_f / m_dot_a
    error = 1
    gamma_new = gamma_init
    while error > 1e-5:
        gamma = gamma_new
        T_5 = Static_Temp(T_5_tot, gamma, 1)
        Cp = findCp((T_4_tot + T_5)/2, f) # Heat capacity is computed over the expansion
        gamma_new  = Cp/ (Cp - R)
        error = np.abs((gamma_new - gamma)/gamma)
    return gamma, T_5

def IterateGammaFromP_ConvNozzle(gamma_n_init, p_8, p_8_tot, T_8_tot, m_dot_a, m_dot_f):
    """For ADAPTED NOZZLES"""
    R = 287.058 
    f = m_dot_f / m_dot_a
    error = 1
    gamma_n_new = gamma_n_init
    while error > 1e-5:
        gamma_n = gamma_n_new
        M8 = MachAdaptedNozzle(p_8_tot, p_8, gamma_n)
        T_8 = Static_Temp(T_8_tot, gamma_n, M8)
        Cp = findCp(T_8_tot , f) # Heat capacity is not computed over the expansion
        gamma_n_new  = Cp/ (Cp - R)
        error = np.abs((gamma_n_new - gamma_n)/gamma_n)
    return M8, gamma_n, T_8

def ConvergingNozzles(p_5_tot, T_5_tot, m_dot_ap, m_dot_f, p_0, v_0):
    """Assumptions: no pressure losses, adiabatic nozzle
        m_dot_ap : air, primary flow
    """
    p_8_tot = p_5_tot 
    T_8_tot = T_5_tot 

    R = 287.058 
    f = m_dot_f / m_dot_ap
    error = 1

    # First, we need to evaluate if the nozzle is choked or not 
    
    gamma8_star, T_8_star = IterateGammaFromT_ConvNozzle(1.4, T_5_tot, T_8_tot, m_dot_ap, m_dot_f)
    NPR_star = NozzlePressureRatio(gamma8_star)

    if (p_8_tot/p_0) >=  NPR_star: 
        # The nozzle is choked and the values computed for NPR* are still valid
        print(f'Nozzle is choked') 
        M8 = 1 
        gamma8 = gamma8_star
        T_8 = T_8_star
        a8 = SpeedSound(gamma8, T_8)
        p8 = StaticPressure(p_8_tot, gamma8, 1)
        v8 = a8 * M8
        A8 = ExhaustArea((m_dot_ap + m_dot_f), v8, T_8, p8) 
        T = (m_dot_ap + m_dot_f ) * v8 - m_dot_ap * v_0 + ( p8 - p_0 ) * A8 
    
    else:
        print(f'Nozzle is adapted')
        # The nozzle is adapted and the values computed for NPR* are NO LONGER valid
        p_8 = p_0
        M8, gamma8, T_8 = IterateGammaFromP_ConvNozzle(1.4, p_8, p_8_tot, T_8_tot, m_dot_ap, m_dot_f)
        print(f'M8, T8', M8, T_8)
        a8 = SpeedSound(gamma8, T_8)
        p8 = StaticPressure(p_8_tot, gamma8, M8)
        v8 = a8 * M8
        T = (m_dot_ap + m_dot_f ) * v8 - m_dot_ap * v_0
        A8 = ExhaustArea((m_dot_ap + m_dot_f), v8, T_8, p8)  # No need to compute it as the pressure difference is 0 
    return T, v8, A8, p8, gamma8, T_8

''' Operating lines '''

def ChokedMassFlow(p_tot, T_tot, gamma, A):
    """
    ṁ* = A · p_tot · sqrt(gamma/(R·T_tot)) · (2/(gamma+1))^((gamma+1)/(2·(gamma-1)))
    """
    R = 287.05
    exponent = (gamma + 1) / (2 * (gamma - 1))
    return A * p_tot * np.sqrt(gamma / (R * T_tot)) * (2 / (gamma + 1))**exponent

def IterateMatchingFan_and_SecondaryNozzle(pi_fs_values, eta_fs_values, m_dot_star_values, p_1_tot, T_1_tot, alpha, A_nozzle, p_a):
    """We assume the nozzle remains choked (as in primary flow), so we can use similarity
       fan_pi and eta_s_fan are tables used to interpolate wrt to m_dot. Select 4 points 
       from operating map to put in tables"""
    p_ref = 101325.0
    T_ref = 288.15
    tol_rel = 1e-4
    max_iter = 100
    NPR_crit = 1.89

    f_eta = interp1d(m_dot_star_values, eta_fs_values, kind='linear', bounds_error=False, fill_value='extrapolate')
    f_pi  = interp1d(m_dot_star_values, pi_fs_values,  kind='linear', bounds_error=False, fill_value='extrapolate')

    def secondary_mass_flow(m_dot_star):
        m_dot_a = Uncorrected_m_dot(m_dot_star, p_1_tot, T_1_tot)
        __, m_dot_s = Primary_Secondary_M_dot(alpha, m_dot_a)
        return m_dot_s

    def nozzle_mass_flow(m_dot_star):
        pi_fs  = float(f_pi(m_dot_star))
        eta_fs = float(f_eta(m_dot_star))
        p_2_tot = pi_fs * p_1_tot
        NPR = p_2_tot / p_a
        assert NPR > NPR_crit, \
            f"Nozzle not choked at m_dot_star={m_dot_star:.1f}, NPR={NPR:.3f} < NPR_crit={NPR_crit:.3f}"
        gamma_12s, Cp_12s, T_2_tot = IterateGamma_Comp(1.4, T_1_tot, pi_fs, eta_fs, 0)
        return ChokedMassFlow(p_2_tot, T_2_tot, gamma_12s, A_nozzle)

    def delta_m(m_dot_star):
        return secondary_mass_flow(m_dot_star) - nozzle_mass_flow(m_dot_star)

    m_star_1 = float(np.min(m_dot_star_values))
    m_star_2 = float(np.max(m_dot_star_values))
    delta_1  = delta_m(m_star_1)
    delta_2  = delta_m(m_star_2)

    if delta_1 * delta_2 > 0:
        raise ValueError(  f"Root not bracketed: Δm1={delta_1:.4f}, Δm2={delta_2:.4f} "  f"at m_dot_star=[{m_star_1:.1f}, {m_star_2:.1f}]")

    m_star_m  = None
    delta_m_m = None
    for i in range(max_iter):
        m_star_m  = (m_star_1 + m_star_2) / 2.0
        delta_m_m = delta_m(m_star_m)

        m_dot_s_m = secondary_mass_flow(m_star_m)
        if abs(delta_m_m) < tol_rel * m_dot_s_m:
            print(f"  Converged in {i+1} iterations  |Δm|={abs(delta_m_m):.2e} kg/s")
            break

        if delta_m_m * delta_1 > 0:
            m_star_1 = m_star_m
            delta_1  = delta_m_m
        else:
            m_star_2 = m_star_m
    else:
        print(f"  No convergence after {max_iter} iterations  |Δm|={abs(delta_m_m):.2e} kg/s")

    pi_fs_conv  = float(f_pi(m_star_m))
    eta_fs_conv = float(f_eta(m_star_m))
    m_dot_real  = m_star_m * (p_1_tot / p_ref) / np.sqrt(T_1_tot / T_ref)
    __, m_dot_s_conv = Primary_Secondary_M_dot(alpha, m_dot_real)
    p_2_tot = pi_fs_conv * p_1_tot
    _, _, T_2_tot = IterateGamma_Comp(1.4, T_1_tot, pi_fs_conv, eta_fs_conv, 0)

    return m_dot_real, m_dot_s_conv, pi_fs_conv, eta_fs_conv, p_2_tot, T_2_tot

def IterateToMatchThrust_OperatingPoint(pi_fs_values, eta_fs_values, m_dot_star_values, p_1_tot, T_1_tot, p_a, v_a, Thrust_ref):
    """Uses the secant method to find the operating point based on a given thrust"""
    p_r = 101325.0
    T_r = 288.15
    tol_rel = 1e-4
    max_iter = 100
    NPR_crit = 1.89

    f_eta = interp1d(m_dot_star_values, eta_fs_values, kind='linear', bounds_error=False, fill_value='extrapolate')
    f_pi  = interp1d(m_dot_star_values, pi_fs_values,  kind='linear', bounds_error=False, fill_value='extrapolate')

    def actual_m_dot(m_dot_star):
        m_dot = Uncorrected_m_dot(m_dot_star, p_1_tot, T_1_tot)
        return m_dot
    
    def nozzle_thrust(m_dot_star):
        pi_fs  = float(f_pi(m_dot_star))
        eta_fs = float(f_eta(m_dot_star))
        p_2_tot = pi_fs * p_1_tot
        gamma2, Cp2, T_2_tot = IterateGamma_Comp(1.4, T_1_tot, pi_fs, eta_fs, 0)
        T, __, __, __, __, __ = ConvergingNozzles(p_2_tot, T_2_tot, actual_m_dot(m_dot_star), 0, p_a, v_a)
        return T
    
    def delta_T(m_dot_star):
        return nozzle_thrust(m_dot_star) - Thrust_ref
    
    # Initatialisation using the extrema from operating map
    m0 = float(np.min(m_dot_star_values))
    m1 = float(np.max(m_dot_star_values))
    T0 = nozzle_thrust(m0)
    T1 = nozzle_thrust(m1)
    print(f"  Init: m*={m0:.4f} T={T0:.4f} N  |  m*={m1:.4f} T={T1:.4f} N")

    m_prev, T_prev = m0, T0
    m_curr, T_curr = m1, T1

    for i in range(1, max_iter + 1):
        dT = T_curr - T_prev
        if abs(dT) < 1e-12:
            print(" dT too small")
            break

        # Linear interpolation to thrust_ref
        m_next = m_curr - (T_curr - Thrust_ref) * (m_curr - m_prev) / dT
        T_next = nozzle_thrust(m_next)
        print(f"  Iter {i}: m*={m_next:.5f}  pi={float(f_pi(m_next)):.5f}"
              f"  eta={float(f_eta(m_next)):.6f}  T={T_next:.4f} N")

        if abs(T_next - Thrust_ref) < tol_rel * Thrust_ref:
            print(f"  Converged in {i} iterations")
            m_prev, T_prev = m_curr, T_curr
            m_curr, T_curr = m_next, T_next
            break

        m_prev, T_prev = m_curr, T_curr
        m_curr, T_curr = m_next, T_next

    # Values at convergence
    gamma = 1.4
    pi_conv  = float(f_pi(m_curr))
    eta_conv = float(f_eta(m_curr))
    p_2_tot  = pi_conv * p_1_tot
    T_2_tot  = T_1_tot * (1 + (pi_conv**((gamma - 1) / gamma) - 1) / eta_conv)
    m_dot    = actual_m_dot(m_curr)
    T, v8, A8, p8, gamma8, T8 = ConvergingNozzles(p_2_tot, T_2_tot, m_dot, 0, p_a, v_a)

    return T, A8, m_dot

''' Space propulsion '''
# ATTENTION, ON DOIT CHANGER R !!!!!

def MachFromPressureRatio(pi_e, gamma):
    exp = (gamma-1)/gamma
    M = np.sqrt(2 * (pi_e**exp - 1) / (gamma - 1))
    return M

def SpeedSoundSpace (gamma, T, R):
    return (gamma * R * T)**0.5

def SpecificImpulse(T, m_dot):
    """Uses g = 9.81"""
    g = 9.81 # m/s**2
    return T/(m_dot*g)

def Mdot_from_MR(MR, m_dot_tot):
    """Mass flow rates for oxidiser and fuel BEFORE combustion"""
    m_dot_O2 = MR * m_dot_tot / (1 + MR)
    m_dot_H2 = m_dot_tot / (1 + MR) 
    return m_dot_O2, m_dot_H2

def Mdot_from_Isp(Isp, T):
    """Uses g = 9.81"""
    g = 9.81
    return T /(Isp*g)

def MInit_res(m_tot_comsumed, m_remaining_in_percent):
    return  m_tot_comsumed / (1 - m_remaining_in_percent)

def VolumeRes(m_init, R, T_r_tot, p_r_tot_init):
    return m_init * R * T_r_tot / p_r_tot_init

def At_from_CT(C_T, p_feed_tot, T):
    return T/(C_T * p_feed_tot)

def MolarMass_R_CombustionGas(m_dot_H2O, m_dot_H2, R_star):
    """Returns molar mass on g/mol and R in J/kgK
    If not given, use R_star = 8.314 J/molK"""
    MH2O = 18           # g/mol
    MH2 = 2             # g/mol

    # We compute the molar flow rates
    Qn_H2O = m_dot_H2O / MH2O 
    Qn_H2 = m_dot_H2 / MH2
    Qn = Qn_H2O + Qn_H2

    M = Qn_H2O / Qn * MH2O + Qn_H2 / Qn * MH2
    R = R_star / M * 1000
    return M, R

def CharacterisicVelocity_TComb(T_comb_tot, R, gamma):
    C = np.sqrt(gamma*R*T_comb_tot) * ((gamma+1)/2)**((gamma+1)/(2*(gamma-1))) / gamma
    return C

def MachFct(M, gamma):
    exp = - (gamma + 1) / (2 * (gamma-1))
    fct = np.sqrt(gamma) * M * (1 + (gamma-1) * M**2 / 2)**exp
    return fct

def AreaRatio(Me, gamma):
    return MachFct(1, gamma)/MachFct(Me, gamma)

def CT_from_Me_AreaRatio(gamma, Me, A_ratio, pi_e):
    """ATTENTION: considers pe = pa"""
    return gamma * A_ratio * Me**2 / pi_e

def FindPressureRatioAtWhichSkirtExpands(A_ratio_1, A_ratio_2, gamma, M1, M2, pi1, pi2, K):
    """Found by equating both thrust coeffs at current pressure ratio, assuming a full isentropic expansion.
       If one nozzle is overexanded, verify that no separation occurs using Sommerfeld criterion (pa/p_tot * p_tot/p_out_nozzle < K)"""
    deno = A_ratio_1 * (1+gamma*M1**2)/pi1 - A_ratio_2 * (1+gamma*M2**2)/pi2
    pi = (A_ratio_1 - A_ratio_2) / deno
    if pi < pi1: 
        print(f'Nozzle 1 overexpanded')
        if pi / pi1 > K:
            print(f'separation')
        else:
            print(f'no separation')
    if pi > pi1: 
        print(f'Nozzle 1 underexpanded')
    if pi < pi2: 
        print(f'Nozzle 2 overexpanded')
        if pi / pi2 > K:
            print(f'separation')
        else:
            print(f'no separation')
    if pi > pi2: 
        print(f'Nozzle 2 underexpanded')
    return pi

def IterateMeTotalCds_FindThrust(gamma, R, p_tot, T_tot, A_ratio, m_dot, p_a):
    """Implements bissection method to finc M_e, the Mach nb at the exhaust since F(M) is always increasing 
       and the throat is choked. A_ratio = Ae/At.

       Set m_dot to 0 if it is not known
       
       Returns T, Isp, C_T (thrust coeff) and pe_star
       """
    def MachFct(M):
        exp = - (gamma + 1) / (2 * (gamma-1))
        fct = np.sqrt(gamma) * M * (1 + (gamma-1) * M**2 / 2)**exp
        return fct
    
    def res(Me):
        res = MachFct(Me) - MachFct(1) / A_ratio
        return res
    
    # Initialisation of bissection method (we know Me>1)
    M_low, M_high = 1.001, 50.0          
    tol           = 1e-8
    max_iter      = 200

    print(f"res(M_low)  = {res(M_low)}")
    print(f"res(M_high) = {res(M_high)}")
    
    assert res(M_low) > 0 and res(M_high) < 0, "Root out of the interval"

    for i in range(max_iter):
        M_mid = 0.5 * (M_low + M_high)
        if res(M_mid) > 0:      # We are left of the root => increase M_low
            M_low = M_mid
        else:                    # We are right of the root => decrease M_high
            M_high = M_mid
        if (M_high - M_low) < tol:
            break

    Me = 0.5 * (M_low + M_high)
    print(f" Me = {Me:.6f}  (converged in {i+1} iterations)")
    pe = StaticPressure(p_tot, gamma, Me)
    Te = Static_Temp(T_tot, gamma, Me)
    ae = SpeedSoundSpace(gamma, Te, R)
    ve = Vel_from_Mach(ae, Me)
    print(f" pe {pe}, Te {Te}, ve{ve}")

    Ae = m_dot * np.sqrt(R*T_tot) / (p_tot * MachFct(Me))
    At = m_dot * np.sqrt(R*T_tot) / (p_tot * MachFct(1))
    print(f" At {At}, Ae {Ae}")

    T = m_dot * ve  + ( pe - p_a ) * Ae 
    Isp = SpecificImpulse(T, m_dot)
    C_T = T / (p_tot * At)

    # if m_dot is not known, we determine C_star and C_T
    if m_dot == 0:
        g = 9.81
        C_T = A_ratio * (1 + gamma * Me**2) * pe / p_tot
        C_star = np.sqrt(R*T_tot)/ MachFct(1)
        Isp = C_T * C_star / g
        
        return Isp, C_T

    else : 
        return T, Isp, C_T, pe

def FindThrustSeparationInNozzle(gamma, R, p_tot, T_tot, m_dot, p_a):
    """
    Sommerfeld criterion: separation at p_sep = p_a / K.
    The effective nozzle exit is the cross-section where p = p_sep.
    Downstream of separation, pressure = p_a.
 
    Returns: T [N], Isp [s], C_T [-], M_sep [-], v_sep [m/s], A_sep [m²]
    """
    def MachFct(M):
        exp = - (gamma + 1) / (2 * (gamma-1))
        fct = np.sqrt(gamma) * M * (1 + (gamma-1) * M**2 / 2)**exp
        return fct

    K = 2.5 
    p_sep = p_a / K
 
    # Find Mach number where static pressure equals p_sep
    def residual_sep(M):
        return StaticPressure(p_tot, gamma, M) - p_sep

    M_lo, M_hi = 1.001, 60.0
    # StaticPressure is decreasing: residual > 0 for small M, < 0 for high M 
    for _ in range(300):
        M_mid = 0.5 * (M_lo + M_hi)
        if residual_sep(M_mid) > 0:
            M_lo = M_mid
        else:
            M_hi = M_mid
        if (M_hi - M_lo) < 1e-10:
            break

    M_sep = 0.5 * (M_lo + M_hi)

    T_sep = Static_Temp(T_tot, gamma, M_sep)
    a_sep = SpeedSoundSpace(gamma, T_sep, R)
    v_sep = Vel_from_Mach(a_sep, M_sep)
 
    # Area at separation from mass-flow function
    A_sep = m_dot * np.sqrt(R*T_tot) / (p_tot * MachFct(M_sep))
    At    = m_dot * np.sqrt(R*T_tot) / (p_tot * MachFct(1))

    T   = m_dot * v_sep + (p_sep - p_a) * A_sep
    Isp = SpecificImpulse(T, m_dot)
    C_T = T / (p_tot * At)
    return T, Isp, C_T, M_sep, v_sep, A_sep