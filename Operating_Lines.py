import numpy as np 
from findCp import findCp
from functions import *
from Convert_unit_and_ISA import *

def ex711():

    '''Data'''
    pi_fp = 1.49
    pi_fs = 1.6
    m_dot_star = 415
    alpha = 5
    
    '''1. Compute m_dot, NPR_s, contribution of sec. flow to thrust, A_s, work absorbed by complete fan'''
    Mc = 0.82
    alt = length_ft(35000)
    T_a, a_a, p_a, rho_a = isa(alt)
    print(f'Ta, pa', T_a, p_a)
    v_a = Vel_from_Mach(a_a, Mc)

    eta_fs = 0.84 # Based on operating map

    T_1_tot, gamma_0 = TotalTemp_gamma(T_a, 0 , Mc)
    tau_r = T_1_tot/T_a
    p_1_tot = tau_r**(gamma_0/(gamma_0-1))*p_a
    print(f'T1tot, p1tot', T_1_tot, p_1_tot)
    m_dot = Uncorrected_m_dot(m_dot_star, p_1_tot, T_1_tot)
    print(f'm_dot', m_dot)

    NPR = pi_fs * p_1_tot / p_a

    m_dot_p, m_dot_s = Primary_Secondary_M_dot(alpha, m_dot)
    print(f'm_dotp, m_dots', m_dot_p, m_dot_s)
    gamma12s, Cp12s, T_2s_tot = IterateGamma_Comp(1.4, T_1_tot, pi_fs, eta_fs, 0)
    p_2s_tot = pi_fs*p_1_tot
    p_5s_tot = p_2s_tot
    T_5s_tot = T_2s_tot
    Ts, __, A_s, __, __, __= ConvergingNozzles(p_5s_tot, T_5s_tot, m_dot_s, 0, p_a, v_a)
    print(f'Ts (kN), As, NPR', Ts/1000, A_s, NPR)

    P_fs = PowerComp(m_dot_s, Cp12s, T_1_tot, T_2s_tot) # Work absorbed by secondary flow
    print(f'P_fs (MW)', P_fs/1e6)

    # Word absorbed by primary flow can be found from operating point on map 
    eta_fp = 0.85
    gamma12p, Cp12p, T_2p_tot = IterateGamma_Comp(1.4, T_1_tot, pi_fp, eta_fp, 0)
    P_fp = PowerComp(m_dot_p, Cp12p, T_1_tot, T_2p_tot) # Work absorbed by secondary flow
    P_tot = P_fs + P_fp
    print(f'P_fp (MW)', P_fp/1e6)
    print(f'P_tot (MW)', P_tot/1e6)

    '''2. For same altitude and flight speed, assuming alpha is unchanged, determine fan operation point 
          for partir rotation speeds 80% and 90% and in particular m_dot_s, Pi_s, secondary flow contribution
          to net thrust and work absorbed by complete fan'''
    print('Part 2: Iteration')
    # We need to find the operating point of the fan iteratively by identifying the m_dot on the fan operating
    # line for which the total cds upstream of the nozzle result in the same nozzle m_dot

    # Assuming the nozzle remains choked, we can use similarity.
    # We should select 4 points from operation map
    m_dot_star_80_values = [275, 300, 325, 339]
    eta_fs_80_values = [0.746, 0.793, 0.839, 0.82]
    pi_fs_80_values = [1.34, 1.37, 1.38, 1.25]

    m_dot_fan, m_dot_s_conv, pi_fs_conv,  eta_fs_conv,  p_2_tot,  T_2_tot = IterateMatchingFan_and_SecondaryNozzle(pi_fs_80_values, eta_fs_80_values, m_dot_star_80_values, p_1_tot, T_1_tot, alpha, A_s, p_a)
    print(m_dot_fan, m_dot_s_conv)
    
def ex712():
    M_a = 0.85
    alt  = length_ft(10000)
    alpha = 20
    p_r = 101325
    T_r = 288
    r_h = 55/2 / 100 # m
    r_t = 80/2 / 100 # m
    m_dot_star_100_values = [74.04,    75.55,   78.27,   79.59,   81.10, 81.86]
    pi_fs_100_values  = [1.3054,  1.3040,  1.2932,  1.2837,  1.2672, 1.2552]
    eta_fs_100_values = [0.9016, 0.91377, 0.92997, 0.93454, 0.93422, 0.92911]

    A_nozzle = (np.pi * r_t**2) - (np.pi * r_h**2)

    '''Using interpolated caracteristics based on values tables, calculate secondary flow contribution to thrust, 
       full power consumption of fan, find both actual and corrected mass flow rates. Is the nozzle choked?'''
    
    # We assume the fan is rotating at its design speed in CORRECTED CDS, such that we can use the characteristic lines. 
    # The challenge is to find the intersection btw operating line of comp and converging nozzles such that m_dot and pi correspond.

    #1. Fan upstream cds
    T_a, a_a, p_a, rho_a = isa(alt)
    print(f'Ta, pa', T_a, p_a)
    v_a = Vel_from_Mach(a_a, M_a)
    T_1 = T_a
    p_1 = p_a
    T_2_tot, gamma = TotalTemp_gamma(T_1, 0, M_a)
    p_2_tot = TotalPressure(p_1, gamma, M_a)

    m_dot_fan, m_dot_s_conv, pi_fs_conv,  eta_fs_conv,  p_3_tot,  T_3_tot = IterateMatchingFan_and_SecondaryNozzle(pi_fs_100_values, eta_fs_100_values, m_dot_star_100_values, p_2_tot, T_2_tot, alpha, A_nozzle, p_a)
    print(m_dot_fan, m_dot_s_conv)

    T, __, __,__, __, __ = ConvergingNozzles(p_3_tot, T_3_tot, m_dot_s_conv, 0, p_a, v_a)
    print(f'Thrust(kN)', T/1000)
    Cp23= findCp((T_3_tot+T_2_tot)/2, 0)
    P_fan = PowerComp(m_dot_fan, Cp23, T_2_tot, T_3_tot)
    print(f'Power fan(MW)', P_fan/1e6)

def ex713():
    """Determine ratio of area of converging nozzles wrt fan frontal aera An/A1 which maximizes 
       propulsive power for a given power on the shaft. Provide an estimate for the optimal 
       nozzle to inlet area ratio An/A1 and the corresponding fan operating point on the characteristic
       in terms of m_dot_star, pi, eta_is. Provide ACTUAL in-flight values for T and m_dot per frontal area
       as well as overall efficiency
   """
    pi_fan_star = 1.54
    eta_f_is_star = 0.92
    alt = length_ft(30000)
    Ma = 0.85
    Ga_star = 192

    G_aOverG_a_star_points = [ 0.894791 ,0.959345 ,1.007327 ,1.032927, 1.050198]
    PiOverPi_star_points = [1.035155, 1.028333, 1.018101, 1.00062, 0.969922]

    PiOverPi_star_points_2 = [ 0.971654, 0.992714 ,1.005178, 1.017584, 1.02435, 1.030144,1.034616]
    EtaOverEtastar_points_2 = [0.905892, 1, 1.008299, 1.002905, 0.981743 ,0.93361 ,0.877178]

    T_a, a_a, p_a, rho_a = isa(alt)
    print(f'Ta, pa', T_a, p_a)
    v_a = Vel_from_Mach(a_a, Ma)

    T_1_tot, gamma_0 = TotalTemp_gamma(T_a, 0 , Ma)
    p_1_tot = TotalPressure(p_a, gamma_0, Ma)
    print(f'T1tot, p1tot', T_1_tot, p_1_tot)

    # We sweep over the table  G_aOverG_a_star_points and compute the propulsive power for ecah possible case
    x = np.linspace(np.min(G_aOverG_a_star_points), np.max(G_aOverG_a_star_points), 100) 

    f_pi = interp1d(G_aOverG_a_star_points, PiOverPi_star_points, kind='linear', bounds_error=False, fill_value='extrapolate')
    f_eta  = interp1d(PiOverPi_star_points_2, EtaOverEtastar_points_2,  kind='linear', bounds_error=False, fill_value='extrapolate')
    
    # Initialisation of empty tables 
    eta_prop = np.zeros(len(x))
    eta_th = np.zeros(len(x))
    A_ratio = np.zeros(len(x))
    A1 = np.zeros(len(x))
    T = np.zeros(len(x))

    for i in range (len(x)):
        pi_fs  = float(f_pi(x[i]))
        if pi_fs > np.max(PiOverPi_star_points_2) or pi_fs < np.min(PiOverPi_star_points_2) :
            continue
        eta_fs = float(f_eta(pi_fs))
        pi_fs_real = pi_fs * pi_fan_star 
        gamma2, Cp2, T_2_tot = IterateGamma_Comp(1.4, T_1_tot, pi_fs_real, eta_fs, 0)
        p_2_tot = pi_fs_real * p_1_tot

        m_dot = Uncorrected_m_dot(x[i] * Ga_star, p_1_tot, T_1_tot)
        T[i], v3, A3, p3, gamma3, T_3 = ConvergingNozzles(p_2_tot, T_2_tot, m_dot, 0, p_a, v_a)
        P_fan = PowerComp(m_dot, Cp2, T_1_tot, T_2_tot)
        print(rho_a)
        eta_prop [i] = PropulsiveEff(m_dot, 0, v3, v_a, T[i])
        eta_th [i] = ThermalEff_from_PowerConsumed(m_dot, v3, v_a, P_fan)
        A1 [i] = m_dot / (rho_a * v_a) 
        A_ratio[i] = A3/A1[i]

    index = np.argmax(eta_prop)
    A1_opt = A1[index]
    eta_overall = OverallEff(eta_th[index], eta_prop[index])
    m_star_opt = x[index] * Ga_star * A1_opt         
    pi_opt = f_pi(x[index]) * pi_fan_star 
    eta_is_opt = f_eta(f_pi(x[index])) * eta_f_is_star  


    print(f'Area ratio', A_ratio[index])
    print(f'Fan operating point:  m_dot_star: {m_star_opt},  pi_f_star: {pi_opt}, eta_s_star: {eta_is_opt}')
    print(f'In-flight values: T/A1: {T[index]/A1_opt},  m/A1: {x[index] * Ga_star}, Overall eff: {eta_overall}')

def ex714():
    '''Secondary outlet duct is converging and adapts itself to guarantee ideal expansion at all operating points'''
    Ma = 0.5
    alt = length_ft(30000)
    n = 37780                   # rpm
    Thrust_ref = 71             # N

    """1. Determine current operating point in terms of corrected n and m_dot, pi and eff
       2. Compute actual m_dot, total cds after fan, static cds in nozzle and exhaust velocity
       3. Determine nozzle outlet area
    """

    T_a, a_a, p_a, rho_a = isa(alt)
    print(f'Ta, pa', T_a, p_a)
    v_a = Vel_from_Mach(a_a, Ma)
    T_1 = T_a
    p_1 = p_a
    T_1_tot, gamma1 = TotalTemp_gamma(T_1, 0, Ma)
    p_1_tot = TotalPressure(p_1, gamma1, Ma)
    print(f'T1tot, p1tot', T_1_tot, p_1_tot)
    n_star = nStar(n, T_1_tot)
    print(f'nstar', n_star)
    # We find the fan operating point on the corresponding characteristic. Curve represented by cubic interpolation points
    m_dot_star_values = [2.5, 2.7, 2.8, 2.9, 3]
    pi_f_values = [1.19, 1.19, 1.18, 1.175, 1.15]
    eta_s_values = [0.62, 0.75, 0.82, 0.83, 0.75]
    
    T, A_conv, m_dot = IterateToMatchThrust_OperatingPoint(pi_f_values, eta_s_values, m_dot_star_values, p_1_tot, T_1_tot, p_a, v_a, Thrust_ref)
   
    print(f'T', T)
    print(f'Aconv', A_conv)
    print(f'm_dot', m_dot)


#ex711()
#ex712()
ex713()
#ex714()