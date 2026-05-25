import numpy as np
from .findCp import findCp
from .atmosphere_isa import Air
from scipy.interpolate import interp1d

def av(x: float, y: float):
    return (x+y)/2

def get_gamma(Cp):
    R = 287.058
    return Cp/(Cp-R)


def find_Qc(hf1, hf2, mdotf, eta_cc):
    return (hf2-hf1)*eta_cc*mdotf


def Heat_Rate(hf1, hf2, mdotf, Pe, eta_cc=1):
    return find_Qc(hf1, hf2, mdotf, eta_cc) / Pe

def adiabatic_convdiv_nozzle(P1,m_dot,f,pa,gamma13 = 1.37,tol = 1e-5):
    """
    Computes flow through a convergent nozzle with choking condition.

    Args:
        P1 (Air): Inlet state.
        m_dot (float): Mass flow rate (kg/s).
        f (float): Fuel-air ratio.
        pa (float): Ambient pressure (Pa).
        gamma13 (float, optional): Heat capacity ratio. Defaults to 1.37.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (Outlet state, choked flag, pressure ratio, area)
    """
    P3 = Air()
    P3.p0 = P1.p0
    P3.T0 = P1.T0
    NPR = P3.p0/pa
    P3.cp = findCp(P3.T0,f)
    P3.get_gamma_from_cp()

    diff = np.inf

    while diff > tol:        
        P3.M = np.sqrt(2/(gamma13 -1)*((P3.p0/pa)**((gamma13-1)/gamma13)-1))
        P3.T = P3.T0/((1+(gamma13-1)/2*P3.M**2))
        Cp13 = findCp(av(P3.T,P1.T0),f)
        new_gamma = get_gamma(Cp13)
        diff = abs(new_gamma - gamma13)/gamma13
        gamma13 = new_gamma
        P3.gamma = gamma13

    P3.p = P3.p0/((1+(gamma13-1)/2*P3.M**2))**(gamma13/(gamma13-1))
    P3.a = np.sqrt(P3.gamma*P3.R*P3.T)
    P3.v = P3.M*P3.a
    P3.rho0 = P3.p0 / (P3.R * P3.T0)
    P3.h0 = P3.cp * P3.T0
    P3.s = P1.s + P3.cp*np.log(P3.T0/P1.T0) - P3.R*np.log(P3.p0/P1.p0)
    P3.rho = P3.p/P3.R/P3.T
    P3.set_static()
    At = m_dot*(1+f)/P3.rho/P3.v

    return P3, NPR, At

def adiabatic_convergent_nozzle(P1,m_dot,f,pa,gamma13 = 1.37,tol = 1e-5):
    """
    Computes flow through a convergent nozzle with choking condition.

    Args:
        P1 (Air): Inlet state.
        m_dot (float): Mass flow rate (kg/s).
        f (float): Fuel-air ratio.
        pa (float): Ambient pressure (Pa).
        gamma13 (float, optional): Heat capacity ratio. Defaults to 1.37.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (Outlet state, choked flag, pressure ratio, area)
    """
    P3 = Air()
    P3.p0 = P1.p0
    P3.T0 = P1.T0
    NPR = P3.p0/pa
    P3.cp = findCp(P3.T0,f)
    P3.get_gamma_from_cp()
    NPR_star = (1 + (P3.gamma-1)/2)**(P3.gamma/(P3.gamma-1))
    diff = np.inf

    while diff > tol:

        if NPR < NPR_star :
            choked = False
            P3.M = np.sqrt(2/(gamma13 -1)*((P3.p0/pa)**((gamma13-1)/gamma13)-1))
            P3.set_static()
            Cp13 = findCp(av(P3.T,P1.T0),f)
            new_gamma = get_gamma(Cp13)
            diff = abs(new_gamma - gamma13)/gamma13
            gamma13 = new_gamma

        else:
            choked = True
            P3.p = P3.p0 / (1 + (gamma13-1)/2)**(gamma13/(gamma13-1))
            P3.T = P3.T0 / (1 + (gamma13-1)/2)
            P3.M = 1

            Cp13 = findCp(av(P3.T,P1.T0),f)
            new_gamma = get_gamma(Cp13)
            diff = abs(new_gamma - gamma13)/gamma13
            gamma13 = new_gamma

    P3.gamma = gamma13
    P3.a = np.sqrt(P3.gamma*P3.R*P3.T)
    P3.v = P3.M*P3.a
    P3.rho0 = P3.p0 / (P3.R * P3.T0)
    P3.h0 = P3.cp * P3.T0
    P3.s = P1.s + P3.cp*np.log(P3.T0/P1.T0) - P3.R*np.log(P3.p0/P1.p0)
    P3.rho = P3.p/P3.R/P3.T
    A = m_dot/P3.rho/P3.v

    return P3, choked, NPR, A

def iterate_M_nozzle(P1,pa,f,gamma = 1.39,tol = 1e-4):
    """Compute iterate M nozzle."""
    P2 = Air(R = P1.R)
    P2.T0 = P1.T0
    P2.p0 = P1.p0
    diff = np.inf
    while diff > tol : 
        P2.M = np.sqrt(2/(gamma-1)*((P1.T0/pa)**((gamma-1)/gamma)-1))
        P2.T = P1.T0/((1+(gamma-1)/2*P1.M**2))
        cp = findCp(av(P2.T,P1.T0), f)
        gamma_new = get_gamma(cp)
        error  = abs(gamma_new-gamma)/gamma
        gamma = gamma_new

    P2.gamma = gamma
    P2.a = np.sqrt(P2.gamma*P2.R*P2.T)
    P2.v = P2.M*P2.a
    
    return P2,cp


def turbofan_thrust(m_dot_p,m_dot_s,f,P0,PP,PS= None,AP= 0,AS = 0):
    """
    Computes turbofan thrust (primary and secondary streams).
    """
    Tp = m_dot_p * (1+f) * PP.v - m_dot_p * P0.v + (PP.p-P0.p)*AP
    if PS is not None:
        Ts = m_dot_s * PS.v - m_dot_s*P0.v + (PS.p-P0.p)*AS
    else:
        Ts = 0
    return Tp,Ts


def simple_thrust(m_dot,f,P0,P,A= 0):
    """
    Computes simple jet thrust.
    """
    Tp = m_dot * (1+f) * P.v - m_dot * P0.v + (P.p-P0.p)*A
    return Tp


def thermal_efficiency(vin,vout,mdot,f,delta_hf):
    """
    Computes thermal efficiency.
    """
    return (mdot*(1+f)*vout**2 - mdot * vin**2)/(2*mdot*f*delta_hf)


def propulsive_efficiency(vin,vout,mdot,f,T):
    """ 
    Computes propulsive efficiency.
    """
    return 2*T*vin/(mdot*(1+f)*vout**2 - mdot * vin**2)

def thermal_efficiency_turbofan(m_dot_p, m_dot_f, m_dot_s, v_eff_p, v_eff_s, m_dot_a, v_a, delta_hf):
    """Compute thermal efficiency turbofan."""
    return (((m_dot_p + m_dot_f)*v_eff_p**2 + m_dot_s*v_eff_s**2 - m_dot_a*v_a**2) / (2*m_dot_f*delta_hf))

def propulsive_efficiency_turbofan(Thrust_tot, v_a, m_dot_p, m_dot_f, m_dot_s, v_eff_p, v_eff_s, m_dot_a):
    """Compute propulsive efficiency turbofan."""
    return ((2*Thrust_tot*v_a) / ((m_dot_p + m_dot_f)*v_eff_p**2 + m_dot_s*v_eff_s**2 - m_dot_a*v_a**2))

def corrected_massflow(mdot_ref,P1,Pref):
    """
    Computes corrected mass flow.
    """
    return mdot_ref / (P1.p/Pref.p) / np.sqrt(Pref.T/P1.T)

def true_massflow_from_corrected(mdot_corr,P1,Pref):
    """
    Computes actual mass flow.
    """
    return mdot_corr * (P1.p0/Pref.p0) * np.sqrt(Pref.T0/P1.T0)

def corrected_rpm(rpm_true,P1,Pref):
    """
    Computes corrected RPM.
    """
    return rpm_true / np.sqrt(P1.T0/Pref.T0)


def true_rpm_from_corrected(rpm_ref,P1,Pref):
    """
    Computes corrected RPM.
    """
    return rpm_ref * np.sqrt(P1.T0/Pref.T0)


def corrected_compr(Pc_ref,P_ref,Pin):
    """
    Computes corrected compressor power.
    """
    return Pc_ref* (Pin.p0/P_ref.p0)*np.sqrt(Pin.T0/P_ref.T0)


def mixing(PP,PS,T0_ref,mdotp,mdots,f,pi_m = 1,tol = 1e-5):
    """
    Computes mixing of primary and secondary streams.
    """
    Pm = Air()
    global_f = f*mdotp /(mdotp + mdots)

    Cps = findCp(av(PS.T0,T0_ref),0)
    Cpp = findCp(av(PP.T0,T0_ref),f)

    diff = np.inf
    Pm.T0 = PS.T0

    while diff > tol:
        Cpnozzle = findCp(av(Pm.T0,T0_ref),global_f)
        new_T0 = T0_ref + (mdots*Cps*(PS.T0-T0_ref) + mdotp * (1+f) * Cpp * (PP.T0 -T0_ref))/((mdotp*(1+f)+mdots)*Cpnozzle)
        diff = abs(Pm.T0 - new_T0) / Pm.T0
        Pm.T0 = new_T0

    Pm.p0 = PP.p0*pi_m
    Pm.cp = findCp(Pm.T0,global_f)
    Pm.get_gamma_from_cp()

    return Pm,global_f


def approx_mixing(PP,PS,T0_ref,mdotp,mdots,f,pi_m = 1,tol = 1e-5):
    """
    Computes mixing of primary and secondary streams.
    """
    Pm = Air()
    global_f = f*mdotp /(mdotp + mdots)
    Pm.T0 = (mdotp*PP.T0 + mdots*PS.T0)/(mdotp+mdots)
    Pm.p0 = PP.p0*pi_m
    Pm.cp = findCp(Pm.T0,global_f)
    Pm.get_gamma_from_cp()
    
    return Pm,global_f


def choking_NPR(gamma):
    """Compute choking NPR."""
    return ((gamma + 1 )/2)**(gamma/(gamma-1))


def ChokedMassFlow(P2, gamma, A):
    """
    ṁ* = A · p_tot · sqrt(gamma/(R·T_tot)) · (2/(gamma+1))^((gamma+1)/(2·(gamma-1)))
    """
    exponent = (gamma + 1) / (2 * (gamma - 1))
    return A * P2.p0 * np.sqrt(gamma / (P2.R * P2.T0)) * (2 / (gamma + 1))**exponent

def PropulsiveEff(m_dot_in, m_dot_f, v_j, v_1, T):
    """Computes propulsive efficiency from mass flow rates, exit velocity vj, inlet velocity v_1 and thrust T"""
    eta_p = (2*T*v_1) / ((m_dot_in + m_dot_f) * v_j**2 - m_dot_in * v_1**2 )
    return eta_p


def ThermalEff_from_Power(m_dot_a, v_j, v_1, P):
    """Compute ThermalEff from Power."""
    eta_th = 0.5 * m_dot_a * (v_j**2 - v_1**2) / P
    return eta_th

def overalleff(eta_t,eta_p):
    """Compute overalleff."""
    return eta_p*eta_t

def f(M,gamma):
    """Return the ideal gas flow function f(M, gamma)."""
    return 1 + (gamma-1)/2 * M**2

def F(M,gamma):
    """Return the isentropic flow function F(M, gamma)."""
    return np.sqrt(gamma)* M * f(M,gamma)**(-(gamma + 1)/ (2 * (gamma - 1)))

def set_Mach(P,gamma = 1.4, tol = 1e-4):
    """Return or compute Mach."""
    diff = np.inf
    for _ in range(5) : 
        P.M = P.v / np.sqrt(gamma * P.R * P.T)  #Mach number
        P.set_total()  # sets total conditions with Mach (namely T0)
        P.cp = findCp(P.T0,0) # finds CP from T0
        P.get_gamma_from_cp() #gets the updated gamma 

def interpolate(x,y):
    """Compute interpolate."""
    return interp1d(x,y, kind='linear', bounds_error=False, fill_value='extrapolate')


def effective_velocity(Pe, mdot, p_a, A_e):
    return Pe.v + (Pe.p - p_a) * A_e / mdot

def intercooler(Pp,Ps,eta_intercooler,mdotp,mdots,tol = 1e-5):
    """Compute intercooler."""
    if Pp.cp is None : 
        Pp.cp = findCp(Pp.T0,0)
    if Ps.cp is None : 
        Ps.cp = findCp(Ps.T0,0)
    
    new_Pp = Air()
    new_Ps = Air()
    cpp = Pp.cp
    cps = Ps.cp

    diff = np.inf
    while diff > tol :

        minimum = min(cpp*mdotp,cps*mdots)
        q = eta_intercooler * minimum * (Pp.T0 - Ps.T0)

        new_Pp.T0 = Pp.T0 - q / cpp / mdotp
        new_Ps.T0 = Ps.T0 + q / cps / mdots

        new_cpp = findCp(av(new_Pp.T0,Pp.T0),0)
        cps = findCp(av(new_Ps.T0,Ps.T0),0)

        diff = (new_cpp - cpp)/cpp
        cpp = new_cpp

    new_Pp.p0 = Pp.p0 
    new_Ps.p0 = Ps.p0
    new_Pp.cp = findCp(new_Pp.T0,0)
    new_Ps.cp = findCp(new_Ps.T0,0)
    new_Pp.get_gamma_from_cp()
    new_Ps.get_gamma_from_cp()

    return new_Pp , new_Ps
