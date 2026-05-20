import numpy as np
from .findCp import findCp
from .atmosphere_isa import Air
from scipy.interpolate import interp1d

def av(x: float, y: float):
    """
    Computes arithmetic mean of two values.

    Args:
        x (float): First value.
        y (float): Second value.

    Returns:
        float: Mean of x and y.
    """
    return (x+y)/2

def get_gamma(Cp):
    """
    Computes specific heat ratio from Cp.

    Args:
        Cp (float): Specific heat capacity at constant pressure (J/kg·K).

    Returns:
        float: Specific heat ratio gamma.
    """
    R = 287.058
    return Cp/(Cp-R)


def find_Qc(hf1, hf2, mdotf, eta_cc):
    """
    Computes combustor heat release.

    Args:
        hf1 (float): Fuel inlet enthalpy (J/kg).
        hf2 (float): Fuel outlet enthalpy (J/kg).
        mdotf (float): Fuel mass flow rate (kg/s).
        eta_cc (float): Combustion efficiency.

    Returns:
        float: Heat released (W).
    """
    return (hf2-hf1)*eta_cc*mdotf


def Heat_Rate(hf1, hf2, mdotf, Pe, eta_cc=1):
    """
    Computes heat rate of the engine.

    Args:
        hf1 (float): Fuel inlet enthalpy (J/kg).
        hf2 (float): Fuel outlet enthalpy (J/kg).
        mdotf (float): Fuel mass flow rate (kg/s).
        Pe (float): Engine power output (W).
        eta_cc (float, optional): Combustion efficiency. Defaults to 1.

    Returns:
        float: Heat rate.
    """
    return find_Qc(hf1, hf2, mdotf, eta_cc) / Pe


def isentropic_eff(h1, h2, h2s, type: str = "compr"):
    """
    Computes isentropic efficiency for compressor or turbine.

    Args:
        h1 (float): Inlet enthalpy (J/kg).
        h2 (float): Actual outlet enthalpy (J/kg).
        h2s (float): Isentropic outlet enthalpy (J/kg).
        type (str, optional): Component type ("compr" or "turb"). Defaults to "compr".

    Returns:
        float: Isentropic efficiency.
    """
    if type == "compr":
        return (h2s-h1)/(h2-h1)

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

        if NPR < NPR_star:
            choked = False
            P3.M = np.sqrt(2/(gamma13 -1)*((P3.p0/pa)**((gamma13-1)/gamma13)-1))
            P3.set_static()
            Cp13 = findCp(av(P3.T,P1.T0),f)
            new_gamma = get_gamma(Cp13)
            diff = abs(new_gamma - gamma13)/gamma13
            gamma13 = new_gamma

        else:
            choked = True
            P3.p = P3.p0 / NPR_star
            P3.T = P3.T0 / (1 + (gamma13-1)/2)
            P3.M = 1

            Cp13 = findCp(av(P3.T,P1.T0),f)
            new_gamma = get_gamma(Cp13)
            diff = abs(new_gamma - gamma13)/gamma13
            gamma13 = new_gamma

    P3.a = np.sqrt(P3.gamma*P3.R*P3.T)
    P3.v = P3.M*P3.a
    P3.rho0 = P3.p0 / (P3.R * P3.T0)
    P3.h0 = P3.cp * P3.T0
    P3.s = P1.s + P3.cp*np.log(P3.T0/P1.T0) - P3.R*np.log(P3.p0/P1.p0)
    P3.rho = P3.p/P3.R/P3.T
    A = m_dot/P3.rho/P3.v

    return P3, choked, NPR, A


def adiabatic_convdiv_nozzle(P1,m_dot,f,pout,gamma12 = 1.37,tol = 1e-5):
    """
    Computes flow through a convergent-divergent nozzle.

    Args:
        P1 (Air): Inlet state.
        m_dot (float): Mass flow rate (kg/s).
        f (float): Fuel-air ratio.
        pout (float): Outlet pressure (Pa).
        gamma12 (float, optional): Heat capacity ratio. Defaults to 1.37.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (Outlet state, area)
    """
    P2 = Air()
    P2.p0 = P1.p0
    P2.T0 = P1.T0
    P2.p = pout
    P2.cp = findCp(P2.T0,f)
    P2.get_gamma_from_cp()

    diff = np.inf

    while diff > tol:
        P2.M = np.sqrt(2/(gamma12-1)*((P2.p0/P2.p)**((gamma12-1)/gamma12)-1))
        P2.T = P2.T0/(1 + (gamma12 - 1) / 2 * P2.M**2)
        Cp12 = findCp(av(P2.T,P1.T0),f)
        new_gamma = get_gamma(Cp12)
        diff = abs(new_gamma-gamma12)/gamma12
        gamma12 = new_gamma

    P2.a = np.sqrt(P2.gamma12*P2.R*P2.T)
    P2.v = P2.M*P2.a
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0
    P2.s = P1.s + P2.cp*np.log(P2.T0/P1.T0) - P2.R*np.log(P2.p0/P1.p0)
    P2.rho = P2.p/P2.R/P2.T
    A = m_dot/P2.rho/P2.v

    return P2, A


def P_compr(m_dot, Cp12, f, P1, P2):
    """
    Computes compressor power from enthalpy rise.
    """
    return m_dot * Cp12 * (P2.T0 - P1.T0)


def P_turb(m_dot, Cp12, f, P1, P2):
    """
    Computes turbine power from enthalpy drop.
    """
    P2.h0 = P1.h0 - Cp12 * (P1.T0 - P2.T0)
    gamma = get_gamma(Cp12)
    P2.p0 = P1.p0 * (P2.T0 / P1.T0)**(gamma/(gamma-1))
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.cp = findCp(P2.T0,f)
    P2.s = P1.s
    return m_dot*Cp12*(P1.T0-P2.T0)


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


def corrected_massflow(mdot_ref,P1,Pref):
    """
    Computes corrected mass flow.
    """
    return mdot_ref / (P1.p0/Pref.p0) / np.sqrt(Pref.T0/P1.T0)

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



def eta_is_from_eta_poly(pi_c, eta_poly, gamma=1.4):
    """
    Converts compressor polytropic efficiency into overall
    isentropic efficiency.

    Relation used:

        eta_is =
            (pi_c^((gamma-1)/gamma) - 1)
            --------------------------------
            (pi_c^((gamma-1)/(gamma*eta_poly)) - 1)

    Args:
        pi_c (float):
            Compressor pressure ratio.

        eta_poly (float):
            Polytropic efficiency.

        gamma (float, optional):
            Specific heat ratio.
            Defaults to 1.4.

    Returns:
        float:
            Isentropic efficiency.
    """

    a = (gamma - 1) / gamma

    return (pi_c**a - 1) / (pi_c**(a / eta_poly) - 1)


def eta_poly_from_eta_is(pi_c, eta_is, gamma=1.4):
    """
    Converts compressor isentropic efficiency into
    polytropic efficiency.

    Args:
        pi_c (float):
            Compressor pressure ratio.

        eta_is (float):
            Isentropic efficiency.

        gamma (float, optional):
            Specific heat ratio.
            Defaults to 1.4.

    Returns:
        float:
            Polytropic efficiency.
    """

    a = (gamma - 1) / gamma

    return (a * np.log(pi_c)) / (np.log(1 + (pi_c**a - 1) / eta_is))

def choking_NPR(gamma):
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
    eta_th = 0.5 * m_dot_a * (v_j**2 - v_1**2) / P
    return eta_th

def overalleff(eta_t,eta_p):
    return eta_p*eta_t