import numpy as np
from .findCp import findCp
from .atmosphere_isa import Air
from .functions import *


def with_eta(P1, m_dot, pi_c, eta_c_is,
                       P2=None, v=0, tol=1e-5):
    """
    Solves compressor outlet state using isentropic efficiency with Cp(T)-dependent properties.

    Args:
        P1 (Air): Inlet state.
        m_dot (float): Mass flow rate (kg/s).
        pi_c (float): Compressor pressure ratio (p02/p01).
        eta_c_is (float): Isentropic efficiency.
        P2 (Air, optional): Outlet state object. Defaults to None.
        v (int, optional): Air model version. Defaults to 0.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (P2, compressor power (W), Cp)
    """
    if P2 == None:
        P2 = Air(v)

    rel_diff = np.inf
    old_cp = P1.get_cp_from_gamma()
    gamma = P1.gamma

    while rel_diff > tol:

        P2_T0s = P1.T0 * pi_c**((gamma - 1) / gamma)
        P2.T0 = (P2_T0s - P1.T0) / eta_c_is + P1.T0
        Cp = findCp(av(P1.T0, P2.T0), 0)
        gamma = get_gamma(Cp)

        rel_diff = np.abs((old_cp - Cp) / old_cp)
        old_cp = Cp

    P2.cp = findCp(P2.T0, 0)
    P2.get_gamma_from_cp()
    P2.p0 = pi_c * P1.p0
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0
    P2.s = P1.s + P2.cp * np.log(P2.T0 / P1.T0) - P2.R * np.log(P2.p0 / P1.p0)

    return P2, m_dot * Cp * (P2.T0 - P1.T0), Cp


def get_eta(P1, P2, m_dot, pi_c):
    """
    Computes compressor isentropic efficiency and power based on inlet/outlet states.

    Args:
        P1 (Air): Inlet state.
        P2 (Air): Outlet state.
        m_dot (float): Mass flow rate (kg/s).
        pi_c (float): Pressure ratio (p02/p01).

    Returns:
        tuple: (isentropic efficiency, compressor power (W))
    """

    P1.get_cp_from_gamma()
    Cp12 = findCp(av(P1.T0, P2.T0), 0)
    gamma12 = get_gamma(Cp12)

    P2_T0s = P1.T0 * pi_c**((gamma12 - 1) / gamma12)
    eta_is = (P2_T0s - P1.T0) / (P2.T0 - P1.T0)

    Pc = m_dot * Cp12 * (P2.T0 - P1.T0)
    P2.p0 = P1.p0 * pi_c

    return eta_is, Pc


def corrected(P1, P2, P1_corr, mdot_corr):
    """
    Computes dynamically scaled compressor outlet state using similarity laws.

    Args:
        P1 (Air): Original inlet state.
        P2 (Air): Original outlet state.
        P1_corr (Air): Corrected inlet state.
        mdot_corr (float): Corrected mass flow rate (kg/s).

    Returns:
        tuple: (corrected outlet state, corrected power, Cp)
    """

    pi_c = P2.p0 / P1.p0
    tau_c = P2.T0 / P1.T0

    P2_corr = Air()
    P2_corr.p0 = P1_corr.p0 * pi_c
    P2_corr.T0 = P1_corr.T0 * tau_c

    Cp12_corr = findCp(av(P1_corr.T0, P2_corr.T0), 0)

    P2_corr.cp = findCp(P2_corr.T0, 0)
    P2_corr.get_gamma_from_cp()
    P2_corr.rho0 = P2_corr.p0 / (P2_corr.R * P2_corr.T0)
    P2_corr.h0 = P2_corr.cp * P2_corr.T0
    P2_corr.s = P1_corr.s + P2_corr.cp * np.log(P2_corr.T0 / P1_corr.T0) - P2_corr.R * np.log(P2_corr.p0 / P1_corr.p0)

    return P2_corr, mdot_corr * Cp12_corr * (P2_corr.T0 - P1_corr.T0), Cp12_corr

def with_eta_poly(P1, m_dot, pi_c, eta_p,P2=None, v=0, tol=1e-5):
    """
    Solves compressor outlet state using isentropic efficiency with Cp(T)-dependent properties.

    Args:
        P1 (Air): Inlet state.
        m_dot (float): Mass flow rate (kg/s).
        pi_c (float): Compressor pressure ratio (p02/p01).
        eta_c_is (float): Isentropic efficiency.
        P2 (Air, optional): Outlet state object. Defaults to None.
        v (int, optional): Air model version. Defaults to 0.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (P2, compressor power (W), Cp)
    """
    if P2 == None:
        P2 = Air(v)

    rel_diff = np.inf
    old_cp = P1.get_cp_from_gamma()
    gamma = P1.gamma

    while rel_diff > tol:

        P2_T0s = P1.T0 * pi_c**((gamma - 1) / gamma)
        eta_c_is = eta_is_from_eta_poly(pi_c,eta_p,gamma)
        P2.T0 = (P2_T0s - P1.T0) / eta_c_is + P1.T0
        Cp = findCp(av(P1.T0, P2.T0), 0)
        gamma = get_gamma(Cp)

        rel_diff = np.abs((old_cp - Cp) / old_cp)
        old_cp = Cp

    P2.cp = findCp(P2.T0, 0)
    P2.get_gamma_from_cp()
    P2.p0 = pi_c * P1.p0
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0
    P2.s = P1.s + P2.cp * np.log(P2.T0 / P1.T0) - P2.R * np.log(P2.p0 / P1.p0)

    return P2, m_dot * Cp * (P2.T0 - P1.T0), Cp

def get_eta_poly(P1,P2,pi_c,f=0):
    Cp12 = findCp(av(P1.T0,P2.T0),f)
    g = get_gamma(Cp12)
    return (g-1)/g * np.log(pi_c)/np.log(P2.T0/P1.T0)

def isentropic_temp(T0a,pi_c,gamma):
    return T0a* pi_c**((gamma-1)/gamma)



def with_eta_poly_per_stage(P1, m_dot, pi_c, eta_p,P2=None, v=0, tol=1e-5):
    """
    Solves compressor outlet state using polytropic efficiency
    with Cp(T)-dependent thermodynamic properties.

    The compressor outlet total temperature is computed from the
    polytropic compression relation:

        T02 / T01 = (p02 / p01)^((gamma - 1)/(gamma * eta_p))

    Since Cp and gamma vary with temperature, the solution is
    iterated until Cp converges.

    Args:
        P1 (Air):
            Inlet stagnation state.

        m_dot (float):
            Mass flow rate [kg/s].

        pi_c (float):
            Compressor total pressure ratio:
                pi_c = p02 / p01

        eta_p (float):
            Compressor polytropic efficiency.

        P2 (Air, optional):
            Outlet state object.
            If None, a new Air object is created.

        v (int, optional):
            Air model version identifier.
            Defaults to 0.

        tol (float, optional):
            Relative Cp convergence tolerance.
            Defaults to 1e-5.

    Returns:
        tuple:
            (
                P2 (Air): outlet stagnation state,
                power (float): compressor power [W],
                Cp (float): converged average Cp [J/kg/K]
            )

    Notes:
        - Uses an average Cp evaluated between inlet and outlet
          temperatures.
        - Assumes ideal gas behavior.
        - Suitable for preliminary gas turbine cycle calculations.
        - For very high pressure ratios, numerical integration of
          differential polytropic compression may be more accurate.
    """

    if P2 is None:
        P2 = Air(v)

    rel_diff = np.inf
    old_cp = P1.get_cp_from_gamma()
    gamma = P1.gamma

    while rel_diff > tol:
        eta_is = (pi_c**((gamma-1)/gamma)-1)/(pi_c**((gamma-1)/(gamma*eta_p))-1)
        T1s = P1.T0 * pi_c**((gamma - 1)/(gamma))
        P2.T0 = P1.T0 + (T1s-P1.T0)/eta_is
        Cp = findCp(av(P1.T0, P2.T0), 0)
        gamma = get_gamma(Cp)
        rel_diff = np.abs((old_cp - Cp) / old_cp)
        old_cp = Cp


    P2.cp = findCp(P2.T0, 0)
    P2.get_gamma_from_cp()
    P2.p0 = pi_c * P1.p0
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0

    P2.s = (P1.s + P2.cp * np.log(P2.T0 / P1.T0)- P2.R * np.log(P2.p0 / P1.p0))

    power = m_dot * Cp * (P2.T0 - P1.T0)

    return P2, power, Cp, eta_is



def IterateMatchingFan_and_SecondaryNozzle(P1,Pref, pi_fs_values, eta_fs_values, m_dot_star_values, BPR, A_nozzle, tol = 1e-4,max_iter = 100):
    """We assume the nozzle remains choked (as in primary flow), so we can use similarity
       fan_pi and eta_s_fan are tables used to interpolate wrt to m_dot. Select 4 points 
       from operating map to put in tables"""
    NPR_crit = 1.89

    f_eta = interp1d(m_dot_star_values, eta_fs_values, kind='linear', bounds_error=False, fill_value='extrapolate')
    f_pi  = interp1d(m_dot_star_values, pi_fs_values,  kind='linear', bounds_error=False, fill_value='extrapolate')

    def nozzle_mass_flow(m_dot_star):
        pi_fs  = float(f_pi(m_dot_star))
        eta_fs = float(f_eta(m_dot_star))
        p_2_tot = pi_fs * P1.p0
        NPR = p_2_tot / Pref.p
        P2s,Pcs,Cp12s = with_eta(P1,m_dot_star,pi_fs,eta_fs)
        gamma12 = get_gamma(Cp12s)
        return ChokedMassFlow(P2s, gamma12, A_nozzle)

    m_star_1 = float(np.min(m_dot_star_values)) 
    m_star_2 = float(np.max(m_dot_star_values))

    delta_1  = true_massflow_from_corrected(m_star_1*BPR/(1+BPR),P1,Pref) - nozzle_mass_flow(m_star_1)
    delta_2  = true_massflow_from_corrected(m_star_2*BPR/(1+BPR),P1,Pref) - nozzle_mass_flow(m_star_2)

    if delta_1 * delta_2 > 0:
        raise ValueError(  f"Root not bracketed: Δm1={delta_1:.4f}, Δm2={delta_2:.4f} "  f"at m_dot_star=[{m_star_1:.1f}, {m_star_2:.1f}]")

    m_star_m  = None
    delta_m_m = None

    for i in range(max_iter):
        m_star_m  = (m_star_1 + m_star_2) / 2.0
        delta_m_m = true_massflow_from_corrected(m_star_m*BPR/(1+BPR),P1,Pref) - nozzle_mass_flow(m_star_m)

        m_dot_m_s = m_star_m *BPR/(1+BPR)
        if abs(delta_m_m) < tol * m_dot_m_s:
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
    m_dot_real  = m_star_m * (P1.p0 / Pref.p) / np.sqrt(P1.T0 / Pref.T)
    m_dot_s_conv = m_dot_real * BPR /(1+ BPR)
    P2 = Air()
    P2,Pc,Cp12 = with_eta(P1,m_dot_real,pi_fs_conv,eta_fs_conv)

    return m_dot_real, m_dot_s_conv, pi_fs_conv, eta_fs_conv, P2