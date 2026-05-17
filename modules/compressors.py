import numpy as np
from .findCp import findCp
from .atmosphere_isa import Air
from .functions import av, get_gamma

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

def with_eta_poly(P1, m_dot, pi_c, eta_p,
                  P2=None, v=0, tol=1e-5):
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

        P2.T0 = P1.T0 * pi_c**((gamma - 1)/(gamma * eta_p))

        Cp = findCp(av(P1.T0, P2.T0), 0)
        gamma = get_gamma(Cp)

        rel_diff = np.abs((old_cp - Cp) / old_cp)
        old_cp = Cp

    P2.cp = findCp(P2.T0, 0)
    P2.get_gamma_from_cp()

    P2.p0 = pi_c * P1.p0
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0

    P2.s = (
        P1.s
        + P2.cp * np.log(P2.T0 / P1.T0)
        - P2.R * np.log(P2.p0 / P1.p0)
    )

    power = m_dot * Cp * (P2.T0 - P1.T0)

    return P2, power, Cp

def get_eta_poly(P1,P2,pi_c,f=0):
    Cp12 = findCp(av(P1.T0,P2.T0),f)
    g = get_gamma(Cp12)
    return (g-1)/g * np.log(pi_c)/np.log(P2.T0/P1.T0)