import numpy as np
from .findCp import findCp
from .atmosphere_isa import Air
from .functions import av, get_gamma

def with_P(P_T, m_dot, f, Pin=None, TIT=None, eta_is=1, tol=1e-5):
    """
    Computes turbine outlet state using power balance with Cp(T)-dependent properties.

    Args:
        P_T (float): Turbine power (W).
        m_dot (float): Mass flow rate (kg/s).
        f (float): Fuel-air ratio.
        Pin (Air, optional): Inlet state object. Defaults to None.
        TIT (float, optional): Turbine inlet temperature (K). Used if Pin is not provided.
        eta_is (float, optional): Isentropic efficiency. Defaults to 1.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (Outlet state, converged Cp, isentropic outlet temperature)
    """
    if Pin is None:
        if TIT is not None:
            Pin = Air()
            Pin.T0 = TIT
            Pin.cp = findCp(TIT, f)
        else:
            print("Error : neither TIT nor inlet condition (P1) were given")

    if Pin.cp is None:
        Pin.cp = findCp(Pin.T0, f)

    Cp_trial = Pin.cp
    diff = np.inf

    while diff > tol:
        T2_0 = Pin.T0 - P_T / (m_dot * (1 + f) * Cp_trial)
        Cp_new = findCp(av(Pin.T0, T2_0), f)
        diff = abs(Cp_new - Cp_trial) / Cp_trial
        Cp_trial = Cp_new
        gamma = get_gamma(Cp_trial)

    T2s = Pin.T0 - (Pin.T0 - T2_0) / eta_is

    P2 = Air()
    P2.T0 = T2_0

    if Pin.p0 is not None:
        P2.p0 = (T2s / Pin.T0)**(gamma / (gamma - 1)) * Pin.p0
        P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.cp = findCp(P2.T0, f)
    P2.get_gamma_from_cp()
    P2.h0 = P2.cp * P2.T0

    if Pin.s is not None and Pin.p0 is not None:
        P2.s = Pin.s + P2.cp * np.log(P2.T0 / Pin.T0) - P2.R * np.log(P2.p0 / Pin.p0)

    return P2, Cp_new, T2s


def get_eta(P1, P2, m_dot, f):
    """
    Computes turbine isentropic efficiency and power output.

    Args:
        P1 (Air): Inlet state object.
        P2 (Air): Outlet state object.
        m_dot (float): Mass flow rate (kg/s).
        f (float): Fuel-air ratio.

    Returns:
        tuple: (isentropic efficiency, turbine power (W), isentropic outlet temperature)
    """

    P1.get_cp_from_gamma()
    Cp12 = findCp(av(P1.T0, P2.T0), f)
    gamma12 = get_gamma(Cp12)

    P2_T0s = P1.T0 * (P2.p0 / P1.p0)**((gamma12 - 1) / gamma12)

    eta_is = (P1.T0 - P2.T0) / (P1.T0 - P2_T0s)

    Pt = m_dot * Cp12 * (P1.T0 - P2.T0)

    return eta_is, Pt, P2_T0s


def with_P_poly(P_T, m_dot, f,
                Pin=None,
                TIT=None,
                eta_p=1,
                tol=1e-5):
    """
    Computes turbine outlet state using power balance with
    polytropic efficiency and Cp(T)-dependent properties.

    The turbine outlet temperature is obtained from the
    shaft power balance, while outlet pressure is computed
    using the polytropic expansion relation.

    Args:
        P_T (float):
            Turbine power output [W].

        m_dot (float):
            Air mass flow rate [kg/s].

        f (float):
            Fuel-air ratio.

        Pin (Air, optional):
            Turbine inlet stagnation state.

        TIT (float, optional):
            Turbine inlet temperature [K].
            Used only if Pin is not provided.

        eta_p (float, optional):
            Turbine polytropic efficiency.
            Defaults to 1.

        tol (float, optional):
            Relative Cp convergence tolerance.
            Defaults to 1e-5.

    Returns:
        tuple:
            (
                P2 (Air): outlet stagnation state,
                Cp (float): converged average Cp,
                pi_t (float): turbine pressure ratio p02/p01
            )

    Notes:
        - Uses average Cp between inlet and outlet temperatures.
        - Assumes ideal-gas behavior.
        - Outlet temperature is determined directly from power balance.
        - Outlet pressure is obtained from the polytropic relation.
    """

    if Pin is None:

        if TIT is not None:
            Pin = Air()
            Pin.T0 = TIT
            Pin.cp = findCp(TIT, f)

        else:
            raise ValueError("Neither TIT nor inlet state Pin was provided.")

    if Pin.cp is None:
        Pin.cp = findCp(Pin.T0, f)

    Cp_trial = Pin.cp
    diff = np.inf
    gamma = Pin.gamma

    while diff > tol:

        T2_0 = (Pin.T0- P_T / (m_dot * (1 + f) * Cp_trial))
        Cp_new = findCp(av(Pin.T0, T2_0), f)
        gamma = get_gamma(Cp_new)
        diff = abs(Cp_new - Cp_trial) / Cp_trial
        Cp_trial = Cp_new

    pi_t = (T2_0 / Pin.T0)**(gamma / (eta_p * (gamma - 1)))
    P2 = Air()
    P2.T0 = T2_0
    P2.cp = findCp(P2.T0, f)
    P2.get_gamma_from_cp()

    if Pin.p0 is not None:
        P2.p0 = pi_t * Pin.p0

    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0

    if Pin.s is not None and Pin.p0 is not None:
        P2.s = (Pin.s+ P2.cp * np.log(P2.T0 / Pin.T0)- P2.R * np.log(P2.p0 / Pin.p0))

    return P2, Cp_new, pi_t