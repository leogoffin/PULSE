import numpy as np
from .findCp import findCp
from .atmosphere_isa import Air
from .functions import av

def with_Qc(P1, Qc, T0_ref, m_dot, m_dot_e, f,
               P2=None, v=0, delta_p=0, tol=1e-5):
    """
    Solve combustor outlet conditions using energy balance with variable Cp.

    Args:
        P1 (Air): Inlet state object.
        Qc (float): Heat released in combustor (W).
        T0_ref (float): Reference temperature for enthalpy split (K).
        m_dot (float): Inlet air mass flow rate (kg/s).
        m_dot_e (float): Exhaust mass flow rate (kg/s).
        f (float): Fuel-air ratio.
        P2 (Air, optional): Outlet state object. Defaults to None.
        v (int, optional): Air model version. Defaults to 0.
        delta_p (float, optional): Pressure loss fraction. Defaults to 0.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        Air: Outlet state object.
    """
    if P2 == None:
        P2 = Air(v)

    rel_diff = np.inf
    old_cp = P2.get_cp_from_gamma()
    Cp1r = findCp(av(P1.T0, T0_ref), 0)

    while rel_diff > tol:

        P2.T0 = (m_dot * Cp1r * (P1.T0 - T0_ref) + Qc) / (m_dot_e * P2.cp) + T0_ref
        Cp = findCp(av(T0_ref, P2.T0), f)
        P2.cp = Cp
        P2.get_gamma_from_cp(Cp)
        rel_diff = np.abs((old_cp - Cp) / old_cp)
        old_cp = Cp

    P2.cp = findCp(P2.T0, 0)
    P2.get_gamma_from_cp()
    P2.h0 = P2.cp * P2.T0
    P2.p0 = P1.p0 * (1 - delta_p)
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.s = P1.s + P2.cp * np.log(P2.T0 / P1.T0) - P2.R * np.log(P2.p0 / P1.p0)

    return P2


def rev_with_f(P2, delta_hf, T0_ref, m_dot, f,
               P1=None, v=0, delta_p=0, tol=1e-5):
    """
    Solves combustor inlet temperature from outlet conditions using energy balance with variable Cp.

    Args:
        P2 (Air): Outlet state object.
        delta_hf (float): Fuel heating value (J/kg).
        T0_ref (float): Reference temperature for enthalpy split (K).
        m_dot (float): Air mass flow rate (kg/s).
        f (float): Fuel-air ratio.
        P1 (Air, optional): Inlet state object. Defaults to None.
        v (int, optional): Air model version. Defaults to 0.
        delta_p (float, optional): Pressure loss fraction. Defaults to 0.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        Air: Inlet state object.
    """
    if P1 == None:
        P1 = Air(v)

    Cp1r = findCp(av(P2.T0, T0_ref), 0)
    Cp2r = findCp(av(P2.T0, T0_ref), f)

    rel_diff = np.inf

    while rel_diff > tol:

        P1.T0 = T0_ref + (m_dot * (1 + f) * Cp2r * (P2.T0 - T0_ref) - m_dot * f * delta_hf) / (m_dot * Cp1r)

        new_Cp1r = findCp(av(T0_ref, P1.T0), 0)
        rel_diff = np.abs((Cp1r - new_Cp1r) / Cp1r)
        Cp1r = new_Cp1r

    P1.cp = findCp(P1.T0, 0)
    P1.get_gamma_from_cp()
    P1.h0 = P1.cp * P1.T0

    if P2.p0 is not None:
        P1.p0 = P2.p0 / (1 - delta_p)
        P1.rho0 = P1.p0 / (P1.R * P1.T0)

    return P1


def with_Tout(P1, T0_ref, m_dot, Delta_hf, T0_out, f=0, eta_cc=1, pi_cc=1, tol=1e-5):
    """
    Computes fuel-air ratio and outlet state for a specified combustor outlet temperature.

    Args:
        P1 (Air): Inlet state object.
        T0_ref (float): Reference temperature (K).
        m_dot (float): Air mass flow rate (kg/s).
        Delta_hf (float): Fuel heating value (J/kg).
        T0_out (float): Target outlet temperature (K).
        f (float, optional): Initial fuel fraction offset. Defaults to 0.
        eta_cc (float, optional): Combustion efficiency. Defaults to 1.
        pi_cc (float, optional): Pressure ratio. Defaults to 1.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        tuple: (Outlet state, fuel mass flow rate, fuel-air ratio)
    """
    Cp1_r = findCp(av(P1.T0, T0_ref), f)
    f_trial = 0.02 + f
    diff = np.inf

    while diff > tol:

        Cp2_r = findCp(av(T0_out, T0_ref), f_trial)

        m_dot_f = m_dot * (1 + f) * (
            Cp2_r * (T0_out - T0_ref) - Cp1_r * (P1.T0 - T0_ref)
        ) / (eta_cc * Delta_hf - Cp2_r * (T0_out - T0_ref))

        new_f = m_dot_f / m_dot + f
        diff = abs(new_f - f_trial) / f_trial
        f_trial = new_f

    P2 = Air(0)
    P2.p0 = P1.p0 * pi_cc
    P2.T0 = T0_out
    P2.cp = findCp(P2.T0, new_f)
    P2.get_gamma_from_cp()
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0
    P2.s = P1.s + P2.cp * np.log(P2.T0 / P1.T0) - P2.R * np.log(P2.p0 / P1.p0)

    return P2, m_dot_f, new_f


def get_eta(P1, T0_ref, m_dot, f, Delta_hf, T0_out):
    """
    Computes combustion efficiency from temperature rise.

    Args:
        P1 (Air): Inlet state object.
        T0_ref (float): Reference temperature (K).
        m_dot (float): Air mass flow rate (kg/s).
        f (float): Fuel-air ratio.
        Delta_hf (float): Fuel heating value (J/kg).
        T0_out (float): Outlet temperature (K).

    Returns:
        float: Combustion efficiency.
    """
    Cp3r = findCp(av(T0_out, T0_ref), f)
    Cp2r = findCp(av(P1.T0, T0_ref), 0)

    eta_cc = (
        (m_dot * (1 + f) * Cp3r * (T0_out - T0_ref)
         - m_dot * Cp2r * (P1.T0 - T0_ref))
        / Delta_hf / (m_dot * f)
    )

    return eta_cc


def with_f(P1, T0_ref, m_dot, Delta_hf, fadd, fin = 0, eta_cc=1, pi_cc=1, tol=1e-5):
    """
    Compute outlet state after adding fuel to the combustor.

    Args:
        P1 (Air): Inlet state object.
        T0_ref (float): Reference temperature (K).
        m_dot (float): Air mass flow rate (kg/s).
        Delta_hf (float): Fuel heating value (J/kg).
        fadd (float): Added fuel-air ratio.
        fin (float, optional): Initial inlet fuel-air ratio. Defaults to 0.
        eta_cc (float, optional): Combustion efficiency. Defaults to 1.
        pi_cc (float, optional): Combustor pressure ratio. Defaults to 1.
        tol (float, optional): Convergence tolerance. Defaults to 1e-5.

    Returns:
        Air: Outlet state object with updated thermodynamic properties.
    """
    Cp1_r = findCp(av(P1.T0, T0_ref), fin)
    ftot = fadd +fin
    diff = np.inf
    P2 = Air()
    P2.T0 = P1.T0*1.4
    while diff > tol:

        Cp2_r = findCp(av(P2.T0, T0_ref), ftot)
        new_T = T0_ref + (fadd * eta_cc * Delta_hf + (1 + fin) * Cp1_r * (P1.T0 - T0_ref))/ ((1 + ftot) * Cp2_r)
        diff = abs(new_T - P2.T0) / P2.T0
        P2.T0 = new_T

    P2.p0 = P1.p0 * pi_cc
    P2.cp = findCp(P2.T0, ftot)
    P2.get_gamma_from_cp()
    P2.rho0 = P2.p0 / (P2.R * P2.T0)
    P2.h0 = P2.cp * P2.T0
    P2.s = P1.s + P2.cp * np.log(P2.T0 / P1.T0) - P2.R * np.log(P2.p0 / P1.p0)

    return P2
