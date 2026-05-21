import numpy as np
from modules.convergences import bisection
from modules.atmosphere_isa import *

def combustion_analysis(
    species_in,
    nu_in,
    M_in,
    quantities_in,
    species_out,
    nu_out,
    M_out,
    basis="mass_flow",
    print_results=True
):
    """
    Generic combustion / chemical reaction analyzer.

    Parameters
    ----------
    species_in : list[str]
        Reactant names.

    nu_in : list[float]
        Stoichiometric coefficients of reactants.

    M_in : list[float]
        Molecular masses of reactants [kg/mol].

    quantities_in : list[float]
        Known quantities of reactants.
        If basis='mass_flow'  -> kg/s
        If basis='molar_flow' -> mol/s

    species_out : list[str]
        Product names.

    nu_out : list[float]
        Stoichiometric coefficients of products.

    M_out : list[float]
        Molecular masses of products [kg/mol].

    basis : str
        "mass_flow" or "molar_flow"

    Returns
    -------
    dict containing:
        - mass_in
        - molar_in
        - mass_out
        - molar_out
        - mass_fractions_out
        - molar_fractions_out
        - mean_molecular_mass_out
    """

    nu_in = np.array(nu_in, dtype=float)
    nu_out = np.array(nu_out, dtype=float)

    M_in = np.array(M_in, dtype=float)
    M_out = np.array(M_out, dtype=float)

    quantities_in = np.array(quantities_in, dtype=float)

    if basis == "mass_flow":
        mass_in = quantities_in
        molar_in = mass_in / M_in

    elif basis == "molar_flow":
        molar_in = quantities_in
        mass_in = molar_in * M_in

    else:
        raise ValueError("basis must be 'mass_flow' or 'molar_flow'")

    extent_candidates = molar_in / nu_in
    extent = np.min(extent_candidates)

    limiting_index = np.argmin(extent_candidates)
    limiting_species = species_in[limiting_index]

    molar_out = nu_out * extent

    # Remaining reactants (excess species)
    molar_remaining = molar_in - nu_in * extent

    # Add remaining reactants to output if non-zero
    extra_species = []
    extra_molar = []
    extra_mass = []
    extra_M = []

    for i in range(len(species_in)):
        if molar_remaining[i] > 1e-12:
            extra_species.append(species_in[i] + " (excess)")
            extra_molar.append(molar_remaining[i])
            extra_mass.append(molar_remaining[i] * M_in[i])
            extra_M.append(M_in[i])

    mass_out = molar_out * M_out


    final_species = list(species_out) + extra_species
    final_molar = np.concatenate([molar_out, np.array(extra_molar)])
    final_mass = np.concatenate([mass_out, np.array(extra_mass)])
    final_M = np.concatenate([M_out, np.array(extra_M)])

    total_mass_out = np.sum(final_mass)
    total_molar_out = np.sum(final_molar)

    mean_M_out = total_mass_out / total_molar_out

    mass_frac_out = final_mass / total_mass_out
    molar_frac_out = final_molar / total_molar_out

    R_star = 8.31
    R = R_star / mean_M_out

    if print_results:

        # Chemical reaction
        reaction_in = " + ".join(
            [f"{nu_in[i]:g} {species_in[i]}" for i in range(len(species_in))]
        )

        reaction_out = " + ".join(
            [f"{nu_out[i]:g} {species_out[i]}" for i in range(len(species_out))]
        )

        print("\n" + "=" * 70)
        print("CHEMICAL REACTION")
        print("=" * 70)
        print(f"{reaction_in}  -->  {reaction_out}")

        print("\n" + "=" * 70)
        print("INPUT STREAM")
        print("=" * 70)

        for i in range(len(species_in)):
            print(
                f"{species_in[i]:10s} | "
                f"ṁ = {mass_in[i]:10.4f} kg/s | "
                f"ṅ = {molar_in[i]:12.4f} mol/s | "
                f"M = {M_in[i]*1000:8.4f} g/mol"
            )

        print("\nLimiting reactant :", limiting_species)
        print(f"Reaction extent   : {extent:.6f} mol/s")

        print("\n" + "=" * 70)
        print("OUTPUT STREAM")
        print("=" * 70)

        for i in range(len(final_species)):
            print(
                f"{final_species[i]:15s} | "
                f"ṁ = {final_mass[i]:10.4f} kg/s | "
                f"ṅ = {final_molar[i]:12.4f} mol/s | "
                f"Y = {mass_frac_out[i]:8.5f} | "
                f"X = {molar_frac_out[i]:8.5f}"
            )

        print("\n" + "=" * 70)
        print("MIXTURE PROPERTIES")
        print("=" * 70)

        print(f"Total mass flow rate  : {total_mass_out:.4f} kg/s")
        print(f"Total molar flow rate : {total_molar_out:.4f} mol/s")
        print(f"Mean molecular mass   : {mean_M_out*1000:.4f} g/mol")
        print(f"Gas constant          : {R:.3f} J/kg K")

        print("=" * 70 + "\n")

    return {
        "species_out": final_species,
        "mass_in": mass_in,
        "molar_in": molar_in,
        "mass_out": final_mass,
        "molar_out": final_molar,
        "mass_fractions_out": mass_frac_out,
        "molar_fractions_out": molar_frac_out,
        "mean_molecular_mass_out": mean_M_out,
        "reaction_extent": extent,
        "limiting_species": limiting_species,
        "R" : R,
    }

def charac_vel(P1):
    return np.sqrt(P1.gamma*P1.R*P1.T0)/P1.gamma * ((P1.gamma +1)/2)**((P1.gamma + 1)/(2*(P1.gamma-1)))


def f(M,gamma):
    return 1 + (gamma-1)/2 * M**2

def F(M,gamma):
    return np.sqrt(gamma)* M * f(M,gamma)**(-(gamma + 1)/ (2 * (gamma - 1)))


def mdot_nozzle(P1,A,M = None,gamma = None):
    if M is None : 
        M = P1.M
    if gamma is None : 
        gamma = P1.gamma
    return P1.p0*A/(np.sqrt(P1.R*P1.T0))*F(M,gamma)
    
def A_nozzle(P1,mdot,M = None,gamma = None):
    if M is None : 
        M = P1.M
    if gamma is None : 
        gamma = P1.gamma
    return mdot/(P1.p0/(np.sqrt(P1.R*P1.T0))*F(M,gamma))

def M_nozzle(P1,mdot,A,gamma = None,Mmin = 0.4,Mmax = 6,tol = 1e-4):
    if gamma is None : 
        gamma = P1.gamma

    a = Mmin
    b = Mmax

    def error(c):
        return  A_nozzle(P1, mdot, c, gamma) - A
    
    fa = error(a)
    fb = error(b)
    c,its =  bisection(error,a,b)
    print(f"Bisection converged in {its} steps")
    print(f"found root at {c:.3e} - error : {error(c):.3e}")
    return c

def space_thrust(mdot,Pe,Ae):
    return mdot*Pe.v + Pe.p*Ae

def exhaust_cds(P1,Me):
    Pe = Air(R = P1.R)
    Pe.p = P1.p0/f(Me,P1.gamma)**(P1.gamma/(P1.gamma-1))
    Pe.T = P1.T0/f(Me,P1.gamma)
    Pe.gamma = P1.gamma
    Pe.v = Me*np.sqrt(P1.gamma*Pe.R*Pe.T)
    Pe.M = Me
    Pe.a = np.sqrt(P1.gamma*Pe.R*Pe.T)
    return Pe

def fully_expanded(Pe,ptot,mdot,A_e,A_t,K=2.5, g = 9.81, doprint=False):
    pa = K*Pe.p
    T = mdot*Pe.v + Pe.p * (1-K)*A_e
    Isp = T/mdot/g
    c_T = T/ptot/A_t
    if doprint:
        print("\nFULLY EXPANDED NOZZLE")
        print("-----------------------")
        print(f"Ambient pressure factor (K) : {K:.4f}")
        print(f"Exit static pressure        : {pa:.3f} Pa")
        print(f"Thrust                     : {T:.6f} N")
        print(f"Specific impulse           : {Isp:.6f} s")
        print(f"Thrust coefficient         : {c_T:.3f}")
    return pa,T,Isp,c_T

def adapted_operation(Pe,ptot,mdot,A_t,K=2.5, g = 9.81, doprint=False):
    pa = Pe.p
    T = mdot*Pe.v
    Isp = T/mdot/g
    c_T = T/ptot/A_t
    if doprint:
        print("\nADAPTED NOZZLE OPERATION")
        print("-------------------------")
        print(f"Exit static pressure        : {pa:.3f} Pa")
        print(f"Thrust                     : {T:.6f} N")
        print(f"Specific impulse           : {Isp:.6f} s")
        print(f"Thrust coefficient         : {c_T:.3f}")
    return pa,T,Isp,c_T

def M_from_pressure_ratio(p0, p, gamma=1.4):
    return np.sqrt(2/(gamma-1)* ((p0/p)**((gamma-1)/gamma) - 1))

def ground_operation(P1,mdot,A_t,K=2.5, g = 9.81, doprint=False):
    P = Air(0)
    P.set_isa(0)
    Pe = Air(R = P1.R)
    Pe.gamma = P1.gamma
    Pe.p = P.p/K
    Pe.M = M_from_pressure_ratio(P1.p0, P.p, P1.gamma)
    Pe.T = P1.T0 / f(Pe.M, P1.gamma)
    Pe.v = Pe.M * np.sqrt(P1.gamma * P1.R * Pe.T)
    Pe.rho = Pe.R*Pe.T/Pe.p
    Asep = mdot/Pe.v/Pe.rho
    T = mdot*Pe.v + (Pe.p - P.p)*Asep
    Isp = T/mdot/g
    c_T = T/P1.p0/A_t
    if doprint:
        print("\nGROUND OPERATION")
        print("----------------")
        print(f"Separation area            : {Asep:.4f} m²")
        print(f"Exit static pressure       : {Pe.p:.3f} Pa")
        print(f"Thrust                     : {T:.6f} N")
        print(f"Specific impulse           : {Isp:.6f} s")
        print(f"Thrust coefficient         : {c_T:.3f}")

    return Pe,T,Isp,c_T