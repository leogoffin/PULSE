import numpy as np
from pystdatm import temperature, pressure, density, speed_of_sound,viscosity

class Air:
    """
    Represents a flowing air/gas parcel at a given velocity.

    Attach any GasModel via .set_state() to compute static and total
    thermodynamic properties, plus Mach, Re, dynamic pressure,
    entropy, and enthalpy.

    Args:
        v (float): Flow velocity (m/s).
    """

    def __init__(self, v: float = 0):
        self.v = v

        # Static (ambient) properties
        self.rho      = None   # kg/m³
        self.a        = None   # m/s
        self.T        = None   # K
        self.p        = None   # Pa
        self.mu       = None   # Pa·s   dynamic viscosity
        self.nu       = None   # m²/s   kinematic viscosity
        self.cp       = None   # J/kg·K
        self.gamma    = None
        self.R        = 287.05 # J/kg·K  specific gas constant
        self.y        = None   # m      altitude
        # Flow quantities
        self.M        = None
        self.q        = None   # Pa     dynamic pressure

        # Stagnation (total) conditions
        self.T0       = None   # K
        self.p0       = None   # Pa
        self.rho0     = None   # kg/m³

        # Calorific quantities  (ref: 0 K, ideal gas)
        self.h        = None   # J/kg   static specific enthalpy
        self.h0       = None   # J/kg   stagnation specific enthalpy
        self.s        = 0      # J/kg·K static specific entropy  (relative to ISA sea-level)


    def get_cp_from_gamma(self,gamma = None):
        if gamma is not None : 
            self.gamma = gamma
        if self.gamma is None : 
            print("Warning : Gamma not defined, will default to 1.4")
            self.gamma = 1.4
        self.cp = self.gamma * self.R / (self.gamma - 1)
        return self.cp

    def get_gamma_from_cp(self,cp = None):
        if cp is not None : 
            self.cp = cp
        if self.cp is None : 
            print("Error : Cp not defined, cannot compute gamma")
            return
        if self.R is None :
            print("Warning : R not defined, will default to 287,05 J/(kg.K)")
        self.gamma = self.cp/(self.cp-self.R)
        return self.gamma

    def set_isa(self,y, dT_ISA=0.0):
        """
        Compute ISA atmospheric properties up to 32 km altitude.

        Based on the International Standard Atmosphere (ISA) model, with an optional
        temperature offset from standard conditions.

        Args:
            h_m (float): Altitude in meters.
            dT_ISA (float, optional): Temperature deviation from ISA standard (°C).
                Positive values indicate warmer-than-standard conditions.
                Defaults to 0.0.

        Returns:
            tuple:
                rho (float): air density (kg/m^3).

                a (float): speed of sound a (m/s).

                T (float): Temperature (K).

                Mu (float): dynamic viscosity mu (Pa.s).

        """
        p0   = 101325.0    # Pa
        T_ref  = 288.15     # K

        T = temperature(y)
        p = pressure(y)
        rho = density(y)
        a = speed_of_sound(y)
        mu = viscosity(T)

        if self.gamma is None : 
            self.gamma = 1.4

        self.p = p
        self.rho = rho
        self.mu = mu
        self.nu = self.mu/self.rho
        self.a = a
        self.T = T
        self.M = self.v/self.a
        self.y = y
        self.q = 1/2*self.rho*self.v**2


        
        if self.cp is None :
            self.get_cp_from_gamma() 

        # Thermodynamics
        self.h = self.cp * self.T
        self.h0 = self.h + 0.5 * self.v**2
        # Entropy (relative to ISA sea level)
        self.s = self.cp * np.log(self.T / T_ref) - self.R * np.log(self.p / p0)
    
    def set_total(self):
        """
        Sets Total properties from air conditions  (isentropic)
        """
        if self.gamma is None : 
            print("Warning : Gamma not defined, will default to 1.4")
            self.gamma = 1.4
        
        if self.M is None :
            if self.T is None : 
                print("Error : Temperature not defined, cannot compute Mach number ")
                return
            if self.v is None : 
                print("Warning : airspeed not defined, will default to 0 m/s ")
                self.v = 0
            self.a = np.sqrt(self.gamma*self.R*self.T)
            self.M = self.v/self.a


        factor  = 1 + (self.gamma - 1) / 2 * self.M**2
        if self.T is not None:
            self.T0 = self.T * factor
        if self.p is not None:
            self.p0   = self.p   * factor ** (self.gamma / (self.gamma - 1))
        if self.rho is not None:
            self.rho0 = self.rho * factor ** (1 / (self.gamma - 1))
        if self.h is not None:
            self.h0 = self.h + 0.5 * self.v**2
        elif self.T0 is not None : 
            self.h0 = self.T0* self.get_cp_from_gamma()

    def set_static(self):
        """
        Sets static properties from stagnation (total) conditions (isentropic)
        """

        if self.gamma is None:
            print("Warning : Gamma not defined, will default to 1.4")
            self.gamma = 1.4

        # Need Mach number
        if self.M is None:
            if self.T is None:
                print("Error : Need either Mach or temperature")
                return
            if self.v is None:
                print("Warning : airspeed not defined, will default to 0 m/s")
                self.v = 0
            self.a = np.sqrt(self.gamma * self.R * self.T)
            self.M = self.v / self.a

        factor = 1 + (self.gamma - 1) / 2 * self.M**2


        if self.T0 is not None:
            self.T = self.T0 / factor

        if self.p0 is not None:
            self.p = self.p0 / factor ** (self.gamma / (self.gamma - 1))

        if self.rho0 is not None:
            self.rho = self.rho0 / factor ** (1 / (self.gamma - 1))

        # Update derived quantities
        if self.T is not None:
            self.a = np.sqrt(self.gamma * self.R * self.T)

        if self.a is not None and self.M is not None:
            self.v = self.M * self.a

        if self.rho0 is not None and self.v is not None:
            self.q = 0.5 * self.rho * self.v**2
        if self.h0 is not None:
            self.h = self.h0 - 0.5 * self.v**2

            
    def set_viscosity(self):
        """
        Computes dynamic (mu) and kinematic (nu) viscosity
        using Sutherland's law.

        Requires:
            self.T (temperature in K)
            self.rho (density) for kinematic viscosity
        """

        if self.T is None:
            print("Error: Temperature not defined, cannot compute viscosity")
            return

        # Sutherland constants for air
        mu_ref = 1.716e-5   # Pa·s
        T_ref  = 273.15     # K
        S      = 110.4      # K

        # Dynamic viscosity (Sutherland's law)
        self.mu = mu_ref * (self.T / T_ref)**1.5 * (T_ref + S) / (self.T + S)

        # Kinematic viscosity
        if self.rho is not None:
            self.nu = self.mu / self.rho
        else:
            self.nu = None
            print("Warning: Density not defined, nu not computed")

        return self.mu, self.nu

    def is_complete(self):
        required = [self.rho, self.T, self.p, self.M]
        return all(v is not None for v in required)

    def RAM(self,RR = 1):
        P1 = Air(self.v)
        P1.T0 = self.T0
        P1.rho0 = self.rho0
        P1.v = 0
        P1.M = 0
        P1.p0 = RR * self.p0
        P1.s = self.s
        P1.h0 = self.h0
        P1.set_static()
        return P1

    def complete_thermo(self):
        """
        Compute the third stagnation property from two known stagnation values.
        
        For ideal gas, stagnation conditions follow p0 = rho0 * R * T0.
        Provide exactly two of: p0 (Pa), T0 (K), rho0 (kg/m³).
        """

        if self.p0 is not None and self.T0 is not None:
            self.rho0 = self.p0 / (self.R * self.T0)
        elif self.p0 is not None and self.rho0 is not None:
            self.T0 = self.p0 / (self.rho0 * self.R)
        elif self.T0 is not None and self.rho0 is not None:
            self.p0 = self.rho0 * self.R * self.T0
        else:
            raise ValueError("Unexpected combination of provided values.")
        self.set_static()

    def print_state(self):
        """
        Pretty print all available flow properties.
        """

        def fmt(val, unit="", precision=3):
            if val is None:
                return "—"
            if abs(val) < 1e-3 or abs(val) > 1e5:
                return f"{val:.{precision}e} {unit}".strip()
            return f"{val:.{precision}f} {unit}".strip()

        print("\n" + "="*50)
        print(" AIR FLOW STATE ")
        print("="*50)
        if self.is_complete :
            print("-- FULL STATE KNOWN --")
        print("\n--- Flow Conditions ---")
        print(f"Velocity (v)        : {fmt(self.v, 'm/s')}")
        print(f"Mach number (M)     : {fmt(self.M)}")
        print(f"Dynamic pressure (q): {fmt(self.q, 'Pa')}")

        print("\n--- Static Properties ---")
        print(f"Density (rho)       : {fmt(self.rho, 'kg/m³')}")
        print(f"Temperature (T)     : {fmt(self.T, 'K')}")
        print(f"Pressure (p)        : {fmt(self.p, 'Pa')}")
        print(f"Speed of sound (a)  : {fmt(self.a, 'm/s')}")
        print(f"Dynamic visc (mu)   : {fmt(self.mu, 'Pa·s', 6)}")
        print(f"Kinematic visc (nu) : {fmt(self.nu, 'm²/s', 6)}")

        print("\n--- Gas Properties ---")
        print(f"cp                  : {fmt(self.cp, 'J/kg·K')}")
        print(f"gamma               : {fmt(self.gamma)}")
        print(f"R                   : {fmt(self.R, 'J/kg·K')}")

        print("\n--- Stagnation (Total) ---")
        print(f"T0                  : {fmt(self.T0, 'K')}")
        print(f"p0                  : {fmt(self.p0, 'Pa')}")
        print(f"rho0                : {fmt(self.rho0, 'kg/m³')}")

        print("\n--- Energy ---")
        print(f"Enthalpy (h)        : {fmt(self.h, 'J/kg')}")
        print(f"Total enthalpy (h0) : {fmt(self.h0, 'J/kg')}")
        print(f"Entropy (s)         : {fmt(self.s, 'J/kg·K')}")

        print("\n--- Misc ---")
        print(f"Altitude (y)        : {fmt(self.y, 'm')}")

        print("="*50 + "\n")



