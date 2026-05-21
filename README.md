# PULSE

The library solves propulsion systems element-by-element, propagating complete thermodynamic states through each component while maintaining a flexible and extensible architecture.

Designed primarily for preliminary engine analysis, and research prototyping, provides reusable building blocks for constructing complete propulsion cycles with physically consistent thermodynamic solutions.

### Core Capabilities

- Compressors
- Turbines
- Combustion chambers
- Intakes
- Nozzles
- Choked flow conditions
- Compressible quasi-1D flow relations
- Thermodynamic state propagation
- Cycle coupling through power and mass-flow balance

The framework is built to allow rapid extension with custom elements, maps, losses, efficiencies, or alternative thermodynamic models.

### Thermodynamic Architecture

At the center of the framework is the Air state object.

The Air class stores the complete thermodynamic state of a flow at a given location in the propulsion system, including variables such as :

- Static and stagnation pressure
- Static and stagnation temperature
- Mach number
- Velocity
- Density
- Enthalpy
- Entropy
- Gas properties
- Composition-dependent quantities

This approach allows each propulsion element to operate independently while maintaining physically consistent state propagation throughout the engine.

### Typical Applications

- Turbojet and turbofan preliminary analysis
- Rocket propulsion calculations
- Thermodynamic cycle studies

### Academic Background

This project was developed in the context of the Aerospace Propulsion - AERO0014 curriculum taught by Koen Hillewaert at the University of Liège.

If this repository contributes to academic work, derivative projects, publications, or teaching material, please acknowledge both:

- This repository and its author
- The Aerospace Propulsion course material by Koen Hillewaert

### Use of AI

Generative AI has been used in accordance to the University AI chart througout this project, namely for debugging, print formating and docstrings.
