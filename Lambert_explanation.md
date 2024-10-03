# Understanding Lambert's Problem

## Introduction

Lambert's Problem is a fundamental concept in orbital mechanics and astrodynamics. It involves determining the orbit connecting two points in space within a specified time. Essentially, given two position vectors and the time of flight between them, Lambert's Problem helps calculate the velocity vectors required to move from the first position to the second in the given time.

## Purpose and Applications

Lambert's Problem is crucial for mission planning and spacecraft navigation. It allows engineers and scientists to:

- **Plan Interplanetary Missions:** Calculate transfer orbits between planets or celestial bodies.
- **Satellite Maneuvers:** Design orbit transfers for satellites, such as moving from Low Earth Orbit (LEO) to Geostationary Earth Orbit (GEO).
- **Trajectory Optimization:** Optimize fuel consumption by finding the most efficient paths.
- **Orbital Rendezvous:** Determine the trajectories needed for spacecraft to meet or dock with other objects.

## How Lambert's Problem Works

### Inputs and Outputs

**Inputs:**

- **Initial Position Vector (r₁):** The starting point in space.
- **Final Position Vector (r₂):** The destination point in space.
- **Time of Flight (Δt):** The time it takes to travel from r₁ to r₂.
- **Gravitational Parameter (μ):** A constant related to the central body's mass (e.g., Earth's gravitational parameter).

**Outputs:**

- **Initial Velocity Vector (v₁):** The velocity required at r₁ to reach r₂ in time Δt.
- **Final Velocity Vector (v₂):** The velocity at r₂ upon arrival.

### Mathematical Background

At its core, Lambert's Problem involves solving the orbital boundary-value problem, which is derived from Kepler's laws of planetary motion. The main challenge is that there isn't a straightforward analytical solution, so numerical methods are often employed.

**Key Concepts:**

- **Conic Sections:** Orbits are conic sections (ellipses, parabolas, hyperbolas) determined by two points and the energy of the orbit.
- **Stumpff Functions:** Special functions used to handle the mathematical expressions involving universal variables.
- **Universal Variable Formulation:** A method that generalizes the solution for all types of orbits.

### Computational Methods

The solution involves iterative numerical techniques:

1. **Initial Guess:** Start with an initial estimate for a parameter (commonly denoted as 'z').
2. **Stumpff Functions Calculation:** Compute the Stumpff functions \( C(z) \) and \( S(z) \) based on 'z'.
3. **Time of Flight Computation:** Calculate the time of flight using the current 'z' value.
4. **Newton-Raphson Iteration:** Adjust 'z' using numerical methods to minimize the difference between the computed and desired time of flight.
5. **Convergence Check:** Repeat the iteration until the solution converges within a specified tolerance.

## Usage in the Code

The provided Python script implements a Lambert Solver with the following components:

- **Vector Operations:** Functions for vector arithmetic necessary for orbital calculations.
- **Stumpff Functions:** Functions to compute \( C(z) \) and \( S(z) \) for any real value of 'z'.
- **LambertSolver Class:**
  - **Initialization:** Takes the gravitational parameter 'μ'.
  - **solve Method:** Implements the iterative method to solve Lambert's Problem.
  - **Predefined Methods:** Functions like `earth_to_moon` simplify common scenarios.
- **Orbital Functions:**
  - **orbital_energy:** Computes the specific orbital energy.
  - **propagate_orbit:** Propagates the orbit using numerical integration (Runge-Kutta method).

### Example Workflow

1. **Input Positions and Time:**
   - Define the initial and final position vectors.
   - Specify the time of flight.
2. **Call the Solver:**
   - Use `solver.solve(r1, r2, dt)` to compute the velocity vectors.
3. **Analyze Results:**
   - Examine the initial and final velocities.
   - Compute orbital energies and check for conservation.
4. **Propagation and Error Checking:**
   - Propagate the orbit to validate the solution.
   - Calculate errors to assess accuracy.

## Limitations and Considerations

- **Numerical Methods:** The solution relies on iterative numerical methods, which may not converge in certain cases.
- **Assumptions:**
  - Two-body problem: Assumes only two bodies are influencing the motion.
  - Instantaneous velocity changes: Real maneuvers may require finite burn times.
- **Accuracy:** Depends on the tolerance settings and the numerical stability of the methods used.
- **Time of Flight:** There may be multiple solutions for certain times of flight (short and long paths).

## Conclusion

Lambert's Problem is essential for planning and executing space missions. Understanding its principles and computational methods enables more efficient and accurate trajectory designs. This script provides a practical implementation for educational purposes, allowing users to explore and visualize orbital transfers.

## References

- **Orbital Mechanics for Engineering Students** by Howard D. Curtis
- **Fundamentals of Astrodynamics** by Roger R. Bate, Donald D. Mueller, and Jerry E. White
- **NASA's Basics of Space Flight:** [https://solarsystem.nasa.gov/basics/](https://solarsystem.nasa.gov/basics/)

---

**Note:** This explanation is intended for educational and public use only. It simplifies complex concepts for better understanding and should not be used for professional or work-related purposes.