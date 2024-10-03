**Note:** This is a **test project** intended for educational and public use only. It should **NOT** be used for work purposes.

## Introduction

Welcome to the **Lambert Solver** project! This Python script provides a collection of functions and classes to perform orbital mechanics calculations, including solving Lambert's problem, propagating orbits, and computing orbital energies.

## Features

- **Vector Operations:** Basic vector arithmetic functions such as addition, subtraction, multiplication, dot product, cross product, and norms.
- **Stumpff Functions:** Calculations of Stumpff functions \( C(z) \) and \( S(z) \) used in orbital mechanics.
- **Lambert Solver:** A class to solve Lambert's problem, determining the orbit connecting two position vectors in a given time.
- **Orbital Calculations:**
  - Compute orbital energy.
  - Propagate orbits using the Runge-Kutta 4th order method.
  - Calculate orbital periods and Hohmann transfer times.
- **Scenarios:** Predefined scenarios for common orbital maneuvers:
  - Earth surface to Low Earth Orbit (LEO)
  - LEO to Geostationary Earth Orbit (GEO) transfer
  - Earth to Moon transfer
  - Custom transfers

## Getting Started

### Prerequisites

- Python 3.x
- No external libraries are required; only the standard math library is used.

### Running the Script

1. **Clone or Download the Repository**

   ```bash
   git clone https://github.com/yourusername/the_small_guy_property.git
   ```

2. **Navigate to the Project Directory**

   ```bash
   cd the_small_guy_property
   ```

3. **Run the Script**

   ```bash
   python main.py
   ```

   Replace `main.py` with the actual name of the Python script.

### Usage

Upon running the script, you'll be presented with a menu to select a scenario:

```
Select a scenario:
1. Earth surface to LEO
2. LEO to GEO transfer
3. Earth to Moon transfer
4. Custom transfer
5. Exit
Enter your choice (1-5):
```

- **Options 1-3:** Run predefined orbital transfer scenarios.
- **Option 4:** Input custom initial and final position vectors, time of flight, and specify the transfer direction.
- **Option 5:** Exit the program.

### Example

**Custom Transfer:**

1. Choose option 4.
2. Input the initial position vector:

   ```
   Enter initial position (x y z in km): 7000 0 0
   ```

3. Input the final position vector:

   ```
   Enter final position (x y z in km): 42164 0 0
   ```

4. Input the time of flight:

   ```
   Enter time of flight (s): 21600
   ```

5. Specify the transfer direction:

   ```
   Clockwise transfer? (y/n): n
   ```

The script will output the results, including initial and final velocities, orbital energies, errors in propagation, and estimated transfer times.

## Disclaimer

This is a **test project** and is intended **solely for educational and public use**. It should **NOT** be used for any work-related or professional purposes. The calculations and methods provided are simplified and may not account for all real-world variables. Always consult with a professional aerospace engineer for accurate orbital mechanics computations.

## Contributing

Contributions are welcome! Feel free to fork the repository and submit pull requests.

## License

This project is open-source and available for public use under the [MIT License](LICENSE).

## Contact

For any questions or feedback, please open an issue in the repository or contact the maintainer.

---

**Remember:** This is a test project for public use. Do not use it for work or professional purposes.