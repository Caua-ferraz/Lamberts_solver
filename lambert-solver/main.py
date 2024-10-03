import math

# Vector operations
def vector_add(a, b):
    return [a[i] + b[i] for i in range(3)]

def vector_subtract(a, b):
    return [a[i] - b[i] for i in range(3)]

def vector_multiply(a, scalar):
    return [a[i] * scalar for i in range(3)]

def vector_dot(a, b):
    return sum(a[i] * b[i] for i in range(3))

def vector_norm(a):
    return math.sqrt(sum(x**2 for x in a))

def vector_cross(a, b):
    return [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]
    ]

# Stumpff functions
def stumpff_c(z):
    if z > 0:
        sz = math.sqrt(z)
        return (1 - math.cos(sz)) / z
    elif z < 0:
        sz = math.sqrt(-z)
        return (1 - math.cosh(sz)) / z
    else:
        return 0.5

def stumpff_s(z):
    if z > 0:
        sz = math.sqrt(z)
        return (sz - math.sin(sz)) / (sz**3)
    elif z < 0:
        sz = math.sqrt(-z)
        return (math.sinh(sz) - sz) / (sz**3)
    else:
        return 1/6

# Lambert Solver
class LambertSolver:
    def __init__(self, mu):
        self.mu = mu  # gravitational parameter

    def solve(self, r1, r2, dt, clockwise=False, max_iterations=1000, tolerance=1e-8):
        r1_norm = vector_norm(r1)
        r2_norm = vector_norm(r2)
        
        # Compute the angle between r1 and r2
        cos_dnu = vector_dot(r1, r2) / (r1_norm * r2_norm)
        cos_dnu = max(min(cos_dnu, 1.0), -1.0)  # Clamp to [-1, 1]
        sin_dnu = math.sqrt(1.0 - cos_dnu**2)
        dnu = math.acos(cos_dnu)

        # Check for zero angle to prevent division by zero
        if abs(sin_dnu) < tolerance:
            raise ValueError("Angle between position vectors is zero or very small; cannot compute transfer orbit.")

        # Determine the direction of motion using the cross product
        cross_r1_r2 = vector_cross(r1, r2)
        if (not clockwise and cross_r1_r2[2] < 0) or (clockwise and cross_r1_r2[2] >= 0):
            sin_dnu = -sin_dnu
            dnu = 2 * math.pi - dnu

        # Compute A
        A = sin_dnu * math.sqrt(r1_norm * r2_norm / (1.0 - cos_dnu))
        
        if A == 0:
            raise ValueError("Angle between position vectors is zero; cannot compute transfer orbit.")

        # Initial guess for z
        z = 0.0

        # Function to compute y(z)
        def compute_y(z):
            C = stumpff_c(z)
            S = stumpff_s(z)
            if C == 0:
                return float('inf')
            y = r1_norm + r2_norm + A * (z * S - 1.0) / math.sqrt(C)
            return y

        # Time of flight function
        def time_of_flight(z):
            y = compute_y(z)
            if y < 0:
                return float('inf')
            C = stumpff_c(z)
            S = stumpff_s(z)
            if C == 0:
                return float('inf')
            chi = math.sqrt(y / C)
            tof = (chi**3 * S + A * chi) / math.sqrt(self.mu)
            return tof

        # Use Newton's method to solve for z
        n = 0
        ratio = 1
        while abs(ratio) > tolerance and n < max_iterations:
            n += 1
            y = compute_y(z)
            if y < 0:
                z += 0.1  # Adjust z to try to find positive y
                continue
            C = stumpff_c(z)
            S = stumpff_s(z)
            if C == 0:
                z += 0.1
                continue
            chi = math.sqrt(y / C)
            tof = (chi**3 * S + A * chi) / math.sqrt(self.mu)
            # Compute derivative numerically
            h = 1e-5
            tof_plus = time_of_flight(z + h)
            tof_minus = time_of_flight(z - h)
            dtof_dz = (tof_plus - tof_minus) / (2 * h)
            if dtof_dz == 0:
                break
            ratio = (tof - dt) / dtof_dz
            z = z - ratio

        if n == max_iterations:
            raise ValueError(f"No convergence after {max_iterations} iterations")

        y = compute_y(z)
        f = 1 - y / r1_norm
        g = A * math.sqrt(y / self.mu)
        gdot = 1 - y / r2_norm

        if abs(g) < tolerance:
            raise ValueError("g is too close to zero, causing division issues")

        # Compute velocities
        v1 = vector_multiply(vector_subtract(r2, vector_multiply(r1, f)), 1 / g)
        v2 = vector_multiply(vector_subtract(vector_multiply(r2, gdot), r1), 1 / g)

        return v1, v2

    def earth_to_position(self, target_position, dt, prograde=True):
        """
        Calculate transfer from Earth's surface to a target position.
        
        :param target_position: Target position vector [x, y, z] in km
        :param dt: Time of flight in seconds
        :param prograde: If True, use prograde transfer; if False, use retrograde
        :return: (v1, v2) Initial and final velocity vectors
        """
        earth_radius = 6371  # km
        earth_surface = [earth_radius, 0, 0]  # Starting from Earth's surface on x-axis
        return self.solve(earth_surface, target_position, dt, clockwise=not prograde)

    def leo_to_geo(self, dt, prograde=True):
        """
        Calculate transfer from Low Earth Orbit (LEO) to Geostationary Earth Orbit (GEO).
        
        :param dt: Time of flight in seconds
        :param prograde: If True, use prograde transfer; if False, use retrograde
        :return: (v1, v2) Initial and final velocity vectors
        """
        leo_altitude = 300  # km
        geo_altitude = 35786  # km
        leo_position = [earth_radius + leo_altitude, 0, 0]
        geo_position = [earth_radius + geo_altitude, 0, 0]
        return self.solve(leo_position, geo_position, dt, clockwise=not prograde)

    def earth_to_moon(self, dt, prograde=True):
        """
        Calculate transfer from Earth's surface to Moon's orbit.
        
        :param dt: Time of flight in seconds
        :param prograde: If True, use prograde transfer; if False, use retrograde
        :return: (v1, v2) Initial and final velocity vectors
        """
        earth_surface = [earth_radius, 0, 0]
        moon_distance = 384400  # km
        moon_position = [moon_distance, 0, 0]
        try:
            return self.solve(earth_surface, moon_position, dt, clockwise=not prograde)
        except ValueError as e:
            raise ValueError(f"Error in Earth to Moon transfer: {str(e)}")

def orbital_energy(r, v, mu):
    r_norm = vector_norm(r)
    v_norm = vector_norm(v)
    return v_norm**2 / 2 - mu / r_norm

def vector_add(*vectors):
    return [sum(components) for components in zip(*vectors)]

def propagate_orbit(r0, v0, dt, mu, num_steps=1000):
    """
    Propagate the orbit using the Runge-Kutta 4th order method.

    :param r0: Initial position vector (km)
    :param v0: Initial velocity vector (km/s)
    :param dt: Time to propagate (s)
    :param mu: Gravitational parameter (km^3/s^2)
    :param num_steps: Number of integration steps
    :return: Final position and velocity vectors (km, km/s)
    """
    r = r0[:]
    v = v0[:]
    h = dt / num_steps  # Time step size

    for _ in range(num_steps):
        # Compute k1
        k1_r = v
        acc1 = vector_multiply(r, -mu / vector_norm(r)**3)
        k1_v = acc1

        # Compute k2
        r_temp = vector_add(r, vector_multiply(k1_r, h / 2))
        v_temp = vector_add(v, vector_multiply(k1_v, h / 2))
        acc2 = vector_multiply(r_temp, -mu / vector_norm(r_temp)**3)
        k2_r = v_temp
        k2_v = acc2

        # Compute k3
        r_temp = vector_add(r, vector_multiply(k2_r, h / 2))
        v_temp = vector_add(v, vector_multiply(k2_v, h / 2))
        acc3 = vector_multiply(r_temp, -mu / vector_norm(r_temp)**3)
        k3_r = v_temp
        k3_v = acc3

        # Compute k4
        r_temp = vector_add(r, vector_multiply(k3_r, h))
        v_temp = vector_add(v, vector_multiply(k3_v, h))
        acc4 = vector_multiply(r_temp, -mu / vector_norm(r_temp)**3)
        k4_r = v_temp
        k4_v = acc4

        # Corrected update for position
        delta_r = vector_multiply(
            vector_add(
                k1_r,
                vector_multiply(k2_r, 2),
                vector_multiply(k3_r, 2),
                k4_r
            ),
            h / 6
        )
        r = vector_add(r, delta_r)

        # Corrected update for velocity
        delta_v = vector_multiply(
            vector_add(
                k1_v,
                vector_multiply(k2_v, 2),
                vector_multiply(k3_v, 2),
                k4_v
            ),
            h / 6
        )
        v = vector_add(v, delta_v)

    return r, v

def orbital_period(r, mu):
    """Calculate orbital period for a circular orbit."""
    return 2 * math.pi * math.sqrt(r**3 / mu)

def hohmann_transfer_time(r1, r2, mu):
    """Calculate Hohmann transfer time."""
    a = (r1 + r2) / 2
    return math.pi * math.sqrt(a**3 / mu)

# Constants
earth_radius = 6371  # km
earth_mu = 398600.4418  # km^3/s^2

def main():
    solver = LambertSolver(earth_mu)

    while True:
        try:
            print("\nSelect a scenario:")
            print("1. Earth surface to LEO")
            print("2. LEO to GEO transfer")
            print("3. Earth to Moon transfer")
            print("4. Custom transfer")
            print("5. Exit")

            choice = int(input("Enter your choice (1-5): "))

            if choice == 5:
                print("Exiting the program.")
                break

            if choice not in [1, 2, 3, 4]:
                raise ValueError("Invalid choice. Please enter a number between 1 and 5.")

            if choice == 1:
                r1 = [earth_radius, 0, 0]  # Earth surface
                theta = math.radians(10)   # 10 degrees in radians
                r2 = [
                    (earth_radius + 400) * math.cos(theta),  # 400 km altitude LEO
                    (earth_radius + 400) * math.sin(theta),
                    0
                ]
                dt = 1000  # Time of flight in seconds
                v1, v2 = solver.solve(r1, r2, dt)

            elif choice == 2:
                r1 = [earth_radius + 400, 0, 0]  # LEO (400 km altitude)
                theta = math.radians(10)  # 10 degrees in radians
                r2 = [
                    (earth_radius + 35786) * math.cos(theta),  # GEO
                    (earth_radius + 35786) * math.sin(theta),
                    0
                ]
                dt = 18000  # Time of flight in seconds (5 hours)
                v1, v2 = solver.solve(r1, r2, dt)

            elif choice == 3:
                r1 = [earth_radius, 0, 0]  # Earth surface
                theta = math.radians(10)  # 10 degrees in radians
                r2 = [
                    384400 * math.cos(theta),  # Moon's orbit
                    384400 * math.sin(theta),
                    0
                ]
                dt = 300000  # Time of flight in seconds (about 3.5 days)
                v1, v2 = solver.solve(r1, r2, dt)

            elif choice == 4:
                r1 = [float(x) for x in input("Enter initial position (x y z in km): ").split()]
                r2 = [float(x) for x in input("Enter final position (x y z in km): ").split()]
                dt = float(input("Enter time of flight (s): "))
                clockwise = input("Clockwise transfer? (y/n): ").lower() == 'y'
                v1, v2 = solver.solve(r1, r2, dt, clockwise)

            print("\nResults:")
            print(f"Initial position (km): {r1}")
            print(f"Final position (km): {r2}")
            print(f"Time of flight (s): {dt}")
            print(f"Initial velocity vector (km/s): {v1}")
            print(f"Final velocity vector (km/s): {v2}")

            # Additional calculations
            e1 = orbital_energy(r1, v1, earth_mu)
            e2 = orbital_energy(r2, v2, earth_mu)
            print(f"Initial orbital energy: {e1:.2f} km^2/s^2")
            print(f"Final orbital energy: {e2:.2f} km^2/s^2")
            print(f"Energy difference: {abs(e1 - e2):.2e} km^2/s^2")

            r2_prop, v2_prop = propagate_orbit(r1, v1, dt, earth_mu)
            r2_error = vector_norm(vector_subtract(r2, r2_prop))
            v2_error = vector_norm(vector_subtract(v2, v2_prop))
            print(f"Position error after propagation: {r2_error:.2f} km")
            print(f"Velocity error after propagation: {v2_error:.2f} km/s")
            print(f"Relative position error: {r2_error / vector_norm(r2):.2e}")
            print(f"Relative velocity error: {v2_error / vector_norm(v2):.2e}")

            t_transfer = math.sqrt(vector_norm(vector_subtract(r2, r1))**3 / (8 * earth_mu)) * math.pi
            print(f"Estimated minimum transfer time: {t_transfer:.2f} s")
            print(f"Input transfer time: {dt:.2f} s")

        except ValueError as e:
            print(f"Error: {str(e)}")
        except Exception as e:
            print(f"An unexpected error occurred: {str(e)}")

        print("\nPress Enter to continue...")
        input()

if __name__ == "__main__":
    main()