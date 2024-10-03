#include <iostream>
#include "../include/lamberts_solver.h"

int main() {
    // Standard gravitational parameter for Earth (km^3/s^2)
    double mu = 398600.4418;

    // Define initial and final position vectors (in kilometers)
    Vector3 r1 = {7000.0, -12124.0, 0.0};
    Vector3 r2 = {12457.0, 0.0, 0.0};

    // Time of flight in seconds
    double delta_t = 3600.0; // 1 hour

    try {
        // Initialize the Lambert solver
        LambertSolver solver(mu);

        // Solve Lambert's problem (prograde)
        auto velocities = solver.solve(r1, r2, delta_t, true);
        Vector3 v1 = velocities.first;
        Vector3 v2 = velocities.second;

        std::cout << "Initial Velocity: (" << v1.x << ", " << v1.y << ", " << v1.z << ") km/s\n";
        std::cout << "Final Velocity:   (" << v2.x << ", " << v2.y << ", " << v2.z << ") km/s\n";
    } catch (const std::exception& e) {
        std::cerr << "Error solving Lambert's problem: " << e.what() << std::endl;
        
        // Add debugging information
        std::cout << "Debug Information:\n";
        std::cout << "mu: " << mu << " km^3/s^2\n";
        std::cout << "r1: (" << r1.x << ", " << r1.y << ", " << r1.z << ") km\n";
        std::cout << "r2: (" << r2.x << ", " << r2.y << ", " << r2.z << ") km\n";
        std::cout << "Time of flight: " << delta_t << " s\n";
        
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}