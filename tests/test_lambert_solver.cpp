#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include "../include/lamberts_solver.h"

// Helper function to compare vectors
bool compare_vectors(const Vector3& a, const Vector3& b, double tolerance = 1e-3) {
    return std::abs(a.x - b.x) < tolerance &&
           std::abs(a.y - b.y) < tolerance &&
           std::abs(a.z - b.z) < tolerance;
}

// Test case structure
struct TestCase {
    Vector3 r1;
    Vector3 r2;
    double delta_t;
    bool is_prograde;
    Vector3 expected_v1;
    Vector3 expected_v2;
    std::string description;
};

// Run a single test case
bool run_test(const TestCase& test, LambertSolver& solver) {
    try {
        std::cout << "\nRunning test case: " << test.description << std::endl;
        std::cout << "Solving Lambert's problem..." << std::endl;
        auto result = solver.solve(test.r1, test.r2, test.delta_t, test.is_prograde);
        std::cout << "Solution found." << std::endl;

        // Print computed velocities
        std::cout << "Computed v1: (" << result.first.x << ", " << result.first.y << ", " << result.first.z << ") km/s" << std::endl;
        std::cout << "Computed v2: (" << result.second.x << ", " << result.second.y << ", " << result.second.z << ") km/s" << std::endl;

        // Compare with expected velocities if provided
        if (test.expected_v1.x != 0 || test.expected_v1.y != 0 || test.expected_v1.z != 0) {
            bool v1_match = compare_vectors(result.first, test.expected_v1);
            bool v2_match = compare_vectors(result.second, test.expected_v2);

            if (v1_match && v2_match) {
                std::cout << "Test passed!" << std::endl;
                return true;
            } else {
                std::cout << "Test failed." << std::endl;
                std::cout << "Expected v1: ("
                          << test.expected_v1.x << ", " << test.expected_v1.y << ", " << test.expected_v1.z
                          << ") km/s, Got: ("
                          << result.first.x << ", " << result.first.y << ", " << result.first.z << ") km/s" << std::endl;
                std::cout << "Expected v2: ("
                          << test.expected_v2.x << ", " << test.expected_v2.y << ", " << test.expected_v2.z
                          << ") km/s, Got: ("
                          << result.second.x << ", " << result.second.y << ", " << result.second.z << ") km/s" << std::endl;
                return false;
            }
        } else {
            // No expected velocities provided
            std::cout << "Test passed (no expected velocities provided)." << std::endl;
            return true;
        }
    } catch (const std::exception& e) {
        std::cerr << "Test threw an exception: " << e.what() << std::endl;
        return false;
    } catch (...) {
        std::cerr << "Test threw an unknown exception" << std::endl;
        return false;
    }
}

int main() {
    // Earth's gravitational parameter (km^3/s^2)
    double mu = 398600.4418;
    LambertSolver solver(mu);

    // Define test cases
    std::vector<TestCase> tests = {
        // Test case 1: Example from Vallado's book (Example 5-5)
        {
            // Initial position vector r1 (km)
            { -6045.0, -3490.0, 2500.0 },
            // Final position vector r2 (km)
            { -3738.0, 3000.0, 5000.0 },
            // Time of flight (s)
            3600.0,
            // Prograde motion
            true,
            // Expected velocities at r1 (km/s)
            { 3.56736, 6.45633, 1.69248 },
            // Expected velocities at r2 (km/s)
            { 5.69075, 3.03686, -1.57662 },
            // Description of the test case
            "Vallado Example 5-5"
        },
        // Test case 2: Another example with known solution
        {
            // Initial position vector r1 (km)
            { 5000.0, 10000.0, 2100.0 },
            // Final position vector r2 (km)
            { -14600.0, 2500.0, 7000.0 },
            // Time of flight (s)
            3600.0,
            // Prograde motion
            true,
            // Expected velocities at r1 (km/s)
            { -5.9925, 1.9254, 3.2456 },
            // Expected velocities at r2 (km/s)
            { -3.3125, -4.1966, -0.38529 },
            // Description of the test case
            "Curtis Example 5.2"
        },
        // Test case 3: Short transfer
        {
            // Initial position vector r1 (km)
            { 7000.0, 0.0, 0.0 },
            // Final position vector r2 (km)
            { 0.0, 7000.0, 0.0 },
            // Time of flight (s)
            1800.0,
            // Prograde motion
            true,
            // Expected velocities at r1 (km/s)
            { -3.3740, 5.6713, 0.0 },
            // Expected velocities at r2 (km/s)
            { -5.6713, -3.3740, 0.0 },
            // Description of the test case
            "Short transfer (quarter orbit)"
        },
        // Test case 4: Long transfer
        {
            // Initial position vector r1 (km)
            { 7000.0, 0.0, 0.0 },
            // Final position vector r2 (km)
            { 0.0, -7000.0, 0.0 },
            // Time of flight (s)
            5400.0,
            // Prograde motion
            true,
            // Expected velocities at r1 (km/s)
            { 1.9441, 3.2745, 0.0 },
            // Expected velocities at r2 (km/s)
            { 3.2745, -1.9441, 0.0 },
            // Description of the test case
            "Long transfer (three-quarter orbit)"
        }
    };

    int passed = 0;
    int total = tests.size();

    for (size_t i = 0; i < tests.size(); ++i) {
        std::cout << "\n==============================" << std::endl;
        if (run_test(tests[i], solver)) {
            passed++;
        }
    }

    std::cout << "\n==============================" << std::endl;
    std::cout << "Passed " << passed << " out of " << total << " tests." << std::endl;

    return (passed == total) ? 0 : 1;
}
