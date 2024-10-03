#include "../include/lamberts_solver.h"
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <sstream>

const double PI = 3.14159265358979323846;

LambertSolver::LambertSolver(double mu_) : mu(mu_) {}

double LambertSolver::norm(const Vector3& v) const {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

Vector3 LambertSolver::cross(const Vector3& a, const Vector3& b) const {
    return Vector3{
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x
    };
}

double LambertSolver::dot(const Vector3& a, const Vector3& b) const {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

std::pair<Vector3, Vector3> LambertSolver::solve(const Vector3& r1, const Vector3& r2, double tof, bool is_prograde) {
    const double tolerance = 1e-8;
    const int max_iterations = 1000;

    double r1_norm = norm(r1);
    double r2_norm = norm(r2);

    if (r1_norm == 0 || r2_norm == 0) {
        throw std::invalid_argument("Position vectors cannot be zero.");
    }

    // Compute the transfer angle
    double cos_dtheta = dot(r1, r2) / (r1_norm * r2_norm);
    cos_dtheta = std::max(-1.0, std::min(1.0, cos_dtheta)); // Clamp to avoid numerical errors
    double sin_dtheta = norm(cross(r1, r2)) / (r1_norm * r2_norm);

    // Determine the direction of motion
    Vector3 cross_prod = cross(r1, r2);
    if (is_prograde) {
        if (cross_prod.z < 0) {
            sin_dtheta = -sin_dtheta;
        }
    } else {
        if (cross_prod.z >= 0) {
            sin_dtheta = -sin_dtheta;
        }
    }

    double dtheta = std::atan2(sin_dtheta, cos_dtheta);
    if (dtheta < 0) {
        dtheta += 2 * PI;
    }

    // Compute A
    double A = sin_dtheta * std::sqrt(r1_norm * r2_norm / (1 - cos_dtheta));

    if (A == 0) {
        throw std::runtime_error("Cannot compute A.");
    }

    // Initial guess for z and bounds
    double z = (dtheta > PI) ? -4 * PI * PI : 4 * PI * PI;
    double z_low = -4 * PI * PI;
    double z_up = 4 * PI * PI;

    double y, C, S, F, dFdz;
    int iteration = 0;
    bool converged = false;

    while (iteration < max_iterations) {
        iteration++;

        C = stumpff_C(z);
        S = stumpff_S(z);

        y = r1_norm + r2_norm + A * (z * S - 1) / std::sqrt(C);

        if (y < 0) {
            // Implement bracketing method
            double z_low = z;
            double z_high = z < 0 ? 0 : z * 2;
            while (y < 0) {
                z = (z_low + z_high) / 2;
                C = stumpff_C(z);
                S = stumpff_S(z);
                y = r1_norm + r2_norm + A * (z * S - 1) / std::sqrt(C);
                if (z < 0) {
                    z_low = z;
                } else {
                    z_high = z;
                }
            }
            continue;
        }

        double sqrt_y = std::sqrt(y);
        F = y * S + A * sqrt_y - std::sqrt(mu) * tof;

        if (std::abs(F) < tolerance) {
            converged = true;
            break;
        }

        dFdz = (y * C + A * S * sqrt_y) / (2.0 * C);

        double damping_factor = 0.5;
        double delta_z = damping_factor * F / dFdz;

        // Limit step size to prevent excessive jumps
        double max_step = 10.0; // Define a reasonable maximum step size
        if (std::abs(delta_z) > max_step) {
            delta_z = (delta_z > 0.0) ? max_step : -max_step;
        }

        z = z - delta_z;

        // Prevent z from becoming excessively large
        if (std::abs(z) > 1e6) { // Adjust the bound as necessary
            std::cerr << "z exceeded reasonable bounds at iteration " << iteration << std::endl;
            throw std::runtime_error("z exceeded reasonable bounds.");
        }

        z = std::max(z_low, std::min(z, z_up));

        // Debugging output
        if (iteration % 10 == 0 || std::abs(F) < tolerance) {
            std::cout << "Iteration " << iteration << ": "
                      << "z = " << z << ", "
                      << "C(z) = " << C << ", "
                      << "S(z) = " << S << ", "
                      << "y(z) = " << y << ", "
                      << "F(z) = " << F << ", "
                      << "dF/dz = " << dFdz << std::endl;
        }
    }

    if (!converged) {
        std::stringstream error_msg;
        error_msg << "Lambert solver did not converge after " << max_iterations << " iterations.\n"
                  << "Final values: F = " << F << ", z = " << z << ", y = " << y << "\n"
                  << "Input parameters: r1 = (" << r1.x << ", " << r1.y << ", " << r1.z << "), "
                  << "r2 = (" << r2.x << ", " << r2.y << ", " << r2.z << "), "
                  << "tof = " << tof << ", is_prograde = " << std::boolalpha << is_prograde;
        throw std::runtime_error(error_msg.str());
    }

    // Compute Lagrange coefficients
    double f = 1 - y / r1_norm;
    double g = A * std::sqrt(y / mu);
    double g_dot = 1 - y / r2_norm;

    // Compute velocity vectors
    Vector3 v1 = {
        (r2.x - f * r1.x) / g,
        (r2.y - f * r1.y) / g,
        (r2.z - f * r1.z) / g
    };

    Vector3 v2 = {
        (g_dot * r2.x - r1.x) / g,
        (g_dot * r2.y - r1.y) / g,
        (g_dot * r2.z - r1.z) / g
    };

    return { v1, v2 };
}

double LambertSolver::stumpff_C(double z) const {
    if (std::abs(z) < 1e-6) {
        return 1.0 / 2.0 - z / 24.0 + z * z / 720.0 - z * z * z / 40320.0;
    } else if (z > 0) {
        return (1.0 - std::cos(std::sqrt(z))) / z;
    } else {
        return (std::cosh(std::sqrt(-z)) - 1.0) / (-z);
    }
}

double LambertSolver::stumpff_S(double z) const {
    if (std::abs(z) < 1e-6) {
        return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - z * z * z / 362880.0;
    } else if (z > 0) {
        return (std::sqrt(z) - std::sin(std::sqrt(z))) / std::pow(z, 1.5);
    } else {
        return (std::sinh(std::sqrt(-z)) - std::sqrt(-z)) / std::pow(-z, 1.5);
    }
}