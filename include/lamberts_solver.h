#ifndef LAMBERTS_SOLVER_H
#define LAMBERTS_SOLVER_H

#include <utility>
#include <stdexcept>
#include <cmath>

struct Vector3 {
    double x;
    double y;
    double z;
};

class LambertSolver {
public:
    /**
     * Constructor
     * @param mu Standard gravitational parameter (km^3/s^2)
     */
    LambertSolver(double mu);

    /**
     * Solves Lambert's problem using the Universal Variable method.
     * @param r1 Initial position vector (km)
     * @param r2 Final position vector (km)
     * @param tof Time of flight (s)
     * @param is_prograde True for prograde, false for retrograde
     * @return Pair of velocity vectors at r1 and r2 (km/s)
     */
    std::pair<Vector3, Vector3> solve(
        const Vector3& r1,
        const Vector3& r2,
        double tof,
        bool is_prograde = true
    );

private:
    double mu; // Standard gravitational parameter

    // Helper functions
    double norm(const Vector3& v) const;
    Vector3 cross(const Vector3& a, const Vector3& b) const;
    double dot(const Vector3& a, const Vector3& b) const;

    // Stumpff functions
    double stumpff_C(double z) const;
    double stumpff_S(double z) const;
};

#endif // LAMBERTS_SOLVER_H
