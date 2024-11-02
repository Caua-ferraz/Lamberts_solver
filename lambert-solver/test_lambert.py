import math
import unittest
from main import LambertSolver, vector_norm, vector_cross, vector_dot

class TestLambertSolver(unittest.TestCase):
    def setUp(self):
        self.mu = 398600.4418  # Earth's gravitational parameter (km³/s²)
        self.solver = LambertSolver(self.mu)
        
    def test_hohmann_transfer(self):
        """Test a simple Hohmann transfer orbit"""
        r1 = [7000.0, 0.0, 0.0]  # 7000 km circular orbit
        r2 = [-42000.0, 100.0, 0.0]  # GEO-like orbit with sufficient y-offset
        
        # Hohmann transfer time (half an orbit)
        dt = math.pi * math.sqrt(((7000.0 + 42000.0)/2.0)**3 / self.mu)
        
        v1, v2 = self.solver.solve(r1, r2, dt)
        
        # Initial circular orbit velocity
        v_circ = math.sqrt(self.mu / r1[0])
        
        # Expected delta-v for Hohmann transfer
        dv1_expected = math.sqrt(self.mu/r1[0]) * (math.sqrt(2*42000.0/(7000.0 + 42000.0)) - 1)
        
        # Compute delta-v from solver's result
        delta_v1 = v1[1] - v_circ
        
        self.assertAlmostEqual(delta_v1, dv1_expected, places=1)  # Compare delta-v

    def test_physical_constraints(self):
        """Test that solutions respect physical constraints"""
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]  # 90-degree transfer
        dt = 3600.0  # 1 hour transfer
        
        v1, v2 = self.solver.solve(r1, r2, dt)
        
        # Check that velocities are reasonable (less than escape velocity)
        escape_velocity = math.sqrt(2 * self.mu / 7000.0)
        self.assertLess(vector_norm(v1), escape_velocity)
        self.assertLess(vector_norm(v2), escape_velocity)
        
    def test_angular_momentum_conservation(self):
        """Test conservation of angular momentum"""
        r1 = [7000.0, 0.0, 0.0]
        r2 = [0.0, 7000.0, 0.0]
        dt = 3600.0
        
        v1, v2 = self.solver.solve(r1, r2, dt)
        
        # Calculate angular momentum vectors at start and end
        h1 = vector_cross(r1, v1)
        h2 = vector_cross(r2, v2)
        
        # Angular momentum magnitude should be conserved
        self.assertAlmostEqual(vector_norm(h1), vector_norm(h2), places=3)
        
    def test_minimum_energy(self):
        """Test that transfer uses reasonable energy"""
        r1 = [7000.0, 0.0, 0.0]
        r2 = [-6999.9, 100.0, 0.0]  # Almost 180-degree transfer with sufficient offset
        
        # Time for minimum energy transfer (π * sqrt(r³/μ))
        dt = math.pi * math.sqrt(7000.0**3 / self.mu)
        
        v1, v2 = self.solver.solve(r1, r2, dt)
        
        # Energy should be negative (bound orbit)
        E = vector_norm(v1)**2 / 2 - self.mu / vector_norm(r1)
        self.assertLess(E, 0)

if __name__ == '__main__':
    unittest.main() 