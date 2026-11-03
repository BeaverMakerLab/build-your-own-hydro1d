"""Run: python -m unittest discover -s code -p 'test_*.py' -v"""
import unittest
import numpy as np
from solver import make_dambreak, simulate, physical_flux, rusanov_flux
from exact_solution import exact_dambreak, intermediate_state


class TestTeachingSolver(unittest.TestCase):
    def test_uniform_still_water(self):
        U0 = np.column_stack((np.full(100, 1.25), np.zeros(100)))
        result = simulate(U0, dx=5., t_final=10.)
        np.testing.assert_allclose(result.state, U0, atol=5e-14, rtol=0)

    def test_uniform_moving_water(self):
        U0 = np.column_stack((np.full(100, 1.25), np.full(100, 1.25*0.6)))
        result = simulate(U0, dx=5., t_final=10.)
        np.testing.assert_allclose(result.state, U0, atol=5e-14, rtol=0)

    def test_constant_flux(self):
        U0 = np.array([[1.3, 0.0]])
        np.testing.assert_allclose(rusanov_flux(U0, U0), physical_flux(U0))

    def test_dambreak_mass_budget_and_positive_depth(self):
        _, dx, U0 = make_dambreak(nx=200)
        res = simulate(U0, dx, 60., snapshot_times=(10.,30.))
        self.assertGreater(res.min_depth, 0.0)
        self.assertLess(abs(res.mass_balance_error), 1e-9)
        self.assertLess(abs(res.mass_final-res.mass_initial), 1e-9)

    def test_exact_plateau_and_convergence(self):
        middle = intermediate_state()
        self.assertAlmostEqual(middle['h_middle'], 1.453840892374573, places=9)
        errors = []
        for nx in (100,200,400):
            x, dx, u0 = make_dambreak(nx=nx)
            res = simulate(u0, dx, 60.)
            exact_h, _ = exact_dambreak(x, 60.)
            errors.append(np.mean(np.abs(res.state[:,0] - exact_h)))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])
        self.assertLess(errors[-1], 0.08)


if __name__ == '__main__':
    unittest.main()
