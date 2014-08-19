import unittest
import decay_mod_unittest as decay_mod
import numpy as np

def exact_discrete_solution(n, I, a, theta, dt):
    dt = float(dt)  # avoid integer division
    factor = (1 - (1-theta)*a*dt)/(1 + theta*dt*a)
    return I*factor**n

class TestDecay(unittest.TestCase):

    def setUp(self):
        """Initialization for each test."""
        pass

    def test_exact_discrete_solution(self):
        """
        Compare result from solver against
        formula for the discrete solution.
        """
        theta = 0.8; a = 2; I = 0.1; dt = 0.8
        N = int(8/dt)  # no of steps
        u, t = decay_mod.solver(I=I, a=a, T=N*dt, dt=dt,
                             theta=theta)
        u_de = np.array([exact_discrete_solution(n, I, a, theta, dt)
                         for n in range(N+1)])
        diff = np.abs(u_de - u).max()
        self.assertAlmostEqual(diff, 0, delta=1E-14)

    def test_solver(self):
        """
        Compare result from solver against
        precomputed arrays for theta=0, 0.5, 1.
        """
        precomputed = {
            't': np.array([ 0. ,  0.5,  1. ,  1.5,  2. ,  2.5,
                            3. ,  3.5,  4. ]),
            0.5: np.array(
                [ 0.8       ,  0.43076923,  0.23195266, 0.12489759,
                  0.06725255,  0.03621291,  0.01949926, 0.0104996 ,
                  0.00565363]),
            0: np.array(
                [  8.00000000e-01,   3.20000000e-01,
                   1.28000000e-01,   5.12000000e-02,
                   2.04800000e-02,   8.19200000e-03,
                   3.27680000e-03,   1.31072000e-03,
                   5.24288000e-04]),
            1: np.array(
                [ 0.8       ,  0.5       ,  0.3125    ,  0.1953125 ,
                  0.12207031,  0.07629395,  0.04768372,  0.02980232,
                  0.01862645]),
            }
        for theta in 0, 0.5, 1:
            u, t = decay_mod.solver(I=0.8, a=1.2, T=4, dt=0.5,
                                 theta=theta)
            diff = np.abs(u - precomputed[theta]).max()
            # Precomputed numbers are known to 8 decimal places
            self.assertAlmostEqual(diff, 0, places=8,
                                   msg='theta=%s' % theta)

    def test_potential_integer_division(self):
        """Choose variables that can trigger integer division."""
        theta = 1; a = 1; I = 1; dt = 2
        N = 4
        u, t = decay_mod.solver(I=I, a=a, T=N*dt, dt=dt, theta=theta)
        u_de = np.array([exact_discrete_solution(n, I, a, theta, dt)
                         for n in range(N+1)])
        diff = np.abs(u_de - u).max()
        self.assertAlmostEqual(diff, 0, delta=1E-14)

    def test_convergence_rates(self):
        """Compare empirical convergence rates to exact ones."""
        # Set command-line arguments directly in sys.argv
        import sys
        sys.argv[1:] = '--I 0.8 --a 2.1 --T 5 '\
                       '--dt 0.4 0.2 0.1 0.05 0.025'.split()
        r = decay_mod.main()
        for theta in r:
            self.assertTrue(r[theta])  # check for non-empty list

        expected_rates = {0: 1, 1: 1, 0.5: 2}
        for theta in r:
            r_final = r[theta][-1]
            # Compare to 1 decimal place
            self.assertAlmostEqual(
                expected_rates[theta], r_final, places=1,
                msg='theta=%s' % theta)

if __name__ == '__main__':
    unittest.main()
