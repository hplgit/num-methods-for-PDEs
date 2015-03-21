import decay
import nose.tools as nt
import numpy as np

def test_third2():
    x = 0.15
    expected = (1/3.)*x
    computed = decay.third(x)
    tol = 1E-15
    success = abs(expected - computed) < tol
    assert success

import numpy as np

def test_exact_discrete_solution2():
    """
    Compare result from solver against
    formula for the discrete solution.
    """
    theta = 0.8; a = 2; I = 0.1; dt = 0.8
    Nt = int(8/dt)  # no of steps
    u, t = decay.solver(I=I, a=a, T=Nt*dt, dt=dt, theta=theta)
    u_de = np.array(
        [decay.exact_discrete_solution(n, I, a, theta, dt)
         for n in range(Nt+1)])
    diff = np.abs(u_de - u).max()
    nt.assert_almost_equal(diff, 0, delta=1E-14)
    tol = 1E-14
    success = diff < tol
    assert success

def test_potential_integer_division2():
    """Choose variables that can trigger integer division."""
    theta = 1; a = 1; I = 1; dt = 2
    Nt = 4
    u, t = decay.solver(I=I, a=a, T=Nt*dt, dt=dt, theta=theta)
    u_de = np.array(
        [decay.exact_discrete_solution(n, I, a, theta, dt)
         for n in range(Nt+1)])
    diff = np.abs(u_de - u).max()
    tol = 1E-14
    nt.assert_almost_equal(diff, 0, delta=1E-14)
