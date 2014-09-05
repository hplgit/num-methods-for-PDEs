import numpy as np
import matplotlib.pyplot as plt
from numpy import exp

def cooling(T0, k, T_s, t_end, dt, theta=0.5):
    """
    Solve T'=-k(T-T_s(t)), T(0)=T0,
    for t in (0,t_end] with steps of dt.
    T_s(t) is a Python function of t.
    theta=0.5 means Crank-Nicolson, 1 is Backward
    Euler, and 0 is Forward Euler scheme.
    """
    dt = float(dt)                  # avoid integer division
    Nt = int(round(t_end/dt))       # no of time intervals
    t_end = Nt*dt                   # adjust to fit time step dt
    u = np.zeros(Nt+1)              # array of u[n] values
    t = np.linspace(0, t_end, Nt+1) # time mesh
    u[0] = T0                       # set initial condition
    for n in range(0, Nt):          # n=0,1,...,Nt-1
        u[n+1] = ((1 - dt*(1 - theta)*k)*u[n] + \
        dt*k*(theta*T_s(t[n+1]) + (1 - theta)*T_s(t[n])))/ \
        (1 + dt*theta*k)
    return u, t


def test_asymptotic():
    """
    Test that any initial condition leads to
    the same asymptotic behavior when T_s=constant.
    """
    plt.figure()
    T_s = 5.
    k = 1.2
    dt = 0.1
    tol = 0.01  # tolerance for testing asymptotic value
    t_end = 7    # make sure t_end is large enough for tol
    T0_values = [0, 2, 4, 5, 6, 8, 10] # test many cases

    for T0 in [0, 2, 4, 5, 6, 8, 10]:
        u, t = cooling(T0, k, lambda t: T_s, t_end, dt)
        plt.plot(t, u)

        assert abs(u[-1] - T_s) < tol, '%s != %s' % (u[-1], T_s)

    plt.legend(['T0=%g' % T0 for T0 in T0_values])
    plt.title('Testing asymptotic behavior T_s=%g' % T_s)
    plt.xlabel('t')
    plt.ylabel('T')
    plt.show()



class Piecewise(object):
    """
    Class for holding a piecewise constant T_s function
    and the exact solution of the problem for this case.
    """
    def __init__(self, T0, k, C0, C1, t_star):
        self.T0 = T0
        self.k = k
        self.C0, self.C1 = C0, C1
        self.t_star = t_star
        self.computed_exact_solution_with_sympy = False


    def __call__(self, t):
        """
        Return value of piecewise constant function.
        t can be float or numpy array.
        """
        if isinstance(t, (float,int)):
            if t <= self.t_star:
                T_s = self.C0
            elif t > self.t_star:
                T_s = self.C1
        else:
            # assume numpy array
            T_s = np.piecewise(t,
                               [t <= self.t_star, t > self.t_star],
                               [self.C0, self.C1])
            # Alternative
            # T_s = np.where(t <= self.t_star, C0, C1)
        return T_s


    def exact_solution_trapezidal_int(self, t):
        """
        Return exact (analytical) solution of the problem.
        t can be float. NOTE: have not tested numpy array.
        Exact solution is produced by the trapezoidal rule
        for integration.
        """
        k, T0, C0, C1 = self.k, self.T0, self.C0, self.C1
        t = np.linspace(0, t, 1000)  # integration points
        dt = t[1]-t[0]
        T = T0*exp(-k*t) + exp(-k*t)*k*dt*np.sum(
            self.T_s(t)*exp(k*t))
        return T

    def exact_solution(self, t, verbose=False):
        """
        Return exact (analytical) solution of the problem.
        t can be float or numpy array.
        Exact solution is produced by sympy.
        """
        if not self.computed_exact_solution_with_sympy:
            self.precompute_exact_solution(verbose=True)

        # self.exact0/2 works with t as numpy array
        if isinstance(t, (float,int)):
            if t < self.t_star:
                return self.exact0(
                    t, self.C0, self.k, self.T0)
            else:
                return self.exact1(
                    t, self.C0, self.C1, self.t_star, self.k, self.T0)
        else:
            # assume numpy array
            return np.where(
                t < self.t_star,
                self.exact0(t, self.C0, self.k, self.T0),
                self.exact1(t, self.C0, self.C1,
                            self.t_star, self.k, self.T0))

    def precompute_exact_solution(self, verbose=False):
        """Compute the exact solution via exact integration in sympy."""

        # sol1: solution for t < t_star, sol2: solution for t > t_star
        import sympy as sp
        T0 = sp.symbols('T0')
        k = sp.symbols('k', positive=True)
        # Piecewise linear T_sunction
        t, t_star, C0, C1 = sp.symbols('t t_star C0 C1')
        T_s = C0
        I = sp.integrate(sp.exp(k*t)*T_s, (t, 0, t))
        sol1 = T0*sp.exp(-k*t) + k*sp.exp(-k*t)*I
        sol1 = sp.simplify(sp.expand(sol1))
        if verbose:
            # Some debugging print
            print 'solution t < t_star:', sol1
            #print sp.latex(sol1)
        T_s = C1
        I = sp.integrate(sp.exp(k*t)*C0, (t, 0, t_star)) + \
            sp.integrate(sp.exp(k*t)*C1, (t, t_star, t))
        sol2 = T0*sp.exp(-k*t) + k*sp.exp(-k*t)*I
        sol2 = sp.simplify(sp.expand(sol2))
        if verbose:
            print 'solution t > t_star:', sol2
            #print sp.latex(sol2)

        # Convert to numerical functions
        self.exact0 = sp.lambdify([t, C0, k, T0],
                                  sol1, modules='numpy')
        self.exact1 = sp.lambdify([t, C0, C1, t_star, k, T0],
                                  sol2, modules='numpy')
        self.computed_exact_solution_with_sympy = True


def check_piecewise():
    """
    Compare exact and numerical solution with piecewise
    constant surrounding temperature.
    """
    T0 = 5.
    k = 1.1
    C0 = 2*T0
    C1 = 0.5*T0
    t_star = 4./k
    t_end = 10

    T_s = Piecewise(T0, k, C0, C1, t_star)

    plt.figure()
    dt_values = [1, 0.1]
    for dt in dt_values:
        u, t = cooling(T0, k, T_s, t_end, dt)
        #u_e = [T_s.exact_solution(t_) for t_ in t]  # scalar version
        plt.plot(t, u)

    u_e = T_s.exact_solution(t)
    plt.plot(t, u_e)
    plt.legend(['CN, dt=%g' % dt for dt in dt_values] + ['exact'])
    plt.title('T(t) for piecewise constant $T_s(t)$')
    plt.xlabel('t')
    plt.ylabel('T')
    plt.savefig('tmp.png')
    plt.show()



def illustrate_instability():
    """Demonstrate oscillations in the
    Crank-Nicolson scheme for the piecewise constant case."""

    T0 = 10000. # Very high so that k can also be large and hence dt smaller.
    k = 6
    C0 = 2*T0
    C1 = 0.5*T0
    t_star = 20./k
    t_end = 10
    dt = 0.4
    dt = 0.8

    T_s = Piecewise(T0, k, C0, C1, t_star)

    u, t = cooling(T0, k, T_s, t_end, dt)
    t_e = np.linspace(0, t_end, 1001)  # fine mesh for exact solution
    u_e = T_s.exact_solution(t_e)

    plt.figure()
    plt.plot(t, u, t_e, u_e)
    plt.legend(['CN', 'exact'])
    plt.title(' Numerical instability with Crank Nicolson ')
    plt.xlabel('t')
    plt.ylabel('T')
    plt.show()


def test_discrete_solution():
    """ Compares the numerical solution with an exact solution of the scheme
    when the T_s is constant """

    T_s = 10
    T0 = 2
    k = 1.2
    t_end = 5
    dt = 0.1

    u, t = cooling(T0, k, lambda t: T_s , t_end, dt)
    dt = t[1] - t[0]
    A = (1 - k*dt/2.)/(1 + k*dt/2.)
    u_discrete_exact = T_s + (T0-T_s)*A**(np.arange(len(t)))
    diff = np.abs(u - u_discrete_exact).max()
    print 'diff between computed and exact discrete solution:', diff
    tol = 1E-13
    success = diff < tol
    assert success, 'diff=%g' % diff

# Visually ok curves but with bugs, see

if __name__=='__main__':
    test_asymptotic()
    check_piecewise()
    test_discrete_solution()
    illustrate_instability()
