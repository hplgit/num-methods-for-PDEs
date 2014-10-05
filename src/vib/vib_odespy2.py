import scitools.std as plt
#import matplotlib.pyplot as plt
import sys
import numpy as np
from scitools.std import StringFunction
from numpy import *  # handy for StringFunction with pi, cos, sin, etc.

class RHS:
    """m*u'' + b*|u'|u' + s(u) = F(t)."""

    def __init__(self, m=1, b=0, s=lambda u: 2*np.pi**2*u,
                 F=lambda t: 0, I=1, V=0, damping='linear'):
        if isinstance(s, str):
            # Turn string formula into Python function
            # (send globals() such that pi, sin, cos, etc from
            # numpy can be used in the string expression)
            s = StringFunction(s, globals=globals())
        if isinstance(F, str):
            F = StringFunction(F, globals=globals())

        self.m, self.b, self.s, self.F = float(m), b, s, F
        self.I, self.V = I, V
        self.damping = damping

    def __call__(self, u, t):
        """Right-hand side function defining the ODE."""
        u, v = u  # u is array of length 2 holding our [u, v]
        if self.damping == 'linear':
            b_term = self.b*v
        elif self.damping == 'quadratic':
            b_term = self.b*np.abs(v)*v
        else:
            b_term = 0
        return [v, (-b_term - self.s(u) + self.F(t))/self.m]


def run_solvers_and_plot(solvers, rhs, T, dt, title=''):
    Nt = int(round(T/float(dt)))
    t_mesh = np.linspace(0, T, Nt+1)
    t_fine = np.linspace(0, T, 8*Nt+1)  # used for very accurate solution

    legends = []
    solver_exact = odespy.RK4(rhs)

    for solver in solvers:
        solver.set_initial_condition([rhs.I, 0])
        u, t = solver.solve(t_mesh)

        solver_name = 'CrankNicolson' if solver.__class__.__name__ == \
                      'MidpointImplicit' else solver.__class__.__name__

        if len(t_mesh) <= 50:
            plt.plot(t, u[:,0])             # markers by default
        else:
            plt.plot(t, u[:,0], '-2')       # no markers
        plt.hold('on')
        legends.append(solver_name)

    # Compare with RK4 on a much finer mesh
    solver_exact.set_initial_condition([rhs.I, 0])
    u_e, t_e = solver_exact.solve(t_fine)

    plt.plot(t_e, u_e[:,0], '-') # avoid markers by spec. line type
    legends.append('exact (RK4, dt=%g)' % (t_fine[1]-t_fine[0]))
    plt.legend(legends, loc='upper right')
    plt.xlabel('t');  plt.ylabel('u')
    plt.title(title)
    plotfilestem = '_'.join(legends)
    plt.savefig('tmp_%s.png' % plotfilestem)
    plt.savefig('tmp_%s.pdf' % plotfilestem)

