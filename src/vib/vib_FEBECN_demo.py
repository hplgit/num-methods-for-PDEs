import odespy
from vib_odespy import run_solvers_and_plot, RHS, \
     VibSolverWrapper4Odespy, plt
from numpy import pi, sin

# Primary ODE: m=1, s(u)=(2*pi)**2*u, such that the period is 1.
# Then we add linear damping and a force term A*sin(w*t) where
# w is half and double of the frequency of the free oscillations.

ODEs = [
    (RHS(b=0.1), 'Small damping, no forcing'),
    (RHS(b=0.4), 'Medium damping, no forcing'),
    (RHS(b=0.4, F=lambda t: 1*sin(0.5*pi*t)),
     'Medium damping, medium forcing w/smaller frequency'),
    (RHS(b=0.4, F=lambda t: 10*sin(0.5*pi*t)),
     'Medium damping, large forcing w/smaller frequency'),
    (RHS(b=1.2, F=lambda t: 10*sin(0.5*pi*t)),
     'Strong damping, large forcing w/smaller frequency'),
    (RHS(b=0.4, F=lambda t: 1*sin(2*pi*t)),
     'Medium damping, medium forcing w/larger frequency'),
    (RHS(b=0.4, F=lambda t: 10*sin(2*pi*t)),
     'Medium damping, large forcing w/larger frequency'),
    (RHS(b=1.2, F=lambda t: 10*sin(2*pi*t)),
     'Strong damping, large forcing w/larger frequency'),
    ]

for rhs, title in ODEs:
    solvers = [
        odespy.ForwardEuler(rhs),
        # Implicit methods must use Newton solver to converge
        odespy.BackwardEuler(rhs, nonlinear_solver='Newton'),
        odespy.CrankNicolson(rhs, nonlinear_solver='Newton'),
        VibSolverWrapper4Odespy(rhs),
        ]

    T = 20          # Period is 1
    dt = 0.05       # 20 steps per period
    filename = 'FEBNCN_' + title.replace(', ', '_').replace('w/', '')
    title = title + ' (dt=%g)' % dt
    plt.figure()
    run_solvers_and_plot(solvers, rhs, T, dt, title=title,
                         filename=filename)

plt.show()
raw_input()
