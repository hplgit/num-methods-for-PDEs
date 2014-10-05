import scitools.std as plt
#import matplotlib.pyplot as plt
from vib_empirical_analysis import minmax, periods, amplitudes
import sys
import odespy
import numpy as np

def f(u, t, w=1):
    u, v = u  # u is array of length 2 holding our [u, v]
    return [v, -w**2*u]

def run_solvers_and_plot(solvers, timesteps_per_period=20,
                         num_periods=1, I=1, w=2*np.pi):
    P = 2*np.pi/w  # duration of one period
    dt = P/timesteps_per_period
    Nt = num_periods*timesteps_per_period
    T = Nt*dt
    t_mesh = np.linspace(0, T, Nt+1)

    legends = []
    for solver in solvers:
        solver.set(f_kwargs={'w': w})
        solver.set_initial_condition([I, 0])
        u, t = solver.solve(t_mesh)

        # Compute energy
        dt = t[1] - t[0]
        E = 0.5*((u[2:,0] - u[:-2,0])/(2*dt))**2 + 0.5*w**2*u[1:-1,0]**2
        # Compute error in energy
        E0 = 0.5*0**2 + 0.5*w**2*I**2
        e_E = E - E0

        solver_name = 'CrankNicolson' if solver.__class__.__name__ == \
                      'MidpointImplicit' else solver.__class__.__name__
        print '*** Relative max error in energy for %s [0,%g] with dt=%g: %.3E' % (solver_name, t[-1], dt, np.abs(e_E).max()/E0)

        # Make plots
        if num_periods <= 80:
            plt.figure(1)
            if len(t_mesh) <= 50:
                plt.plot(t, u[:,0])             # markers by default
            else:
                plt.plot(t, u[:,0], '-2')       # no markers
            plt.hold('on')
            legends.append(solver_name)
            plt.figure(2)
            if len(t_mesh) <= 50:
                plt.plot(u[:,0], u[:,1])        # markers by default
            else:
                plt.plot(u[:,0], u[:,1], '-2')  # no markers
            plt.hold('on')

        if num_periods > 20:
            minima, maxima = minmax(t, u[:,0])
            p = periods(maxima)
            a = amplitudes(minima, maxima)
            plt.figure(3)
            plt.plot(range(len(p)), 2*np.pi/p, '-')
            plt.hold('on')
            plt.figure(4)
            plt.plot(range(len(a)), a, '-')
            plt.hold('on')

    # Compare with exact solution plotted on a very fine mesh
    t_fine = np.linspace(0, T, 10001)
    u_e = I*np.cos(w*t_fine)
    v_e = -w*I*np.sin(w*t_fine)

    if num_periods < 80:
        plt.figure(1)
        plt.plot(t_fine, u_e, '-') # avoid markers by spec. line type
        legends.append('exact')
        plt.legend(legends, loc='upper left')
        plt.xlabel('t');  plt.ylabel('u')
        plt.title('Time step: %g' % dt)
        plt.savefig('vib_%d_%d_u.png' % (timesteps_per_period, num_periods))
        plt.savefig('vib_%d_%d_u.pdf' % (timesteps_per_period, num_periods))

        plt.figure(2)
        plt.plot(u_e, v_e, '-') # avoid markers by spec. line type
        plt.legend(legends, loc='lower right')
        plt.xlabel('u(t)');  plt.ylabel('v(t)')
        plt.title('Time step: %g' % dt)
        plt.savefig('vib_%d_%d_pp.png' % (timesteps_per_period, num_periods))
        plt.savefig('vib_%d_%d_pp.pdf' % (timesteps_per_period, num_periods))
        del legends[-1]  # fig 3 and 4 does not have exact value

    if num_periods > 20:
        plt.figure(3)
        plt.legend(legends, loc='center right')
        plt.title('Empirically estimated periods')
        plt.savefig('vib_%d_%d_p.pdf' % (timesteps_per_period, num_periods))
        plt.savefig('vib_%d_%d_p.png' % (timesteps_per_period, num_periods))
        plt.figure(4)
        plt.legend(legends, loc='center right')
        plt.title('Empirically estimated amplitudes')
        plt.savefig('vib_%d_%d_a.pdf' % (timesteps_per_period, num_periods))
        plt.savefig('vib_%d_%d_a.png' % (timesteps_per_period, num_periods))

# Define different sets of experiments
solvers_theta = [
    odespy.ForwardEuler(f),
    # Implicit methods must use Newton solver to converge
    odespy.BackwardEuler(f, nonlinear_solver='Newton'),
    odespy.CrankNicolson(f, nonlinear_solver='Newton'),
    ]

solvers_RK = [odespy.RK2(f), odespy.RK4(f)]
solvers_accurate = [odespy.RK4(f),
                    odespy.CrankNicolson(f, nonlinear_solver='Newton'),
                    odespy.DormandPrince(f, atol=0.001, rtol=0.02)]
solvers_CN = [odespy.CrankNicolson(f, nonlinear_solver='Newton')]

if __name__ == '__main__':
    timesteps_per_period = 20
    solver_collection = 'theta'
    num_periods = 1
    try:
        # Example: python vib_odespy.py 30 accurate 50
        timesteps_per_period = int(sys.argv[1])
        solver_collection = sys.argv[2]
        num_periods = int(sys.argv[3])
    except IndexError:
        pass # default values are ok
    solvers = eval('solvers_' + solver_collection)  # list of solvers
    run_solvers_and_plot(solvers,
                         timesteps_per_period,
                         num_periods)
    #plt.show()
    #raw_input()
