import scitools.std as plt
#import matplotlib.pyplot as plt
import sys
import odespy
import numpy as np

class RHS:
    """m*u'' + b*u' + k*u = A_F*cos(w_F*t)."""

    def __init__(self, m=1, b=0, k=2*np.pi, A_F=0.01, w_F=1.5,
                 I=1, V=0, damping='linear'):
        self.m, self.b, self.k = float(m), b, k
        self.A_F, self.w_F = A_F, w_F
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
        return [v, (-b_term - self.k*u +
                    self.A_F*np.cos(self.w_F*t))/self.m]

    def exact(self, t):
        # Valid for linear s(u)
        k, b, m, A_F, w_F, I, V = self.k, self.b, self.m, \
                                  self.A_F, self.w_F, self.I, self.V
        b_crit = 2*np.sqrt(k*m)
        w_e = np.sqrt(k/m)
        zeta = b/b_crit
        zeta1p = zeta + np.sqrt(zeta**2 - 1)
        zeta1m = zeta - np.sqrt(zeta**2 - 1)
        zeta2 = np.sqrt(zeta**2 - 1)

        if zeta > 1:
            # No oscillations
            sol1 = (V + w_e*zeta1p*I)/(2*w_e*zeta2)*np.exp(-w_e*zeta1m*t)
            sol2 = (V + w_e*zeta1m*I)/(2*w_e*zeta2)*np.exp(-w_e*zeta1p*t)
            u_h = sol1 - sol2
        elif zeta == 1:
            u_h = (I + (V + w_e*I)*t)*np.exp(-w_e*t)
        else:
            # Oscillations
            A = np.sqrt(I**2 + ((V + zeta*w_e*I)/(w_e*zeta2))**2)
            phi = np.arctan((V + zeta*w_e*I)/(I*w_e*zeta2))
            u_h = A*np.exp(-zeta*w_e*t)*np.cos(zeta2*w_e*t - phi)

        # Excitation: F=F_0*cos(w_F*t)
        # F_0 and w_F must be extracted from F......?
        phi_0 = np.arctan(b*w_F/(k - m*w_F**2))
        A_0 = A_F/np.sqrt((k - m*w_F**2)**2 + (b*w_F)**2)
        u_p = A_0*np.cos(omega*t - phi_0)

        # Test: all special cases...
        return u_h + u_p

# NOT FINISHED
def run_solvers_and_plot(solvers, timesteps_per_period=20,
                         num_periods=1, b=0):
    w = 2*np.pi    # frequency of undamped free oscillations
    P = 2*np.pi/w  # duration of one period
    dt = P/timesteps_per_period
    Nt = num_periods*timesteps_per_period
    T = Nt*dt
    t_mesh = np.linspace(0, T, Nt+1)
    t_fine = np.linspace(0, T, 8*Nt+1)  # used for very accurate solution

    legends = []
    solver_exact = odespy.RK4(f)

    for solver in solvers:
        solver.set_initial_condition([solver.users_f.I, 0])
        u, t = solver.solve(t_mesh)

        solver_name = 'CrankNicolson' if solver.__class__.__name__ == \
                      'MidpointImplicit' else solver.__class__.__name__

        # Make plots (plot last 10 periods????)
        if num_periods <= 80:
            plt.figure(1)
            if len(t_mesh) <= 50:
                plt.plot(t, u[:,0])             # markers by default
            else:
                plt.plot(t, u[:,0], '-2')       # no markers
            plt.hold('on')
            legends.append(solver_name)

    # Compare with exact solution plotted on a very fine mesh
    #t_fine = np.linspace(0, T, 10001)
    #u_e = solver.users_f.exact(t_fine)
    # Compare with RK4 on a much finer mesh
    solver_exact.set_initial_condition([solver.users_f.I, 0])
    u_e, t_e = solver_exact.solve(t_fine)

    if num_periods < 80:
        plt.figure(1)
        plt.plot(t_e, u_e[:,0], '-') # avoid markers by spec. line type
        legends.append('exact (RK4)')
        plt.legend(legends, loc='upper left')
        plt.xlabel('t');  plt.ylabel('u')
        plt.title('Time step: %g' % dt)
        plt.savefig('vib_%d_%d_u.png' % (timesteps_per_period, num_periods))
        plt.savefig('vib_%d_%d_u.pdf' % (timesteps_per_period, num_periods))
        plt.savefig('vib_%d_%d_u.eps' % (timesteps_per_period, num_periods))

#f = RHS(b=0.4, A_F=1, w_F=2)
f = RHS(b=0.4, A_F=0, w_F=2)
f = RHS(b=0.4, A_F=1, w_F=np.pi)
f = RHS(b=0.4, A_F=20, w_F=2*np.pi) # qualitatively wrong FE, almost ok BE, smaller T
#f = RHS(b=0.4, A_F=20, w_F=0.5*np.pi)  # cool, FE almost there, BE good

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
    plt.show()
    raw_input()

