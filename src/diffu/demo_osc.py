from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve, use_solver
from numpy import linspace, zeros
import time, os, sys, shutil

def solver(I, a, L, Nx, F, T, theta=0.5, u_L=0, u_R=0,
           user_action=None):
    """
    Solve the diffusion equation u_t = a*u_xx on (0,L) with
    boundary conditions u(0,t) = u_L and u(L,t) = u_R,
    for t in (0,T]. Initial condition: u(x,0) = I(x).

    Method: (implicit) theta-rule in time.

    Nx is the total number of mesh cells; mesh points are numbered
    from 0 to Nx.
    F is the dimensionless number a*dt/dx**2 and implicitly specifies the
    time step. No restriction on F.
    T is the stop time for the simulation.
    I is a function of x.

    user_action is a function of (u, x, t, n) where the calling code
    can add visualization, error computations, data analysis,
    store solutions, etc.

    The coefficient matrix is stored in a scipy data structure for
    sparse matrices. Input to the storage scheme is a set of
    diagonals with nonzero entries in the matrix.
    """
    import time
    t0 = time.clock()

    x = linspace(0, L, Nx+1)    # mesh points in space
    dx = x[1] - x[0]
    dt = F*dx**2/a
    Nt = int(round(T/float(dt)))
    print 'Number of time steps:', Nt
    t = linspace(0, T, Nt+1)    # mesh points in time

    u   = zeros(Nx+1)   # solution array at t[n+1]
    u_1 = zeros(Nx+1)   # solution at t[n]

    # Representation of sparse matrix and right-hand side
    diagonal = zeros(Nx+1)
    lower    = zeros(Nx+1)
    upper    = zeros(Nx+1)
    b        = zeros(Nx+1)

    # Precompute sparse matrix (scipy format)
    Fl = F*theta
    Fr = F*(1-theta)
    diagonal[:] = 1 + 2*Fl
    lower[:] = -Fl  #1
    upper[:] = -Fl  #1
    # Insert boundary conditions
    # (upper[1:] and lower[:-1] are the active alues)
    upper[0:2] = 0
    lower[-2:] = 0
    diagonal[0] = 1
    diagonal[Nx] = 1

    diags = [0, -1, 1]
    A = spdiags([diagonal, lower, upper], diags, Nx+1, Nx+1)
    #print A.todense()

    # Set initial condition
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Time loop
    for n in range(0, Nt):
        b[1:-1] = u_1[1:-1] + Fr*(u_1[:-2] - 2*u_1[1:-1] + u_1[2:])
        b[0] = u_L; b[-1] = u_R  # boundary conditions
        u[:] = spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, n+1)

        # Switch variables before next step
        u_1, u = u, u_1

    t1 = time.clock()
    return u, x, t, t1-t0


# Case: initial discontinuity

theta2name = {1: 'BE', 0: 'FE', 0.5: 'CN'}


class PlotU:
    def __init__(self, theta, F, Nx, L):
        self.theta, self.F, self.Nx, self.L  = theta, F, Nx, L
        self._make_plotdir()

    def __call__(self, u, x, t, n):
        from scitools.std import plot, savefig
        umin = -0.1; umax = 1.1  # axis limits for plotting
        title = 'Method: %s, F=%g, t=%f' % \
                (theta2name[self.theta], self.F, t[n])
        plot(x, u, 'r-',
             axis=[0, self.L, umin, umax],
             title=title)
        savefig(os.path.join(self.plotdir, 'frame_%04d.png' % n))

        # Pause the animation initially, otherwise 0.2 s between frames
        if n == 0:
            time.sleep(2)
        else:
            time.sleep(0.2)

    def _make_plotdir(self):
        self.plotdir = '%s_F%g' % (theta2name[self.theta], self.F)
        if os.path.isdir(self.plotdir):
            shutil.rmtree(self.plotdir)
        os.mkdir(self.plotdir)

    def make_movie(self):
        """Go to plot directory and make movie files."""
        orig_dir = os.getcwd()
        os.chdir(self.plotdir)
        cmd = 'avconv -r 1 -i frame_%04d.png -vcodec libvpx movie.webm'
        cmd = 'avconv -r 1 -i frame_%04d.png -vcodec flv movie.flv'
        cmd = 'avconv -r 1 -i frame_%04d.png -vcodec libtheora movie.ogg'
        os.system(cmd)
        cmd = 'scitools movie output_file=index.html frame*.png'
        os.system(cmd)
        os.chdir(orig_dir)

def I(x):
    return 0 if x > 0.5 else 1


def run_command_line_args():
    # Command-line arguments: Nx F theta
    Nx = int(sys.argv[1])
    F = float(sys.argv[2])
    theta = float(sys.argv[3])
    plot_u = PlotU(theta, F, Nx, L=1)
    u, x, t, cpu = solver(I, a=1, L=1, Nx=Nx, F=F, T=T,
                          theta=theta, u_L=1, u_R=0,
                          user_action=plot_u)
    plot_u.make_movie()

def run_BE_CN_FE():
    # cases: list of (Nx, C, theta, T) values
    cases = [(7, 5, 0.5, 3), (15, 0.5, 0, 0.25),
             (15, 0.5, 1, 0.12),]
    for Nx, F, theta, T in cases:
        print 'theta=%g, F=%g, Nx=%d' % (theta, F, Nx)
        plot_u = PlotU(theta, F, Nx, L=1)
        u, x, t, cpu = solver(I, a=1, L=1, Nx=Nx, F=F, T=T,
                              theta=theta, u_L=1, u_R=0,
                              user_action=plot_u)
        plot_u.make_movie()
        raw_input('Type Return to proceed with next case: ')

if __name__ == '__main__':
    if len(sys.argv) == 1:
        # No command-line arguments: run predefined cases
        run_BE_CN_FE()
    elif len(sys.argv) == 4:
        run_command_line_args()
