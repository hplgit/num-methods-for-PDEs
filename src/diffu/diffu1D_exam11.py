from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve, use_solver
from numpy import linspace, zeros
import time

def solver(I, a, L, Nx, Fo, T, theta=0.5, u_L=0, u_R=0,
           user_action=None):
    """
    Solve the diffusion equation u_t = a*u_xx on (0,L) with
    boundary conditions u(0,t) = u_L and u(L,t) = u_R,
    for t in (0,T]. Initial condition: u(x,0) = I(x).

    Method: (implicit) theta-rule in time.

    Nx is the total number of mesh cells; mesh points are numbered
    from 0 to Nx.
    Fo is the dimensionless number a*dt/dx**2 and implicitly specifies the
    time step. No restriction on Fo.
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

    x = linspace(0, L, Nx+1)   # mesh points in space
    dx = x[1] - x[0]
    dt = Fo*dx**2/a
    Nt = int(round(T/float(dt)))
    print 'Nt:', Nt
    t = linspace(0, T, Nt+1)    # mesh points in time

    u   = zeros(Nx+1)   # solution array at t[n+1]
    u_1 = zeros(Nx+1)   # solution at t[n]

    # Representation of sparse matrix and right-hand side
    diagonal = zeros(Nx+1)
    lower    = zeros(Nx+1)
    upper    = zeros(Nx+1)
    b        = zeros(Nx+1)

    # Precompute sparse matrix (scipy format)
    Fol = Fo*theta
    For = Fo*(1-theta)
    diagonal[:] = 1 + 2*Fol
    lower[:] = -Fol  #1
    upper[:] = -Fol  #1
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
        b[1:-1] = u_1[1:-1] + For*(u_1[:-2] - 2*u_1[1:-1] + u_1[2:])
        b[0] = u_L; b[-1] = u_R  # boundary conditions
        u[:] = spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, n+1)

        # Switch variables before next step
        u_1, u = u, u_1

    t1 = time.clock()
    return u, x, t, t1-t0


# Case: initial discontinuity


def plot_u(u, x, t, n):
    from scitools.std import plot
    umin = -0.1; umax = 1.1  # axis limits for plotting
    plot(x, u, 'r-', axis=[0, L, umin, umax], title='t=%f' % t[n])

    # Pause the animation initially, otherwise 0.2 s between frames
    if t[n] == 0:
        time.sleep(2)
    else:
        time.sleep(0.2)


L = 1
a = 1

def I(x):
    return 0 if x > L/2. else 1

# Command-line arguments: Nx Fo theta
import sys
Nx = 15
Fo = 0.5
theta = 0
T = 3
#theta = 1
#Nx = int(sys.argv[1])
#Fo = float(sys.argv[2])
#theta = float(sys.argv[3])

cases = [(7, 5, 0.5, 3), (15, 0.5, 0, 0.5),]
for Nx, Fo, theta, T in cases:
    print 'theta=%g, Fo=%g, Nx=%d' % (theta, Fo, Nx)
    u, x, t, cpu = solver(I, a, L, Nx, Fo, T,
                          theta=theta, u_L=1, u_R=0,
                          user_action=plot_u)
    raw_input('CR: ')
