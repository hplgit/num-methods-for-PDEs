"""
Solve the diffusion equation

    u_t = (a(x)*u_x)_x + f(x,t)

on (0,L) with boundary conditions u(0,t) = u_L and u(L,t) = u_R,
for t in (0,T]. Initial condition: u(x,0) = I(x).

The following naming convention of variables are used.

===== ==========================================================
Name  Description
===== ==========================================================
Nx    The total number of mesh cells; mesh points are numbered
      from 0 to Nx.
T     The stop time for the simulation.
I     Initial condition (Python function of x).
a     Variable coefficient (constant).
L     Length of the domain ([0,L]).
x     Mesh points in space.
t     Mesh points in time.
n     Index counter in time.
u     Unknown at current/new time level.
u_1   u at the previous time level.
dx    Constant mesh spacing in x.
dt    Constant mesh spacing in t.
===== ==========================================================

``user_action`` is a function of ``(u, x, t, n)``, ``u[i]`` is the
solution at spatial mesh point ``x[i]`` at time ``t[n]``, where the
calling code can add visualization, error computations, data analysis,
store solutions, etc.
"""

from scipy.sparse import spdiags
from scipy.sparse.linalg import spsolve, use_solver
from numpy import linspace, zeros, random
import time, sys


def solver_theta(I, a, L, Nx, D, T, theta=0.5, u_L=1, u_R=0,
                 user_action=None):
    """
    The a variable is an array of length Nx+1 holding the values of
    a(x) at the mesh points.

    Method: (implicit) theta-rule in time.

    Nx is the total number of mesh cells; mesh points are numbered
    from 0 to Nx.
    D = dt/dx**2 and implicitly specifies the time step.
    T is the stop time for the simulation.
    I is a function of x.

    user_action is a function of (u, x, t, n) where the calling code
    can add visualization, error computations, data analysis,
    store solutions, etc.
    """
    import time
    t0 = time.clock()

    x = linspace(0, L, Nx+1)   # mesh points in space
    dx = x[1] - x[0]
    dt = D*dx**2
    #print 'dt=%g' % dt
    Nt = int(round(T/float(dt)))
    t = linspace(0, T, Nt+1)   # mesh points in time

    u   = zeros(Nx+1)   # solution array at t[n+1]
    u_1 = zeros(Nx+1)   # solution at t[n]

    """
    Basic formula in the scheme:

    0.5*(a[i+1] + a[i])*(u[i+1] - u[i]) -
    0.5*(a[i] + a[i-1])*(u[i] - u[i-1])

    0.5*(a[i+1] + a[i])*u[i+1]
    0.5*(a[i] + a[i-1])*u[i-1]
    -0.5*(a[i+1] + 2*a[i] + a[i-1])*u[i]
    """

    Dl = 0.5*D*theta
    Dr = 0.5*D*(1-theta)

    # Representation of sparse matrix and right-hand side
    diagonal = zeros(Nx+1)
    lower    = zeros(Nx+1)
    upper    = zeros(Nx+1)
    b        = zeros(Nx+1)
    # "Active" values: diagonal[:], upper[1:], lower[:-1]

    # Precompute sparse matrix (scipy format)
    diagonal[1:-1] = 1 + Dl*(a[2:] + 2*a[1:-1] + a[:-2])
    lower[0:-2] = -Dl*(a[1:-1] + a[:-2])
    upper[2:]   = -Dl*(a[2:] + a[1:-1])
    # Insert boundary conditions
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
        b[1:-1] = u_1[1:-1] + Dr*(
            (a[2:] + a[1:-1])*(u_1[:-2] - u_1[1:-1]) -
            (a[2:] + a[0:-2])*(u_1[1:-1] - u_1[:-2]))
        b[0] = u_L; b[-1] = u_R  # boundary conditions
        u[:] = spsolve(A, b)

        if user_action is not None:
            user_action(u, x, t, n+1)

        # Switch variables before next step
        u_1, u = u, u_1

    t1 = time.clock()
    return u, x, t, t1-t0


def viz(I, a, L, Nx, D, T, umin, umax, theta, u_L, u_R,
        animate=True):

    from scitools.std import plot
    def plot_u(u, x, t, n):
        plot(x, u, 'r-', axis=[0, L, umin, umax], title='t=%f' % t[n])
        if t[n] == 0:
            time.sleep(3)
        else:
            time.sleep(0.5)

    user_action = plot_u if animate else lambda u,x,t,n: None

    u, x, t, cpu = solver_theta(
        I, a, L, Nx, D, T, theta, u_L, u_R, user_action=user_action)
    return u, x, cpu

def fill_a(a_consts, L, Nx):
    """
    *a_consts*: ``[[x0, a0], [x1, a1], ...]`` is a
    piecewise constant function taking the value ``a0`` in ``[x0,x1]``,
    ``a1`` in ``[x1,x2]``, and so forth.

    Return a finite difference function ``a`` on a uniform mesh with
    Nx+1 points in [0, L] where the function takes on the piecewise
    constant values of *a_const*. That is,

    ``a[i] = a_consts[s][1]`` if ``x[i]`` is in subdomain
    ``[a_consts[s][0], a_consts[s+1][0]]``.
    """
    a = zeros(Nx+1)
    x = linspace(0, L, Nx+1)
    s = 0  # subdomain counter
    for i in range(len(x)):
        if s < len(a_consts)-1 and x[i] > a_consts[s+1][0]:
            s += 1
        a[i] = a_consts[s][1]
    return a

def u_exact_stationary(x, a, u_L, u_R):
    """
    Return stationary solution of a 1D variable coefficient
    Laplace equation: (a(x)*v'(x))'=0, v(0)=u_L, v(L)=u_R.

    v(x) = u_L + (u_R-u_L)*(int_0^x 1/a(c)dc / int_0^L 1/a(c)dc)
    """
    Nx = x.size - 1
    g = zeros(Nx+1)    # integral of 1/a from 0 to x
    dx = x[1] - x[0]   # assumed constant
    i = 0
    g[i] = 0.5*dx/a[i]
    for i in range(1, Nx):
        g[i] = g[i-1] + dx/a[i]
    i = Nx
    g[i] = g[i-1] + 0.5*dx/a[i]
    v = u_L + (u_R - u_L)*g/g[-1]
    return v



