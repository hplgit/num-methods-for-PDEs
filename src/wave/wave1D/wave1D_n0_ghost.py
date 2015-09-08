#!/usr/bin/env python
# As wave1D_dn0.py, but using ghost cells and index sets.
"""
1D wave equation with homogeneous Neumann conditions::

  u, x, t, cpu = solver(I, V, f, c, L, dt, C, T, user_action)

Function solver solves the wave equation

   u_tt = c**2*u_xx + f(x,t) on

(0,L) with du/dn=0 on x=0 and x = L.

dt is the time step.
T is the stop time for the simulation.
C is the Courant number (=c*dt/dx).
dx is computed from dt and C.

I and f are functions: I(x), f(x,t).
user_action is a function of (u, x, t, n) where the calling code
can add visualization, error computations, data analysis,
store solutions, etc.

Function viz::

  viz(I, V, f, c, L, dt, C, T, umin, umax, animate=True)

calls solver with a user_action function that can plot the
solution on the screen (as an animation).
"""
import numpy as np

def solver(I, V, f, c, L, dt, C, T, user_action=None):
    """
    Solve u_tt=c^2*u_xx + f on (0,L)x(0,T].
    u(0,t)=U_0(t) or du/dn=0 (U_0=None),
    u(L,t)=U_L(t) or du/dn=0 (u_L=None).
    """
    Nt = int(round(T/dt))
    t = np.linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = np.linspace(0, L, Nx+1)       # Mesh points in space
    C2 = C**2; dt2 = dt*dt            # Help variables in the scheme

    # Wrap user-given f, V
    if f is None or f == 0:
        f = (lambda x, t: 0)
    if V is None or V == 0:
        V = (lambda x: 0)

    u   = np.zeros(Nx+3)   # Solution array at new time level
    u_1 = np.zeros(Nx+3)   # Solution at 1 time level back
    u_2 = np.zeros(Nx+3)   # Solution at 2 time levels back

    Ix = range(1, u.shape[0]-1)
    It = range(0, t.shape[0])

    import time;  t0 = time.clock()  # for measuring CPU time

    # Load initial condition into u_1
    for i in Ix:
        u_1[i] = I(x[i-Ix[0]])  # Note the index transformation in x
    # Ghost values set according to du/dx=0
    i = Ix[0]
    u_1[i-1] = u_1[i+1]
    i = Ix[-1]
    u_1[i+1] = u_1[i-1]

    if user_action is not None:
        # Make sure to send the part of u that corresponds to x
        user_action(u_1[Ix[0]:Ix[-1]+1], x, t, 0)

    # Special formula for the first step
    for i in Ix:
        u[i] = u_1[i] + dt*V(x[i-Ix[0]]) + \
               0.5*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt2*f(x[i-Ix[0]], t[0])
    # Ghost values set according to du/dx=0
    i = Ix[0]
    u[i-1] = u[i+1]
    i = Ix[-1]
    u[i+1] = u[i-1]

    if user_action is not None:
        # Make sure to send the part of u that corresponds to x
        user_action(u[Ix[0]:Ix[-1]+1], x, t, 1)

    # Update data structures for next step
    #u_2[:] = u_1;  u_1[:] = u  # safe, but slower
    u_2, u_1, u = u_1, u, u_2

    for n in range(1, Nt):
        for i in Ix:
            u[i] = - u_2[i] + 2*u_1[i] + \
                   C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                   dt2*f(x[i-Ix[0]], t[n])
        # Ghost values set according to du/dx=0
        i = Ix[0]
        u[i-1] = u[i+1]
        i = Ix[-1]
        u[i+1] = u[i-1]

        if user_action is not None:
            # Make sure to send the part of u that corresponds to x
            if user_action(u[Ix[0]:Ix[-1]+1], x, t, n+1):
                break

        # Update data structures for next step
        #u_2[:] = u_1;  u_1[:] = u  # safe, but slower
        u_2, u_1, u = u_1, u, u_2

    # Important to correct the mathematically wrong u=u_2 above
    # before returning u
    u = u_1
    cpu_time = t0 - time.clock()
    return u[1:-1], x, t, cpu_time


from wave1D_u0 import viz
from wave1D_n0 import plug
# Cannot just import test_plug because wave1D_n0.test_plug will
# then call wave1D.solver, not the solver above

def test_plug():
    """
    Check that an initial plug is correct back after one period,
    if C=1.
    """
    L = 1.0
    I = lambda x: 0 if abs(x-L/2.0) > 0.1 else 1

    Nx = 10
    c = 0.5
    C = 1
    dt = C*(L/Nx)/c
    nperiods = 4
    T = L/c*nperiods  # One period: c*T = L
    u, x, t, cpu = solver(
        I=I, V=None, f=None, c=c, L=L,
        dt=dt, C=C, T=T, user_action=None)
    u_0 = np.array([I(x_) for x_ in x])
    diff = np.abs(u - u_0).max()
    tol = 1E-13
    assert diff < tol

if __name__ == '__main__':
    test_plug()
