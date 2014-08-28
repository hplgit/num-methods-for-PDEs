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
from scitools.std import *

def solver(I, V, f, c, L, dt, C, T, user_action=None):
    """
    Solve u_tt=c^2*u_xx + f on (0,L)x(0,T].
    u(0,t)=U_0(t) or du/dn=0 (U_0=None), u(L,t)=U_L(t) or du/dn=0 (u_L=None).
    """
    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = linspace(0, L, Nx+1)       # Mesh points in space
    C2 = C**2; dt2 = dt*dt         # Help variables in the scheme

    # Wrap user-given f, V
    if f is None or f == 0:
        f = (lambda x, t: 0)
    if V is None or V == 0:
        V = (lambda x: 0)

    u   = zeros(Nx+3)   # Solution array at new time level
    u_1 = zeros(Nx+3)   # Solution at 1 time level back
    u_2 = zeros(Nx+3)   # Solution at 2 time levels back

    Ix = range(1, u.shape[0]-2)
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
        u[i] = u_1[i] + dt*V(x[i]) + \
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
    u_2[:], u_1[:] = u_1, u

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
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u[1:-1], x, t, cpu_time


def viz(I, V, f, c, L, dt, C, T, umin, umax, animate=True):
    """Run solver and visualize u at each time level."""
    import scitools.std as plt, time, glob, os

    def plot_u(u, x, t, n):
        """user_action function for solver."""
        plt.plot(x, u, 'r-',
                 xlabel='x', ylabel='u',
                 axis=[0, L, umin, umax],
                 title='t=%f' % t[n])
        # Let the initial condition stay on the screen for 2
        # seconds, else insert a pause of 0.2 s between each plot
        time.sleep(2) if t[n] == 0 else time.sleep(0.2)
        plt.savefig('frame_%04d.png' % n)  # for movie making

    # Clean up old movie frames
    for filename in glob.glob('frame_*.png'):
        os.remove(filename)

    user_action = plot_u if animate else None
    u, x, t, cpu = solver(I, V, f, c, L, dt, C, T, user_action)
    if animate:
        plt.movie('frame_*.png', encoder='html', fps=4,
                  output_file='movie.html')
    return cpu

import nose.tools as nt


def plug(C=1, Nx=50, animate=True, T=2):
    """Plug profile as initial condition."""

    def I(x):
        if abs(x-L/2.0) > 0.1:
            return 0
        else:
            return 1

    cpu = viz(I, None, None, c, L, Nx, C, T,
              umin=-1.1, umax=1.1, animate=animate)


def test_plug():
    """Check that an initial plug is correct back after one period."""
    L = 1
    I = lambda x: 0 if abs(x-L/2.0) > 0.1 else 1

    u_s, x, t, cpu = solver(
        I=I, V=None, f=None, c=0.5, L=L,
        Nx=50, C=1, T=4, user_action=None)
    u_v, x, t, cpu = solver(
        I=I, V=None, f=None, c=0.5, L=L,
        Nx=50, C=1, T=4, user_action=None)
    diff = abs(u_s - u_v).max()
    nt.assert_almost_equal(diff, 0, places=13)
    u_0 = array([I(x_) for x_ in x])
    diff = abs(u_s - u_0).max()
    nt.assert_almost_equal(diff, 0, places=13)

if __name__ == '__main__':
    import sys
    from scitools.misc import function_UI
    cmd = function_UI([test_plug, plug], sys.argv)
    eval(cmd)
