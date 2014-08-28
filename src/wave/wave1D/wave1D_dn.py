#!/usr/bin/env python
"""
1D wave equation with Dirichlet or Neumann conditions::

  u, x, t, cpu = solver(I, V, f, c, U_0, U_L, L, dt, C, T,
                        user_action, version='scalar')

Function solver solves the wave equation

   u_tt = c**2*u_xx + f(x,t) on

(0,L) with u=U_0 or du/dn=0 on x=0, and u=u_L or du/dn=0
on x = L. If U_0 or U_L equals None, the du/dn=0 condition
is used, otherwise U_0(t) and/or U_L(t) are used for Dirichlet cond.
Initial conditions: u=I(x), u_t=V(x).

T is the stop time for the simulation.
dt is the desired time step.
C is the Courant number (=c*dt/dx).

I, f, U_0, U_L are functions: I(x), f(x,t), U_0(t), U_L(t).
U_0 and U_L can also be 0, or None, where None implies
du/dn=0 boundary condition. f and V can also be 0 or None
(equivalent to 0).

user_action is a function of (u, x, t, n) where the calling code
can add visualization, error computations, data analysis,
store solutions, etc.

Function viz::

  viz(I, V, f, c, U_0, U_L, L, dt, C, T, umin, umax,
      version='scalar', animate=True)

calls solver with a user_action function that can plot the
solution on the screen (as an animation).
"""
from scitools.std import *

def solver(I, V, f, c, U_0, U_L, L, dt, C, T,
           user_action=None, version='scalar'):
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

    # Wrap user-given f, I, V, U_0, U_L if None or 0
    if f is None or f == 0:
        f = (lambda x, t: 0) if version == 'scalar' else \
            lambda x, t: zeros(x.shape)
    if I is None or I == 0:
        I = (lambda x: 0) if version == 'scalar' else \
            lambda x: zeros(x.shape)
    if V is None or V == 0:
        V = (lambda x: 0) if version == 'scalar' else \
            lambda x: zeros(x.shape)
    if U_0 is not None:
        if isinstance(U_0, (float,int)) and U_0 == 0:
            U_0 = lambda t: 0
        # else: U_0(t) is a function
    if U_L is not None:
        if isinstance(U_L, (float,int)) and U_L == 0:
            U_L = lambda t: 0
        # else: U_L(t) is a function

    u   = zeros(Nx+1)   # Solution array at new time level
    u_1 = zeros(Nx+1)   # Solution at 1 time level back
    u_2 = zeros(Nx+1)   # Solution at 2 time levels back

    Ix = range(0, Nx+1)
    It = range(0, Nt+1)

    import time;  t0 = time.clock()  # CPU time measurement

    # Load initial condition into u_1
    for i in Ix:
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for the first step
    for i in Ix[1:-1]:
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt2*f(x[i], t[0])

    i = Ix[0]
    if U_0 is None:
        # Set boundary values du/dn = 0
        # x=0: i-1 -> i+1 since u[i-1]=u[i+1]
        # x=L: i+1 -> i-1 since u[i+1]=u[i-1])
        ip1 = i+1
        im1 = ip1  # i-1 -> i+1
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(u_1[im1] - 2*u_1[i] + u_1[ip1]) + \
               0.5*dt2*f(x[i], t[0])
    else:
        u[0] = U_0(dt)

    i = Ix[-1]
    if U_L is None:
        im1 = i-1
        ip1 = im1  # i+1 -> i-1
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(u_1[im1] - 2*u_1[i] + u_1[ip1]) + \
               0.5*dt2*f(x[i], t[0])
    else:
        u[i] = U_L(dt)

    if user_action is not None:
        user_action(u, x, t, 1)

    # Update data structures for next step
    u_2[:], u_1[:] = u_1, u

    for n in It[1:-1]:
        # Update all inner points
        if version == 'scalar':
            for i in Ix[1:-1]:
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                       dt2*f(x[i], t[n])

        elif version == 'vectorized':
            u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
                      C2*(u_1[0:-2] - 2*u_1[1:-1] + u_1[2:]) + \
                      dt2*f(x[1:-1], t[n])
        else:
            raise ValueError('version=%s' % version)

        # Insert boundary conditions
        i = Ix[0]
        if U_0 is None:
            # Set boundary values
            # x=0: i-1 -> i+1 since u[i-1]=u[i+1] when du/dn=0
            # x=L: i+1 -> i-1 since u[i+1]=u[i-1] when du/dn=0
            ip1 = i+1
            im1 = ip1
            u[i] = - u_2[i] + 2*u_1[i] + \
                   C2*(u_1[im1] - 2*u_1[i] + u_1[ip1]) + \
                   dt2*f(x[i], t[n])
        else:
            u[0] = U_0(t[n+1])

        i = Ix[-1]
        if U_L is None:
            im1 = i-1
            ip1 = im1
            u[i] = - u_2[i] + 2*u_1[i] + \
                   C2*(u_1[im1] - 2*u_1[i] + u_1[ip1]) + \
                   dt2*f(x[i], t[n])
        else:
            u[i] = U_L(t[n+1])

        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Update data structures for next step
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time


def viz(I, V, f, c, U_0, U_L, L, dt, C, T, umin, umax,
        version='scalar', animate=True):
    """Run solver and visualize u at each time level."""
    import scitools.std as plt, time, glob, os
    if callable(U_0):
        bc_left = 'u(0,t)=U_0(t)'
    elif U_0 is None:
        bc_left = 'du(0,t)/dx=0'
    else:
        bc_left = 'u(0,t)=0'
    if callable(U_L):
        bc_right = 'u(L,t)=U_L(t)'
    elif U_L is None:
        bc_right = 'du(L,t)/dx=0'
    else:
        bc_right = 'u(L,t)=0'

    def plot_u(u, x, t, n):
        """user_action function for solver."""
        plt.plot(x, u, 'r-',
                 xlabel='x', ylabel='u',
                 axis=[0, L, umin, umax],
                 title='t=%.3f, %s, %s' % (t[n], bc_left, bc_right))
        # Let the initial condition stay on the screen for 2
        # seconds, else insert a pause of 0.2 s between each plot
        time.sleep(2) if t[n] == 0 else time.sleep(0.2)
        plt.savefig('frame_%04d.png' % n)  # for movie making

    # Clean up old movie frames
    for filename in glob.glob('frame_*.png'):
        os.remove(filename)

    user_action = plot_u if animate else None
    u, x, t, cpu = solver(I, V, f, c, U_0, U_L, L, dt, C, T,
                          user_action, version)
    if animate:
        plt.movie('frame_*.png', encoder='html', fps=4,
                  output_file='movie.html')
        # Make other movie formats: Flash, Webm, Ogg, MP4
        codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                         libtheora='ogg')
        fps = 6
        filespec = 'frame_%04d.png'
        movie_program = 'avconv'  # or 'ffmpeg'
        for codec in codec2ext:
            ext = codec2ext[codec]
            cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
                  '-vcodec %(codec)s movie.%(ext)s' % vars()
            print cmd
            os.system(cmd)
    return cpu

import nose.tools as nt

def test_quadratic():
    """
    Check the scalar and vectorized versions work for
    a quadratic u(x,t)=x(L-x)(1+t/2) that is exactly reproduced.
    We simulate in [0, L/2] and apply a symmetry condition
    at the end x=L/2.
    """
    u_exact = lambda x, t: x*(L-x)*(1+0.5*t)
    I = lambda x: u_exact(x, 0)
    V = lambda x: 0.5*u_exact(x, 0)
    f = lambda x, t: 2*(1+0.5*t)*c**2
    U_0 = lambda t: u_exact(0, t)
    U_L = None
    L = 2.5
    c = 1.5
    C = 0.75
    Nx = 3  # Very coarse mesh for this exact test
    dt = C*(L/Nx)/c
    T = 18  # long time integration

    def assert_no_error(u, x, t, n):
        u_e = u_exact(x, t[n])
        diff = abs(u - u_e).max()
        nt.assert_almost_equal(diff, 0, places=13)

    solver(I, V, f, c, U_0, U_L, L/2, dt, C, T,
           user_action=assert_no_error, version='scalar')
    solver(I, V, f, c, U_0, U_L, L/2, dt, C, T,
           user_action=assert_no_error, version='vectorized')


def plug(C=1, Nx=50, animate=True, version='scalar', T=2, loc=0.5,
         bc_left='u=0', ic='u'):
    """Plug profile as initial condition."""
    L = 1.
    c = 1

    def I(x):
        if abs(x-loc) > 0.1:
            return 0
        else:
            return 1

    u_L = 0 if bc_left == 'u=0' else None
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    if ic == 'u':
        # u(x,0)=plug, u_t(x,0)=0
        cpu = viz(lambda x: 0 if abs(x-loc) > 0.1 else 1,
                  None, None, c, u_L, None, L, dt, C, T,
                  umin=-1.1, umax=1.1, version=version, animate=animate)
    else:
        # u(x,0)=0, u_t(x,0)=plug
        cpu = viz(None, lambda x: 0 if abs(x-loc) > 0.1 else 1,
                  None, c, u_L, None, L, dt, C, T,
                  umin=-0.25, umax=0.25, version=version, animate=animate)

def gaussian(C=1, Nx=50, animate=True, version='scalar', T=1, loc=5,
             bc_left='u=0', ic='u'):
    """Gaussian function as initial condition."""
    L = 10.
    c = 10
    sigma = 0.5

    def G(x):
        return 1/sqrt(2*pi*sigma)*exp(-0.5*((x-loc)/sigma)**2)

    u_L = 0 if bc_left == 'u=0' else None
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    umax = 1.1*G(loc)
    if ic == 'u':
        # u(x,0)=Gaussian, u_t(x,0)=0
        cpu = viz(G, None, None, c, u_L, None, L, dt, C, T,
                  umin=-umax, umax=umax, version=version, animate=animate)
    else:
        # u(x,0)=0, u_t(x,0)=Gaussian
        cpu = viz(None, G, None, c, u_L, None, L, dt, C, T,
                  umin=-umax/6, umax=umax/6, version=version, animate=animate)

def test_plug():
    """Check that an initial plug is correct back after one period."""
    L = 1
    c = 0.5
    dt = (L/10)/c  # Nx=10
    I = lambda x: 0 if abs(x-L/2.0) > 0.1 else 1

    u_s, x, t, cpu = solver(
        I=I,
        V=None, f=None, c=0.5, U_0=None, U_L=None, L=L,
        dt=dt, C=1, T=4, user_action=None, version='scalar')
    u_v, x, t, cpu = solver(
        I=I,
        V=None, f=None, c=0.5, U_0=None, U_L=None, L=L,
        dt=dt, C=1, T=4, user_action=None, version='vectorized')
    diff = abs(u_s - u_v).max()
    nt.assert_almost_equal(diff, 0, places=13)
    u_0 = array([I(x_) for x_ in x])
    diff = abs(u_s - u_0).max()
    nt.assert_almost_equal(diff, 0, places=13)

def guitar(C=1, Nx=50, animate=True, version='scalar', T=2):
    """Triangular initial condition for simulating a guitar string."""
    L = 1.
    c = 1
    x0 = 0.8*L
    dt = L/Nx/c  # choose the stability limit (if C<1, dx gets larger)
    I = lambda x: x/x0 if x < x0 else 1./(1-x0)*(1-x)

    cpu = viz(I, None, None, c, U_0, U_L, L, dt, C, T,
              umin=-1.1, umax=1.1, version=version, animate=True)
    print 'CPU time: %s version =' % version, cpu


def moving_end(C=1, Nx=50, reflecting_right_boundary=True, T=2,
               version='vectorized'):
    """
    Sinusoidal variation of u at the left end.
    Right boundary can be reflecting or have u=0, according to
    reflecting_right_boundary.
    """
    L = 1.
    c = 1
    dt = L/Nx/c  # choose the stability limit (if C<1, dx gets larger)
    I = lambda x: 0

    def U_0(t):
        return (0.25*sin(6*pi*t) if ((t < 1./6) or (0.5 + 3./12 <= t <= 0.5 + 4./12 + 0.0001) or (1.5 <= t <= 1.5 + 1./3 + 0.0001)) else 0)

    if reflecting_right_boundary:
        U_L = None
    else:
        U_L = 0
    umax = 1.1*0.5
    cpu = viz(I, None, None, c, U_0, U_L, L, dt, C, T,
              umin=-umax, umax=umax, version=version, animate=True)
    print 'CPU time: %s version =' % version, cpu


def sincos(C=1):
    """Test of exact analytical solution (sine in space, cosine in time)."""
    L = 10.0
    c = 1
    T = 5
    Nx = 80
    dt = (L/Nx)/c  # choose the stability limit with given Nx

    def u_exact(x, t):
        m = 3.0
        return cos(m*pi/L*t)*sin(m*pi/(2*L)*x)

    I = lambda x: u_exact(x, 0)
    U_0 = lambda t: u_exact(0, t)
    U_L = None # Neumann condition

    cpu = viz(I, None, None, c, U_0, U_L, L, dt, C, T,
              umin=-1.1, umax=1.1, version='scalar', animate=True)

    # Convergence study
    def action(u, x, t, n):
        e = abs(u - exact(x, t[n])).max()
        errors_in_time.append(e)

    E = []
    dt = []
    Nx_values = [10, 20, 40, 80, 160]
    for Nx in Nx_values:
        errors_in_time = []
        dt = (L/Nx)/c
        solver(I, None, None, c, U_0, U_L, L, dt, C, T,
               user_action=action, version='scalar')
        E.append(max(errors_in_time))
        _dx = L/Nx
        _dt = C*_dx/c
        dt.append(_dt)
        print dt[-1], E[-1]
    return dt, E

if __name__ == '__main__':
    import sys
    from scitools.misc import function_UI
    cmd = function_UI([test_quadratic, test_plug, plug,
                       gaussian, sincos, guitar,
                       moving_end, sincos], sys.argv)
    eval(cmd)
