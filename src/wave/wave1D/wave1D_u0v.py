#!/usr/bin/env python
"""
1D wave equation with u=0 at the boundary.
The solver function here offers scalar and vectorized versions.
See wave1D_u0_s.py for documentation. The only difference
is that function solver takes an additional argument "version":
version='scalar' implies explicit loops over mesh point,
while version='vectorized' provides a vectorized version.
"""
from numpy import *

def solver(I, V, f, c, L, dt, C, T, user_action=None,
           version='vectorized'):
    """Solve u_tt=c^2*u_xx + f on (0,L)x(0,T]."""
    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)   # Mesh points in time
    dx = dt*c/float(C)
    Nx = int(round(L/dx))
    x = linspace(0, L, Nx+1)       # Mesh points in space
    C2 = C**2                      # Help variable in the scheme
    if f is None or f == 0:
        f = (lambda x, t: 0) if version == 'scalar' else \
            lambda x, t: zeros(x.shape)
    if V is None or V == 0:
        V = (lambda x: 0) if version == 'scalar' else \
            lambda x: zeros(x.shape)

    u   = zeros(Nx+1)   # Solution array at new time level
    u_1 = zeros(Nx+1)   # Solution at 1 time level back
    u_2 = zeros(Nx+1)   # Solution at 2 time levels back

    import time;  t0 = time.clock()  # for measuring CPU time

    # Load initial condition into u_1
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for first time step
    n = 0
    for i in range(1, Nx):
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
               0.5*dt**2*f(x[i], t[n])
    u[0] = 0;  u[Nx] = 0

    if user_action is not None:
        user_action(u, x, t, 1)

    # Switch variables before next step
    u_2[:], u_1[:] = u_1, u

    for n in range(1, Nt):
        # Update all inner points at time t[n+1]

        if version == 'scalar':
            for i in range(1, Nx):
                u[i] = - u_2[i] + 2*u_1[i] + \
                       C2*(u_1[i-1] - 2*u_1[i] + u_1[i+1]) + \
                       dt**2*f(x[i], t[n])
        elif version == 'vectorized':   # (1:-1 slice style)
            f_a = f(x, t[n])  # Precompute in array
            u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
                C2*(u_1[0:-2] - 2*u_1[1:-1] + u_1[2:]) + \
                dt**2*f_a[1:-1]
        elif version == 'vectorized2':  # (1:Nx slice style)
            f_a = f(x, t[n])  # Precompute in array
            u[1:Nx] =  - u_2[1:Nx] + 2*u_1[1:Nx] + \
                C2*(u_1[0:Nx-1] - 2*u_1[1:Nx] + u_1[2:Nx+1]) + \
                dt**2*f_a[1:Nx]

        # Insert boundary conditions
        u[0] = 0;  u[Nx] = 0
        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Switch variables before next step
        u_2[:], u_1[:] = u_1, u

    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time

def viz(I, V, f, c, L, dt, C, T, umin, umax, animate=True,
        version='vectorized'):
    """Run solver and visualize u at each time level."""
    import scitools.std as plt, time, glob, os
    #num_frames = 100 # max no of frames in movie

    def plot_u(u, x, t, n):
        """user_action function for solver."""
        try:
            every = t.size/num_frames
        except NameError:
            every = 1  # plot every frame
        if n % every == 0:
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
    u, x, t, cpu = solver(I, V, f, c, L, dt, C, T,
                          user_action, version)
    if not animate:
        return cpu

    # Make movie files
    fps = 4  # Frames per second
    plt.movie('frame_*.png', encoder='html', fps=fps,
              output_file='movie.html')
    # Ex: avconv -r 4 -i frame_%04d.png -vcodec libtheora movie.ogg
    codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                     libtheora='ogg')
    filespec = 'frame_%04d.png'
    movie_program = 'avconv'  # or 'ffmpeg'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
              '-vcodec %(codec)s movie.%(ext)s' % vars()
        os.system(cmd)
    return cpu

import nose.tools as nt

def test_quadratic():
    """
    Check the scalar and vectorized versions work for
    a quadratic u(x,t)=x(L-x)(1+t/2) that is exactly reproduced.
    """
    # The following function must work for x as array or scalar
    u_exact = lambda x, t: x*(L - x)*(1 + 0.5*t)
    I = lambda x: u_exact(x, 0)
    V = lambda x: 0.5*u_exact(x, 0)
    # f is a scalar (zeros_like(x) works for scalar x too)
    f = lambda x, t: zeros_like(x) + 2*c**2*(1 + 0.5*t)

    L = 2.5
    c = 1.5
    C = 0.75
    Nx = 3  # Very coarse mesh for this exact test
    dt = C*(L/Nx)/c
    T = 18

    def assert_no_error(u, x, t, n):
        u_e = u_exact(x, t[n])
        diff = abs(u - u_e).max()
        nt.assert_almost_equal(diff, 0, places=13)

    solver(I, V, f, c, L, dt, C, T,
           user_action=assert_no_error, version='scalar')
    solver(I, V, f, c, L, dt, C, T,
           user_action=assert_no_error, version='vectorized')

def guitar(C):
    """Triangular wave (pulled guitar string)."""
    L = 0.75
    x0 = 0.8*L
    a = 0.005
    freq = 440
    wavelength = 2*L
    c = freq*wavelength
    omega = 2*pi*freq
    num_periods = 1
    T = 2*pi/omega*num_periods
    # Choose dt the same as the stability limit for Nx=50
    dt = L/50./c

    def I(x):
        return a*x/x0 if x < x0 else a/(L-x0)*(L-x)

    umin = -1.2*a;  umax = -umin
    cpu = viz(I, 0, 0, c, L, dt, C, T, umin, umax, animate=True)

def run_efficiency_experiments():
    L = 1
    x0 = 0.8*L
    a = 1
    c = 2
    T = 8
    C = 0.9
    umin = -1.2*a;  umax = -umin

    def I(x):
        return a*x/x0 if x < x0 else a/(L-x0)*(L-x)

    intervals = []
    speedup = []
    for Nx in [50, 100, 200, 400, 800]:
        dt = C*int(round(float(L)/Nx))/float(c)
        print 'solving scalar Nx=%d' % Nx,
        cpu_s = viz(I, 0, 0, c, L, dt, C, T, umin, umax,
                    animate=False, version='scalar')
        print cpu_s
        print 'solving vectorized Nx=%d' % Nx,
        cpu_v = viz(I, 0, 0, c, L, dt, C, T, umin, umax,
                    animate=False, version='vectorized')
        print cpu_v
        intervals.append(Nx)
        speedup.append(cpu_s/float(cpu_v))
        print 'Nx=%3d: cpu_v/cpu_s: %.3f' % (Nx, 1./speedup[-1])
    print 'Nx:', intervals
    print 'Speed-up:', speedup

if __name__ == '__main__':
    test_quadratic()  # verify
    import sys
    try:
        C = float(sys.argv[1])
        print 'C=%g' % C
    except IndexError:
        C = 0.85
    guitar(C)
    #run_efficiency_experiments()
