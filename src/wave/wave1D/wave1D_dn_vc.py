#!/usr/bin/env python
"""
1D wave equation with Dirichlet or Neumann conditions
and variable wave velocity::

 u, x, t, cpu = solver(I, V, f, c, U_0, U_L, L, dt, C, T,
                       user_action=None, version='scalar',
                       stability_safety_factor=1.0)

Solve the wave equation u_tt = (c**2*u_x)_x + f(x,t) on (0,L) with
u=U_0 or du/dn=0 on x=0, and u=u_L or du/dn=0
on x = L. If U_0 or U_L equals None, the du/dn=0 condition
is used, otherwise U_0(t) and/or U_L(t) are used for Dirichlet cond.
Initial conditions: u=I(x), u_t=V(x).

T is the stop time for the simulation.
dt is the desired time step.
C is the Courant number (=max(c)*dt/dx).
stability_safety_factor enters the stability criterion:
C <= stability_safety_factor (<=1).

I, f, U_0, U_L, and c are functions: I(x), f(x,t), U_0(t),
U_L(t), c(x).
U_0 and U_L can also be 0, or None, where None implies
du/dn=0 boundary condition. f and V can also be 0 or None
(equivalent to 0). c can be a number or a function c(x).

user_action is a function of (u, x, t, n) where the calling code
can add visualization, error computations, data analysis,
store solutions, etc.
"""
import time, glob, shutil
from scitools.std import *

def solver(I, V, f, c, U_0, U_L, L, dt, C, T,
           user_action=None, version='scalar',
           stability_safety_factor=1.0):
    """Solve u_tt=(c^2*u_x)_x + f on (0,L)x(0,T]."""
    if isinstance(c, (float,int)):
        c = zeros(x.shape) + c
    elif callable(c):
        # Call c(x) and fill array c
        c_ = zeros(x.shape)
        for i in range(Nx+1):
            c_[i] = c(x[i])
        c = c_

    Nt = int(round(T/dt))
    t = linspace(0, Nt*dt, Nt+1)      # Mesh points in time
    dx = dt*c.max()/(stability_safety_factor*C)
    Nx = int(round(L/dx))
    x = linspace(0, L, Nx+1)          # Mesh points in space

    q = c**2
    C2 = (dt/dx)**2; dt2 = dt*dt      # Help variables in the scheme

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
    if U_L is not None:
        if isinstance(U_L, (float,int)) and U_L == 0:
            U_L = lambda t: 0

    u   = zeros(Nx+1)   # Solution array at new time level
    u_1 = zeros(Nx+1)   # Solution at 1 time level back
    u_2 = zeros(Nx+1)   # Solution at 2 time levels back

    import time;  t0 = time.clock()  # CPU time measurement

    Ix = range(0, Nx+1)
    It = range(0, Nt+1)

    # Load initial condition into u_1
    for i in range(0,Nx+1):
        u_1[i] = I(x[i])

    if user_action is not None:
        user_action(u_1, x, t, 0)

    # Special formula for the first step
    for i in Ix[1:-1]:
        u[i] = u_1[i] + dt*V(x[i]) + \
        0.5*C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i]) - \
                0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
        0.5*dt2*f(x[i], t[0])

    i = Ix[0]
    if U_0 is None:
        # Set boundary values (x=0: i-1 -> i+1 since u[i-1]=u[i+1]
        # when du/dn = 0, on x=L: i+1 -> i-1 since u[i+1]=u[i-1])
        ip1 = i+1
        im1 = ip1  # i-1 -> i+1
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
        0.5*dt2*f(x[i], t[0])
    else:
        u[i] = U_0(dt)

    i = Ix[-1]
    if U_L is None:
        im1 = i-1
        ip1 = im1  # i+1 -> i-1
        u[i] = u_1[i] + dt*V(x[i]) + \
               0.5*C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
        0.5*dt2*f(x[i], t[0])
    else:
        u[i] = U_L(dt)

    if user_action is not None:
        user_action(u, x, t, 1)

    # Update data structures for next step
    #u_2[:] = u_1;  u_1[:] = u  # safe, but slower
    u_2, u_1, u = u_1, u, u_2

    for n in It[1:-1]:
        # Update all inner points
        if version == 'scalar':
            for i in Ix[1:-1]:
                u[i] = - u_2[i] + 2*u_1[i] + \
                    C2*(0.5*(q[i] + q[i+1])*(u_1[i+1] - u_1[i])  - \
                        0.5*(q[i] + q[i-1])*(u_1[i] - u_1[i-1])) + \
                dt2*f(x[i], t[n])

        elif version == 'vectorized':
            u[1:-1] = - u_2[1:-1] + 2*u_1[1:-1] + \
            C2*(0.5*(q[1:-1] + q[2:])*(u_1[2:] - u_1[1:-1]) -
                0.5*(q[1:-1] + q[:-2])*(u_1[1:-1] - u_1[:-2])) + \
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
                   C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            dt2*f(x[i], t[n])
        else:
            u[i] = U_0(t[n+1])

        i = Ix[-1]
        if U_L is None:
            im1 = i-1
            ip1 = im1
            u[i] = - u_2[i] + 2*u_1[i] + \
                   C2*(0.5*(q[i] + q[ip1])*(u_1[ip1] - u_1[i])  - \
                       0.5*(q[i] + q[im1])*(u_1[i] - u_1[im1])) + \
            dt2*f(x[i], t[n])
        else:
            u[i] = U_L(t[n+1])

        if user_action is not None:
            if user_action(u, x, t, n+1):
                break

        # Update data structures for next step
        #u_2[:] = u_1;  u_1[:] = u  # safe, but slower
        u_2, u_1, u = u_1, u, u_2

    # Important to correct the mathematically wrong u=u_2 above
    # before returning u
    u = u_1
    cpu_time = t0 - time.clock()
    return u, x, t, cpu_time


import nose.tools as nt

def test_quadratic():
    """
    Check the scalar and vectorized versions work for
    a quadratic u(x,t)=x(L-x)(1+t/2) that is exactly reproduced,
    provided c(x) is constant.
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
        #print n, diff
        #print u
        #print u_e
        nt.assert_almost_equal(diff, 0, places=13)

    solver(I, V, f, c, U_0, U_L, L/2, dt, C, T,
           user_action=assert_no_error, version='scalar',
           stability_safety_factor=1)
    solver(I, V, f, c, U_0, U_L, L/2, dt, C, T,
           user_action=assert_no_error, version='vectorized',
           stability_safety_factor=1)

def test_plug():
    """Check that an initial plug is correct back after one period."""
    L = 1.
    c = 0.5
    dt = (L/10)/c  # Nx=10
    I = lambda x: 0 if abs(x-L/2.0) > 0.1 else 1

    u_s, x, t, cpu = solver(
        I=I,
        V=None, f=None, c=c, U_0=None, U_L=None, L=L,
        dt=dt, C=1, T=4, user_action=None, version='scalar')
    u_v, x, t, cpu = solver(
        I=I,
        V=None, f=None, c=c, U_0=None, U_L=None, L=L,
        dt=dt, C=1, T=4, user_action=None, version='vectorized')
    diff = abs(u_s - u_v).max()
    nt.assert_almost_equal(diff, 0, places=13)
    u_0 = array([I(x_) for x_ in x])
    diff = abs(u_s - u_0).max()
    nt.assert_almost_equal(diff, 0, places=13)


class PlotSolution:
    """
    Class for the user_action function in solver.
    Visualizes the solution only.
    """
    def __init__(self,
                 casename='tmp',    # Prefix in filenames
                 umin=-1, umax=1,   # Fixed range of y axis
                 pause_between_frames=None,  # Movie speed
                 backend='matplotlib',       # or 'gnuplot'
                 screen_movie=True, # Show movie on screen?
                 title='',          # Extra message in title
                 every_frame=1):    # Show every_frame frame
        self.casename = casename
        self.yaxis = [umin, umax]
        self.pause = pause_between_frames
        module = 'scitools.easyviz.' + backend + '_'
        exec('import %s as plt' % module)
        self.plt = plt
        self.screen_movie = screen_movie
        self.title = title
        self.every_frame = every_frame

        # Clean up old movie frames
        for filename in glob('frame_*.png'):
            os.remove(filename)

    def __call__(self, u, x, t, n):
        if n % self.every_frame != 0:
            return
        title = 't=%.3g' % t[n]
        if self.title:
            title = self.title + ' ' + title
        self.plt.plot(x, u, 'r-',
                     xlabel='x', ylabel='u',
                     axis=[x[0], x[-1],
                           self.yaxis[0], self.yaxis[1]],
                     title=title,
                     show=self.screen_movie)
        # pause
        if t[n] == 0:
            time.sleep(2)  # let initial condition stay 2 s
        else:
            if self.pause is None:
                pause = 0.2 if u.size < 100 else 0
            time.sleep(pause)

        self.plt.savefig('%s_frame_%04d.png' % (self.casename, n))

    def make_movie_file(self):
        """
        Create subdirectory based on casename, move all plot
        frame files to this directory, and generate
        an index.html for viewing the movie in a browser
        (as a sequence of PNG files).
        """
        # Make HTML movie in a subdirectory
        directory = self.casename
        if os.path.isdir(directory):
            shutil.rmtree(directory)   # rm -rf directory
        os.mkdir(directory)            # mkdir directory
        # mv frame_*.png directory
        for filename in glob('frame_*.png'):
            os.rename(filename, os.path.join(directory, filename))
        os.chdir(directory)        # cd directory
        self.plt.movie('frame_*.png', encoder='html',
                       output_file='index.html', fps=4)

        # Make other movie formats: Flash, Webm, Ogg, MP4
        codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                         libtheora='ogg')
        filespec = 'frame_%04d.png'
        movie_program = 'avconv'  # or 'ffmpeg'
        for codec in codec2ext:
            ext = codec2ext[codec]
            cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s '\
                  '-vcodec %(codec)s movie.%(ext)s' % vars()
            os.system(cmd)
        os.chdir(os.pardir)  # move back to parent directory

def demo_BC_plug(C=1, Nx=40, T=4):
    """Demonstrate u=0 and u_x=0 boundary conditions with a plug."""
    action = PlotSolution('plug', -1.3, 1.3, every_frame=1,
                          title='u(0,t)=0, du(L,t)/dn=0.')
    L = 1.
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    solver(I=lambda x: 0 if abs(x-L/2.0) > 0.1 else 1,
           V=0, f=0, c=1, U_0=lambda t: 0, U_L=None, L=L,
           dt=dt, C=C, T=T,
           user_action=action, version='vectorized',
           stability_safety_factor=1)
    action.make_movie_file()

def demo_BC_gaussian(C=1, Nx=80, T=4):
    """Demonstrate u=0 and u_x=0 boundary conditions with a bell function."""
    action = PlotSolution('gaussian', -1.3, 1.3, every_frame=1,
                          title='u(0,t)=0, du(L,t)/dn=0.')
    L = 1.
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    solver(I=lambda x: exp(-0.5*((x-0.5)/0.05)**2),
           V=0, f=0, c=1, U_0=lambda t: 0, U_L=None, L=L,
           dt=dt, C=C, T=T,
           user_action=action, version='vectorized',
           stability_safety_factor=1)
    action.make_movie_file()

def moving_end(C=1, Nx=50, reflecting_right_boundary=True,
               version='vectorized'):
    L = 1.
    c = 1
    dt = (L/Nx)/c  # choose the stability limit with given Nx
    T = 3
    I = lambda x: 0
    V = 0
    f = 0

    def U_0(t):
        retur 0.25*sin(6*pi*t) if t < 1./3 else 0

    if reflecting_right_boundary:
        U_L = None
        bc_right = 'du(L,t)/dx=0'
    else:
        U_L = 0
        bc_right = 'u(L,t)=0'

    action = PlotSolution('moving_end', -1, 1, every_frame=4,
             title='u(0,t)=0.25*sin(6*pi*t) if t < 1/3 else 0, '
                          + bc_right)
    solver(I, V, f, c, U_0, U_L, L, dt, C, T,
           user_action=action, version=version,
           stability_safety_factor=1)


class PlotMediumAndSolution(PlotSolution):
    def __init__(self, medium, **kwargs):
        """Mark medium in plot: medium=[x_L, x_R]."""
        self.medium = medium
        PlotSolution.__init__(self, **kwargs)

    def __call__(self, u, x, t, n):
        if n % self.every_frame != 0:
            return
        # Plot u and mark medium x=x_L and x=x_R
        x_L, x_R = self.medium
        umin, umax = self.yaxis
        title = 'Nx=%d, t=%f' % (x.size-1, t[n])
        if self.title:
            title = self.title + ' ' + title
        self.plt.plot(x, u, 'r-',
                      [x_L, x_L], [umin, umax], 'k--',
                      [x_R, x_R], [umin, umax], 'k--',
                      xlabel='x', ylabel='u',
                      axis=[x[0], x[-1], umin, umax],
                      title=title)
        if t[n] == 0:
            time.sleep(2)  # let initial condition stay 2 s
        # No sleep - this is used for large meshes
        self.plt.savefig('frame_%04d.png' % n)  # for movie making


def pulse(C=1, Nx=200, animate=True, version='vectorized', T=2,
          loc='center', pulse_tp='gaussian', slowness_factor=2,
          medium=[0.7, 0.9], every_frame=1, sigma=0.05):
    """
    Various peaked-shaped initial conditions on [0,1].
    Wave velocity is decreased by the slowness_factor inside
    medium. The loc parameter can be 'center' or 'left',
    depending on where the initial pulse is to be located.
    The sigma parameter governs the width of the pulse.
    """
    # Use scaled parameters: L=1 for domain length, c_0=1
    # for wave velocity outside the domain.
    L = 1.0
    c_0 = 1.0
    if loc == 'center':
        xc = L/2
    elif loc == 'left':
        xc = 0

    if pulse_tp in ('gaussian','Gaussian'):
        def I(x):
            return exp(-0.5*((x-xc)/sigma)**2)
    elif pulse_tp == 'plug':
        def I(x):
            return 0 if abs(x-xc) > sigma else 1
    elif pulse_tp == 'cosinehat':
        def I(x):
            # One period of a cosine
            w = 2
            a = w*sigma
            return 0.5*(1 + cos(pi*(x-xc)/a)) \
                   if xc - a <= x <= xc + a else 0

    elif pulse_tp == 'half-cosinehat':
        def I(x):
            # Half a period of a cosine
            w = 4
            a = w*sigma
            return cos(pi*(x-xc)/a) \
                   if xc - 0.5*a <= x <= xc + 0.5*a else 0
    else:
        raise ValueError('Wrong pulse_tp="%s"' % pulse_tp)

    def c(x):
        return c_0/slowness_factor \
               if medium[0] <= x <= medium[1] else c_0

    umin=-0.5; umax=1.5*I(xc)
    casename = '%s_Nx%s_sf%s' % \
               (pulse_tp, Nx, slowness_factor)
    action = PlotMediumAndSolution(
        medium, casename=casename, umin=umin, umax=umax,
        every_frame=every_frame, screen_movie=animate)

    dt = (L/Nx)/c  # choose the stability limit with given Nx
    # Lower C will then use this dt, but smaller Nx
    solver(I=I, V=None, f=None, c=c, U_0=None, U_L=None,
           L=L, dt=dt, C=C, T=T,
           user_action=action, version=version,
           stability_safety_factor=1)



if __name__ == '__main__':
    pass
