from numpy import zeros, linspace

def solver_v1(I, w, dt, T):
    """
    Solve u'=v, v' = - w**2*u for t in (0,T], u(0)=I and v(0)=0,
    by a central finite difference method with time step dt.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = zeros(Nt+1)
    v = zeros(Nt+1)
    t = linspace(0, Nt*dt, Nt+1)  # mesh for u
    t_v = t + dt/2                # mesh for v

    u[0] = I
    v[0] = 0 - 0.5*dt*w**2*u[0]
    for n in range(1, Nt+1):
        u[n] = u[n-1] + dt*v[n-1]
        v[n] = v[n-1] - dt*w**2*u[n]
    return u, t, v, t_v

class HalfInt:
    """
    Class for allowing to write n+half and mean n,
    while n-half is n-1. Used for nice notation in staggered
    meshes.
    """
    def __radd__(self, other):
        return other

    def __rsub__(self, other):
        return other - 1

half = HalfInt()  # singleton object

def solver(I, w, dt, T):
    """
    Solve u'=v, v' = - w**2*u for t in (0,T], u(0)=I and v(0)=0,
    by a central finite difference method with time step dt.
    """
    dt = float(dt)
    Nt = int(round(T/dt))
    u = zeros(Nt+1)
    v = zeros(Nt+1)
    t = linspace(0, Nt*dt, Nt+1)  # mesh for u
    t_v = t + dt/2                # mesh for v

    u[0] = I
    v[0+half] = 0 - 0.5*dt*w**2*u[0]
    for n in range(1, Nt+1):
        print n, n+half, n-half
        u[n] = u[n-1] + dt*v[n-half]
        v[n+half] = v[n-half] - dt*w**2*u[n]
    return u, t, v, t_v

import nose.tools as nt

def test_staggered():
    I = 1.2; w = 2.0; T = 5; dt = 2/w
    u, t, v, t_v = solver(I, w, dt, T)
    from vib_undamped import solver as solver2
    u2, t2 = solver2(I, w, dt, T)
    error = abs(u - u2).max()
    nt.assert_almost_equal(error, 0, places=14)

if __name__ == '__main__':
    test_staggered()
