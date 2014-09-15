"""
Examples on approximating functions by global basis functions,
using the approx1D.py module.
"""
from approx1D import *
from Lagrange import *
import matplotlib.pyplot as plt
#import scitools.std as plt
import sympy as sp
import sys
x = sp.Symbol('x')


def sines(x, N):
    return [sp.sin(sp.pi*(i+1)*x) for i in range(N+1)]

def cosines(x, N):
    return [sp.cos(sp.pi*i*x) for i in range(N+1)]

def sines_cosines(x, N):
    c = [sp.cos(sp.pi*i*x) for i in range(N+1)]
    s = [sp.sin(sp.pi*i*x) for i in range(1, N+1)]
    return c + s

def taylor(x, N):
    return [x**i for i in range(N+1)]


# ----------------------------------------------------------------------

def run_linear_leastsq_parabola():
    f = 10*(x-1)**2 - 1
    psi = [1, x]
    Omega = [1, 2]
    u, c = least_squares(f, psi, Omega)
    comparison_plot(f, u, Omega, 'parabola_ls_linear')

def run_taylor_leastsq_parabola_illconditioning(N=2):
    """
    Test Taylor approx to a parabola and exact symbolics vs
    ill-conditioned numerical approaches.
    """
    f = 10*(x-1)**2 - 1
    u, c = least_squares(f, psi=[x**i for i in range(N+1)], Omega=[1, 2])
    # Note: in least_squares there is extra code for numerical solution
    # of the systems
    print 'f:', sp.expand(f)
    print 'u:', sp.expand(u)
    comparison_plot(f, u, [1, 2], 'parabola_ls_taylor%d' % N)

def run_sines_leastsq_parabola(boundary_term=False):
    for N in (4, 12):
        f = 10*(x-1)**2 - 1
        psi = sines(x, N)
        Omega = [0, 1]
        if boundary_term:
            f0 = 9; f1 = -1
            b = f0*(1-x) + x*f1
            u, c = least_squares_orth(f-b, psi, Omega)
            u = u + b
        else:
            u, c = least_squares_orth(f, psi, Omega)
        plt.figure()
        comparison_plot(f, u, Omega, 'parabola_ls_sines%d%s' %
                        (N, '_wfterm' if boundary_term else ''))


def run_sine_by_powers(N):
    f = sp.sin(x)
    psi = taylor(x, N)
    Omega=[0, 2*sp.pi]
    u = least_squares(f, psi, Omega)
    comparison_plot(f, u, Omega)

def run_Lagrange_poly(N):
    # Test of symbolic and numeric evaluation of Lagrange polynomials
    psi, points = Lagrange_polynomials_01(x, N)
    print psi
    print points
    x = 0.5
    psi, points = Lagrange_polynomials_01(x, N)
    print psi
    print points


def run_Lagrange_leastsq_sin(N, ymin=-1.2, ymax=1.2):
    # Least-squares use of Lagrange polynomials
    f = sp.sin(2*sp.pi*x)
    psi, points = Lagrange_polynomials_01(x, N)
    Omega=[0, 1]
    u = least_squares(f, psi, Omega)
    comparison_plot(f, u, Omega, filename='Lagrange_ls_sin_%d' % (N+1),
                    plot_title='Least squares approximation by '\
                    'Lagrange polynomials of degree %d' % N,
                    ymin=ymin, ymax=ymax)

def run_Lagrange_leastsq_abs(N):
    """Least-squares with of Lagrange polynomials for |1-2x|."""
    f = sp.abs(1-2*x)
    # This f will lead to failure of sympy integrate, fallback on numerical int.
    psi, points = Lagrange_polynomials_01(x, N)
    Omega=[0, 1]
    u = least_squares(f, psi, Omega)
    comparison_plot(f, u, Omega, filename='Lagrange_ls_abs_%d' % (N+1),
                    plot_title='Least squares approximation by '\
                    'Lagrange polynomials of degree %d' % N)

def run_linear_interp1_parabola():
    f = 10*(x-1)**2 - 1
    psi = [1, x]
    Omega = [1, 2]
    points = [1 + sp.Rational(1,3), 1 + sp.Rational(2,3)]
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega, 'parabola_interp1_linear')


def run_linear_interp2_parabola():
    # as run_linear_interp1_parabola, but other interpolation points
    f = 10*(x-1)**2 - 1
    psi = [1, x]
    Omega = [1, 2]
    points = [1, 2]
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega, 'parabola_interp2_linear')

def run_quadratic_interp_parabola():
    f = 10*(x-1)**2 - 1
    psi = [1, x, x**2]
    Omega = [1, 2]
    points = [1, 1.2, 2]
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega, 'parabola_interp3_quadratic')


def run_poly_interp_sin(N):
    f = sp.sin(sp.pi*x)
    psi = taylor(x, N)
    Omega = [1, 2]
    import numpy as np
    points = np.linspace(1, 2, N+1)
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega, 'sin_interp_poly%d' % N)


def run_Lagrange_interp_sin(N, ymin=-1.2, ymax=1.2):
    f = sp.sin(2*sp.pi*x)
    psi, points = Lagrange_polynomials_01(x, N)
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega=[0, 1],
                    filename='Lagrange_interp_sin_%d' % (N+1),
                    plot_title='Interpolation by Lagrange polynomials '\
                    'of degree %d' % N,
                    ymin=ymin, ymax=ymax)

def run_Lagrange_interp_poly(n, N):
    f = x**n
    psi, points = Lagrange_polynomials_01(x, N)
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega=[0, 1],
                    filename='Lagrange_interp_p%d_%d' % (n, N+1),
                    plot_title='Interpolation by Lagrange polynomials '\
                    'of degree %d' % N)

def run_Lagrange_interp_abs(N, ymin=None, ymax=None):
    f = abs(1-2*x)
    psi, points = Lagrange_polynomials_01(x, N)
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega=[0, 1],
                    filename='Lagrange_interp_abs_%d' % (N+1),
                    plot_title='Interpolation by Lagrange polynomials '\
                    'of degree %d' % N, ymin=ymin, ymax=ymax)
    # Make figures of Lagrange polynomials (psi)
    plt.figure()
    xcoor = np.linspace(0, 1, 1001)
    legends = []
    for i in (2, (N+1)/2+1):
        fn = sp.lambdify([x], psi[i])
        ycoor = fn(xcoor)
        plt.plot(xcoor, ycoor)
        legends.append(r'$\psi_%d$' % i)
        plt.hold('on')
    plt.legend(legends)
    plt.plot(points, [0]*len(points), 'ro')
    #if ymin is not None and ymax is not None:
    #    axis([xcoor[0], xcoor[-1], ymin, ymax])
    plt.savefig('Lagrange_basis_%d.pdf' % (N+1))
    plt.savefig('Lagrange_basis_%d.png' % (N+1))

def run_Lagrange_interp_abs_Cheb(N, ymin=None, ymax=None):
    f = sp.Abs(1-2*x)
    fn = sp.lambdify([x], f)
    psi, points= Lagrange_polynomials(x, N, [0, 1],
                                      point_distribution='Chebyshev')
    u = interpolation(f, psi, points)
    comparison_plot(f, u, Omega=[0, 1],
                    filename='Lagrange_interp_abs_Cheb_%d' % (N+1),
                    plot_title='Interpolation by Lagrange polynomials '\
                    'of degree %d' % N, ymin=ymin, ymax=ymax)
    print 'Interpolation points:', points

    # Make figures of Lagrange polynomials (psi)
    plt.figure()
    xcoor = np.linspace(0, 1, 1001)
    legends = []
    for i in (2, (N+1)/2+1):
        fn = sp.lambdify([x], psi[i])
        ycoor = fn(xcoor)
        plt.plot(xcoor, ycoor)
        legends.append(r'$\psi_%d$' % i)
        plt.hold('on')
    plt.legend(legends)
    plt.plot(points, [0]*len(points), 'ro')
    #if ymin is not None and ymax is not None:
    #    axis([xcoor[0], xcoor[-1], ymin, ymax])
    plt.savefig('Lagrange_basis_Cheb_%d.pdf' % (N+1))
    plt.savefig('Lagrange_basis_Cheb_%d.png' % (N+1))

def run_Lagrange_interp_abs_conv(N=[3, 6, 12, 24]):
    f = sp.abs(1-2*x)
    f = sp.sin(2*sp.pi*x)
    fn = sp.lambdify([x], f, modules='numpy')
    resolution = 50001
    import numpy as np
    xcoor = np.linspace(0, 1, resolution)
    fcoor = fn(xcoor)
    Einf = []
    E2 = []
    h = []
    for _N in N:
        psi, points = Lagrange_polynomials_01(x, _N)
        u = interpolation(f, psi, points)
        un = sp.lambdify([x], u, modules='numpy')
        ucoor = un(xcoor)
        e = fcoor - ucoor
        Einf.append(e.max())
        E2.append(np.sqrt(np.sum(e*e/e.size)))
        h.append(1./_N)
    print Einf
    print E2
    print h
    print N
    # Assumption: error = CN**(-N)
    print 'convergence rates:'
    for i in range(len(E2)):
        C1 = E2[i]/(N[i]**(-N[i]/2))
        C2 = Einf[i]/(N[i]**(-N[i]/2))
        print N[i], C1, C2
    # Does not work properly...


if __name__ == '__main__':
    functions = \
        [eval(fname) for fname in dir() if fname.startswith('run_')]
    from scitools.misc import function_UI
    cmd = function_UI(functions, sys.argv)
    eval(cmd)
