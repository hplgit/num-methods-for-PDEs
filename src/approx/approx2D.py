"""
Approximation of functions by linear combination of basis functions in
function spaces and the least squares method or the collocation method
for determining the coefficients. 2D version.
"""
import sympy as sym
from scitools.std import surfc, savefig, linspace, title, figure, \
     hot, axis, ndgrid, subplot


def least_squares(f, psi, Omega, symbolic=True):
    """
    Given a function f(x,y) on a rectangular domain
    Omega=[[xmin,xmax],[ymin,ymax]],
    return the best approximation to f(x,y) in the space V
    spanned by the functions in the list psi.
    """
    N = len(psi) - 1
    A = sym.zeros((N+1, N+1))
    b = sym.zeros((N+1, 1))
    x, y = sym.symbols('x y')
    print '...evaluating matrix...'
    for i in range(N+1):
        for j in range(i, N+1):
            print '(%d,%d)' % (i, j)

            integrand = psi[i]*psi[j]
            if symbolic:
                I = sym.integrate(integrand,
                                 (x, Omega[0][0], Omega[0][1]),
                                 (y, Omega[1][0], Omega[1][1]))
            if not symbolic or isinstance(I, sym.Integral):
                # Could not integrate symbolically, use numerical int.
                print 'numerical integration of', integrand
                integrand = sym.lambdify([x,y], integrand)
                I = sym.mpmath.quad(integrand,
                                   [Omega[0][0], Omega[0][1]],
                                   [Omega[1][0], Omega[1][1]])
            A[i,j] = A[j,i] = I
        integrand = psi[i]*f
        if symbolic:
            I = sym.integrate(integrand,
                             (x, Omega[0][0], Omega[0][1]),
                             (y, Omega[1][0], Omega[1][1]))
        if not symbolic or isinstance(I, sym.Integral):
            # Could not integrate symbolically, use numerical int.
            print 'numerical integration of', integrand
            integrand = sym.lambdify([x,y], integrand)
            I = sym.mpmath.quad(integrand,
                               [Omega[0][0], Omega[0][1]],
                               [Omega[1][0], Omega[1][1]])
        b[i,0] = I
    print
    print 'A:\n', A, '\nb:\n', b
    if symbolic:
        c = A.LUsolve(b)  # symbolic solve
    else:
        c = sym.mpmath.lu_solve(A, b)  # numerical solve
    print 'coeff:', c

    # c is a sympy Matrix object, numbers are in c[i,0]
    u = sum(c[i,0]*psi[i] for i in range(len(psi)))
    print 'approximation:', u
    print 'f:', sym.expand(f)
    print sym.latex(A, mode='plain')
    print sym.latex(b, mode='plain')
    print sym.latex([c[i,0] for i in range(len(psi))], mode='plain')
    return u


def comparison_plot(f, u, Omega, plotfile='tmp'):
    """Compare f(x,y) and u(x,y) for x,y in Omega in a plot."""
    x, y = sym.symbols('x y')

    f = sym.lambdify([x,y], f, modules="numpy")
    u = sym.lambdify([x,y], u, modules="numpy")
    # When doing symbolics, Omega can easily contain symbolic expressions,
    # assume .evalf() will work in that case to obtain numerical
    # expressions, which then must be converted to float before calling
    # linspace below
    for r in range(2):
        for s in range(2):
            if not isinstance(Omega[r][s], (int,float)):
                Omega[r][s] = float(Omega[r][s].evalf())

    resolution = 41  # no of points in plot
    xcoor = linspace(Omega[0][0], Omega[0][1], resolution)
    ycoor = linspace(Omega[1][0], Omega[1][1], resolution)
    xv, yv = ndgrid(xcoor, ycoor)
    # Vectorized functions expressions does not work with
    # lambdify'ed functions without the modules="numpy"
    exact  = f(xv, yv)
    approx = u(xv, yv)
    figure()
    surfc(xv, yv, exact, title='f(x,y)',
          colorbar=True, colormap=hot(), shading='flat')
    if plotfile:
        savefig('%s_f.pdf' % plotfile, color=True)
        savefig('%s_f.png' % plotfile)
    figure()
    surfc(xv, yv, approx, title='f(x,y)',
          colorbar=True, colormap=hot(), shading='flat')
    if plotfile:
        savefig('%s_u.pdf' % plotfile, color=True)
        savefig('%s_u.png' % plotfile)

if __name__ == '__main__':
    print 'Module file not meant for execution.'
