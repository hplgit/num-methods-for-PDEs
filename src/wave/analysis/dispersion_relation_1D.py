
def tilde_w(c, k, dx, dt):
    C = dt*c/dx
    return 2/dt*asin(C*sin(k*dx/2))

def tilde_c(c, k, dx, dt):
    return tilde_w(c, k, dx, dt)/k

def r(C, p):
    return 1/(C*p)*asin(C*sin(p))  # important with 1, not 1. for sympy

def makeplot():
    n = 16
    p = linspace(0.001, pi/2, n)
    legends = []
    for C in 1.0, 0.95, 0.8, 0.3:
        plot(p, r(C, p))
        hold('on')
        legends.append('C=%g' % C)
    title('Numerical divided by exact wave velocity')
    legend(legends, fancybox=True, loc='lower left')
    axis([p[0], p[-1], 0.6, 1.1])
    xlabel('p')
    ylabel('velocity ratio')
    savefig('tmp.pdf')
    savefig('tmp.eps')
    savefig('tmp.png')
    show()

def sympy_analysis():
    C, p = symbols('C p')
    # Turn series expansion into polynomial and Python function
    # representations
    rs = r(C, p).series(p, 0, 7).removeO()
    print 'series representation of r(C, p):', rs
    print 'factored series representation of r(C, p):', factor(rs)
    rs1 = factor((rs - 1).extract_leading_order(p)[0][0])
    print 'leading order of the error:', rs1
    rs_poly = poly(rs)
    print 'polynomial representation:', rs_poly
    rs_pyfunc = lambdify([C, p], rs)  # can be used for plotting
    # Know that rs_pyfunc is correct (=1) when C=1, check that
    print rs_pyfunc(1, 0.1), rs_pyfunc(1, 0.76)

    # Alternative method for extracting terms in a series expansion:
    import itertools
    rs = [t for t in itertools.islice(r(C, p).lseries(p), 4)]
    print rs
    rs = [factor(t) for t in rs]
    print rs
    rs = sum(rs)
    print rs

    # true error
    x, t, k, w, c, dx, dt = symbols('x t k w c dx dt')
    u_n = cos(k*x - tilde_w(c, k, dx, dt)*t)
    u_e = cos(k*x - w*t)
    e = u_e-u_n
    # sympy cannot do this series expansion
    #print e.series(dx, 0, 4)

if __name__ == '__main__':
    #from scitools.std import *
    #makeplot()
    from sympy import * # erases sin and other math functions from numpy
    sympy_analysis()
