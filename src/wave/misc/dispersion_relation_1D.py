
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
    rs = r(C, p).series(p, 0, 7)
    print rs

    for term in r(C, p).lseries(p): # get rid of the O() term
        print term                  # only two terms
    rs_factored = [factor(term) for term in rs.lseries(p)]
    rs_factored = sum(rs_factored)
    print rs_factored
    rs_f = poly(rs_factored)  # convert to polynomial form
    print rs_f
    print rs.extract_leading_order(p)

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

