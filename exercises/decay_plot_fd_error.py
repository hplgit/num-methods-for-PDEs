import sympy as sp

def D_f(u, dt):
    return (u(t + dt) - u(t))/dt

def D_b(u, dt):
    return (u(t) - u(t - dt))/dt

def D_c(u, dt):
    return (u(t + dt) - u(t - dt))/(2*dt)

a, t, dt, p = sp.symbols('a t dt p')

def u(t):
    return sp.exp(-a*t)

dudt = sp.diff(u(t), t)

from numpy import logspace, exp
from matplotlib.pyplot import plot, semilogx, legend, show, loglog

# Map operator function name to logical names
operator2name = dict(D_f='forward', D_b='backward', D_c='central')
legends = []
for operator in D_f, D_b, D_c:
    E = operator(u, dt)/dudt
    # Expand, set p=a*dt, simplify
    E = sp.expand(E)
    E = E.subs(a*dt, p)
    E = sp.simplify(E)
    print '%s E:' % operator2name[operator.__name__], E
    print 'Taylor series:', E.series(p, 0, 3)
    latex_expr = sp.latex(E)

    E = sp.lambdify([p], E, modules='numpy')
    p_values = logspace(-6, -0.5, 101)
    y = E(p_values)
    semilogx(p_values, y)
    legends.append(operator2name[operator.__name__] +
                   ': $' + latex_expr + '$')
legend(legends, loc='lower left')
show()
