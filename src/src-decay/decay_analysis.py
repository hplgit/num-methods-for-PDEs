"""
Theoretical estimates of errors in the theta rule via Taylor
series expansions and use of sympy for symbolic computations.
"""

# This file is aimed at being run as an interactive session via
# scitools file2interactive:
#
# Terminal> scitools file2interactive decay_analysis.py
#

import sys
from sympy import *

p = Symbol('p')
A_e = exp(-p)

# Demo on Taylor polynomials
A_e.series(p, 0, 6)

"""
# NOTE: rsolve can solve recurrence relations:
a, dt, I, n = symbols('a dt I n')
u = Function('u')
f = u(n+1) - u(n) + dt*a*u(n+1)
rsolve(f, u(n), {u(0): I})
# However, 0 is the answer!
# Experimentation shows that we cannot have symbols dt, a in the
# recurrence equation, just n or numbers.
# Even if we worked with scaled equations, dt is in there,
# rsolve cannot be applied.
"""

# Numerical amplification factor
theta = Symbol('theta')
A = (1-(1-theta)*p)/(1+theta*p)

# Interactive session for demonstrating subs
A.subs(theta, 1)                  # A for Backward Euler
A.subs(theta, Rational(1,2))      # Crank-Nicolson
A.subs(theta, 0).series(p, 0, 4)  # Taylor-expanded A for Forward Euler
A.subs(theta, 1).series(p, 0, 4)  # Taylor-expanded A for Backward Euler
A.subs(theta, Rational(1,2)).series(p, 0, 4) # Taylor-expanded A for C-N
A_e.series(p, 0, 4)               # Taylor-expanded exact A

# Error in amplification factors
half = Rational(1, 2)
FE = A_e.series(p, 0, 4) - A.subs(theta, 0).series(p, 0, 4)
BE = A_e.series(p, 0, 4) - A.subs(theta, 1).series(p, 0, 4)
CN = A_e.series(p, 0, 4) - A.subs(theta, half).series(p, 0, 4)
FE
BE
CN

# Ratio of amplification factors
FE = 1 - (A.subs(theta, 0)/A_e).series(p, 0, 4)
BE = 1 - (A.subs(theta, 1)/A_e).series(p, 0, 4)
CN = 1 - (A.subs(theta, half)/A_e).series(p, 0, 4)
FE
BE
CN

print "Error in solution:"
n, a, dt, t, T = symbols('n a dt t T')
u_e = exp(-p*n)
u_n = A**n
error = u_e.series(p, 0, 4) - u_n.subs(theta, 0).series(p, 0, 4)
print error
FE = error
error = error.subs('n', 't/dt').subs(p, 'a*dt')
#error = error.extract_leading_order(dt)[0][0]  # as_leading_term is simpler
error = error.as_leading_term(dt)
print 'Global error at a point t:', error
# error = sum(error.as_ordered_terms()[:-1]) # get rid of O() term, keep the rest
error_error_L2 = sqrt(integrate(error**2, (t, 0, T)))
print 'L2 error:', simplify(error_error_L2)
#error_error_L2 = error_error_L2.series(dt, 0, 3).as_leading_term(dt)  # series breaks down


sys.exit(0)
#BE = u_e.series(p, 0, 4) - u_n.subs(theta, 1).series(p, 0, 4)
#CN = u_e.series(p, 0, 4) - u_n.subs(theta, half).series(p, 0, 4)
FE
BE
CN
simplify(FE)
simplify(BE)
simplify(CN)





