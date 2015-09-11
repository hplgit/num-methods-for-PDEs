from sympy import *
x, w, b = symbols('x w b')
eq = x**2 + b*x + w**2
#s = solve(eq == 0, x)
#u = expand(exp(s[1]), complex=True)
#u = re(u)  # im(u)

# not smart enough for this:
"""
t = Symbol('t', real=True)
w_tilde = Symbol('w_tilde', real=True)
dt = Symbol('dt', real=True)

B = exp(I*w_tilde*t)
Q = (B.subs(t, t+dt) - 2*B + B.subs(t, t-dt))
print Q
Q = expand(Q/B, complex=True)
print Q
Q = simplify(Q)
print Q
"""

# Taylor series expansion of the numerical frequency
dt, w, t, T = symbols('dt w t T')
w_tilde_e = 2/dt*asin(w*dt/2)
w_tilde_series = w_tilde_e.series(dt, 0, 4)
print 'w_tilde series expansion:', w_tilde_series
print 'Error in frequency, leading order term:', \
      (w-w_tilde_series).as_leading_term(dt)
# Get rid of O() term
w_tilde_series = w_tilde_series.removeO()
print 'w_tilde series without O() term:', w_tilde_series

# The error mesh function (I=1)
#u_e = cos(w*t) - cos(w_tilde_e*t)  # problems with /dt around dt=0
error = cos(w*t) - cos(w_tilde_series*t)
print 'The global error:', error
print 'Series expansion of the global error:', error.series(dt, 0, 6)
print 'Series expansion of the global error:', error
error = error.series(dt, 0, 6).as_leading_term(dt)
print 'Leading order of the global error:', error
error_L2 = sqrt(integrate(error**2, (t, 0, T)))
print error_L2
#print error_L2.series(dt, 0, 2)  # break down
"""
error_L2 = simplify(error_L2.series(dt, 0, 4).as_leading_term(dt))
print 'L2 error:', error_L2
"""
