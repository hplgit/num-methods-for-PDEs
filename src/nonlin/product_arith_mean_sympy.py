from sympy import *
t, dt = symbols('t dt')
P, Q = symbols('P Q', cls=Function)

# Target expression P(t_{n+1/2})*Q(t_{n+1/2})
# Simpler: P(0)*Q(0)
# Arithmetic means of each factor:
# 1/4*(P(-dt/2) + P(dt/2))*(Q(-dt/2) + Q(dt/2))
# Arithmetic mean of the product:
# 1/2*(P(-dt/2)*Q(-dt/2) + P(dt/2)*Q(dt/2))
# Let's Taylor expand to compare

target = P(0)*Q(0)
num_terms = 6
P_p = P(t).series(t, 0, num_terms).subs(t, dt/2)
print P_p
P_m = P(t).series(t, 0, num_terms).subs(t, -dt/2)
print P_m
Q_p = Q(t).series(t, 0, num_terms).subs(t, dt/2)
print Q_p
Q_m = Q(t).series(t, 0, num_terms).subs(t, -dt/2)
print Q_m

product_mean = Rational(1,2)*(P_m*Q_m + P_p*Q_p)
product_mean = simplify(expand(product_mean))
product_mean_error = product_mean - target

factor_mean = Rational(1,2)*(P_m + P_p)*Rational(1,2)*(Q_m + Q_p)
factor_mean = simplify(expand(factor_mean))
factor_mean_error = factor_mean - target

print 'product_mean_error:', product_mean_error
print 'factor_mean_error:', factor_mean_error
