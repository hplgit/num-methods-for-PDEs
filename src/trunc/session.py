# To be executed by scitools file2interactive
from truncation_errors import TaylorSeries
from sympy import *
u, dt = symbols('u dt')
u_Taylor = TaylorSeries(u, 4)
u_Taylor(dt)
FE = (u_Taylor(dt) - u)/dt
FE
simplify(FE)

from truncation_errors import DiffOp
u = Symbol('u')
diffop = DiffOp(u, independent_variable='t')
diffop['geometric_mean']
diffop['Dtm']
diffop.operator_names()

def decay():
    u, a = symbols('u a')
    diffop = DiffOp(u, independent_variable='t',
                    num_terms_Taylor_series=3)
    D1u = diffop.D(1)   # symbol for du/dt
    ODE = D1u + a*u     # define ODE
    # Define schemes
    FE = diffop['Dtp'] + a*u
    CN = diffop['Dt' ] + a*u
    BE = diffop['Dtm'] + a*u
    theta = diffop['barDt'] + a*diffop['weighted_arithmetic_mean']
    theta = simplify(expand(theta))
    # Residuals (truncation errors)
    R = {'FE': FE-ODE, 'BE': BE-ODE, 'CN': CN-ODE, 'theta': theta-ODE}
    return R

print 'decay:', decay()

