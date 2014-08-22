# generate interactive demo session running the script
# below as input to scitools file2interactive
import sympy as sp
import sys

# Integration associated with sine expansion for -u''=2
i, j = sp.symbols('i j', integer=True)
x, L = sp.symbols('x L')
# Cannot do this one unless the arguments are i*x and j*x
#A_ij = sp.integrate(sp.sin((i+1)*sp.pi*x/L)*sp.sin((j+1)*sp.pi*x/L),
#                    (x, 0, L))
A_ii = sp.integrate(sp.sin((i+1)*sp.pi*x/L)**2, (x, 0, L))
print A_ii
f = 2
a = 2*L/(sp.pi**2*(i+1)**2)
c_i = a*sp.integrate(f*sp.sin((i+1)*sp.pi*x/L), (x, 0, L))
c_i = sp.simplify(c_i)
print c_i
print sp.latex(c_i, mode='plain')
#sys.exit(0)

x, x_m, h, X = sp.symbols('x x_m h X')

from fe_approx1D_numint import *
c = approximate(sp.sin(x), symbolic=True, d=1, N_e=4, numint='Trapezoidal',
                Omega=[0,sp.pi])
print c
c = approximate(sp.sin(x), symbolic=True, d=1, N_e=4, numint='Simpson',
                Omega=[0,sp.pi])
print c
#sys.exit(0)
from fe_approx1D import *

# "Hand"-integration of element matrix and vector
A_00 = sp.integrate(h/8*(1-X)**2, (X, -1, 1))
print A_00
print sp.latex(A_00, mode='plain')
A_10 = sp.integrate(h/8*(1+X)*(1-X), (X, -1, 1))
print A_10
print sp.latex(A_10, mode='plain')
A_11 = sp.integrate(h/8*(1+X)**2, (X, -1, 1))
print A_11
print sp.latex(A_11, mode='plain')
x = x_m + h/2*X
b_0 = sp.integrate(h/4*x*(1-x)*(1-X), (X, -1, 1))
b_1 = sp.integrate(h/4*x*(1-x)*(1+X), (X, -1, 1))
print b_0
print b_1
print sp.latex(b_0, mode='plain')
print sp.latex(b_1, mode='plain')

phi = basis(d=1)
phi
element_matrix(phi, Omega_e=[0.1, 0.2], symbolic=True)
element_matrix(phi, Omega_e=[0.1, 0.2], symbolic=False)

h, x = sp.symbols('h x')
nodes = [0, h, 2*h]
elements = [[0, 1], [1, 2]]
phi = basis(d=1)
f = x*(1-x)
A, b = assemble(nodes, elements, phi, f, symbolic=True)
A
b
c = A.LUsolve(b)
c
fn = sp.lambdify([x], f)
[fn(xc) for xc in nodes]

# The corresponding numerical computations, as done by sympy and
# still based on symbolic integration, goes as follows:

nodes = [0, 0.5, 1]
elements = [[0, 1], [1, 2]]
phi = basis(d=1)
x = sp.Symbol('x')
f = x*(1-x)
A, b = assemble(nodes, elements, phi, f, symbolic=False)
A
b
c = A.LUsolve(b)
c

d=1; N_e=8; Omega=[0,1]  # 8 linear elements on [0,1]
phi = basis(d)
f = x*(1-x)
nodes, elements = mesh_symbolic(N_e, d, Omega)
A, b = assemble(nodes, elements, phi, f, symbolic=True)
A

from fe_approx1D_numint import *
c = approximate(sp.sin(x), symbolic=True, d=1, N_e=4, numint='Trapezoidal',
                Omega=[0,sp.pi])
print c
c = approximate(sp.sin(x), symbolic=True, d=1, N_e=4, numint='Simpson',
                Omega=[0,sp.pi])
print c

# The integration does not work with sin(pi*x), but works fine with
# sin(x) on [0,pi] instead.
#approximate(sp.sin(sp.pi*x), symbolic=True, d=1, N_e=3, numint=None,
#            Omega=[0,1])
c = approximate(sp.sin(x), symbolic=True, d=1, N_e=2, numint=None,
                Omega=[0,sp.pi])
print sp.simplify(c[1,0].subs('h', sp.pi/2))

