import sympy as sp

# Computing with Dirichlet conditions: -u''=2 and sines
x, L = sp.symbols('x L')
e_Galerkin = x*(L-x) - 8*L**2*sp.pi**(-3)*sp.sin(sp.pi*x/L)
e_colloc = x*(L-x) - 2*L**2*sp.pi**(-2)*sp.sin(sp.pi*x/L)

# Verify max error for x=L/2
dedx_Galerkin = sp.diff(e_Galerkin, x)
print dedx_Galerkin.subs(x, L/2)
dedx_colloc = sp.diff(e_colloc, x)
print dedx_colloc.subs(x, L/2)

# Compute max error: x=L/2, evaluate numerical, and simplify
print 'Max error Galerkin/least.sq.:', \
      sp.simplify(e_Galerkin.subs(x, L/2).evalf(n=3))
print 'Max error colloc.:', \
      sp.simplify(e_colloc.subs(x, L/2).evalf(n=3))
import sys
sys.exit(0)


# Computing with Neumann and Dirichlet conditions: -u''=2
x, C, D = sp.symbols('x C D')
i, j = sp.symbols('i j', integer=True)

integrand = (i+1)*(j+1)*(1-x)**(i+j)
A_ij = sp.integrate(integrand, (x, 0, 1))
A_ij = sp.simplify(A_ij)
print A_ij
psi_i = (1-x)**(i+1)
integrand = 2*psi_i - D*(i+1)*(1-x)**i
b_i = sp.integrate(integrand, (x, 0, 1)) - C*psi_i.subs(x, 0)
b_i = sp.factor(sp.simplify(b_i))
print b_i
print sp.expand(2 - (2+i)*(D+C))
