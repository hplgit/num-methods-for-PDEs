import sympy as sp

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

