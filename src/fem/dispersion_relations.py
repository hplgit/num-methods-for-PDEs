

def A_FE(p, C, lumped=False):
    f = 1 if lumped else 1 - (2./3)*sin(p)**2
    return 1 - 4*C*sin(p)**2/f

def A_BE(p, C, lumped=False):
    f = 1 if lumped else 1 - (2./3)*sin(p)**2
    return 1/(1 + 4*C*sin(p)**2/f)  # 1 for sympy

def A_CN(p, C, lumped=False):
    f = 1 if lumped else 1 - (2./3)*sin(p)**2
    return (1 - 2*C*sin(p)**2/f)/(1 + 2*C*sin(p)**2/f)

def A_exact(p, C):
    return exp(-C*p**2)

methods = {'FE': A_FE, 'BE': A_BE}
C_values = {'2': 2, '1/2': 0.5}       # coarse mesh
#C_values = {'1/6': 1./6, '1/12': 1./12}  # fine mesh
from scitools.std import *
n = 16
p = linspace(0, pi/2, n)
for method in methods:
    figure()
    legends = []
    for C_name in C_values:
        C = C_values[C_name]
        A_e = A_exact(p, C)
        func = eval('A_' + method)
        A_FEM = func(p, C, False)
        A_FDM = func(p, C, True)
        legends.append('C=%s, FEM' % C_name)
        legends.append('C=%s, FDM' % C_name)
        plot(p, A_FEM, p, A_FDM)
        hold('on')
    plot(p, A_e)
    legends.append('exact')
    legend(legends, loc='lower left')
    if method == 'FE':
        plot([0, pi/2], [-1, -1], 'k--')
    title('Method: %s' % method)
    savefig('tmp%s_%s.png' % (len(C_values), method))
    savefig('tmp%s_%s.pdf' % (len(C_values), method))
raw_input()
