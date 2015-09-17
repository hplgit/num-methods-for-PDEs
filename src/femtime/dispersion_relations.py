

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

methods = {'FE': A_FE, 'BE': A_BE, 'CN': A_CN}
C_values_coarse = {'2': 2, '1/2': 0.5}       # coarse mesh
C_values_fine = {'1/6': 1./6, '1/12': 1./12}  # fine mesh
from scitools.std import *

n = 16
p = linspace(0, pi/2, n)
for method in methods:
    for C_values in [C_values_coarse, C_values_fine]:
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
        xlabel('$p$'); ylabel(r'$A/A_{\rm e}$')
        filestem = 'diffu_A_factors_%s_%s' % \
                   ('coarse' if C_values == C_values_coarse else 'fine', method)
        savefig(filestem + '.png')
        savefig(filestem + '.pdf')
raw_input()
