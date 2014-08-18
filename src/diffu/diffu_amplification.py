from numpy import *
#from matplotlib.pyplot import *
from scitools.easyviz.matplotlib_ import *

def A_exact(Fo, p):
    return exp(-4*Fo*p**2)

def A_FE(Fo, p):
    return 1 - 4*Fo*sin(p)**2

def A_BE(Fo, p):
    return 1/(1 + 4*Fo*sin(p)**2)

def A_FoN(Fo, p):
    return (1 - 2*Fo*sin(p)**2)/(1 + 2*Fo*sin(p)**2)

def compare_plot(Fo, p):
    figure()
    plot(p, A_BE(Fo, p),
         p, A_exact(Fo, p),
         p, A_CN(Fo, p),
         p, A_FE(Fo, p),)
    legend(['BE', 'exact', 'CN', 'FE'])
    title('Fo=%g' % Fo)
    print 'Fo:', Fo
    if 0.2 >= Fo > 0.02:
        axis([p[0], p[-1], 0.3, 1])
    elif Fo <= 0.02:
        axis([p[0], p[-1], 0.75, 1])
    else:
        axis([p[0], p[-1], -1.2, 1])
    xlabel('$p=k\Delta x$')
    savefig('A_Fo%s.pdf' % (str(Fo).replace('.', '')))
    savefig('A_Fo%s.png' % (str(Fo).replace('.', '')))


p = linspace(0, pi/2, 101)
#for Fo in 20, 2, 0.5, 0.25, 0.1, 0.01:
#    compare_plot(Fo, p)

from sympy import *
Fo, p, dx, dt = symbols('Fo p dx dt')
#A_err_FE = A_FE(Fo, p)/A_exact(Fo, p)
A_err_FE = A_exact(Fo, p) - A_FE(Fo, p)
A_err_FE = A_FE(Fo, p)/A_exact(Fo, p)
#print 'Error in A, FE:', A_err_FE.series(Fo, 0, 6)
A_err_FE = A_err_FE.subs('Fo', 'dt/dx**2').subs('p', 'dx')
print 'Error in A, FE:', A_err_FE.series(dt, 0, 6)
print latex(A_err_FE.series(Fo, 0, 6))
A_err_BE = A_exact(Fo, p) - A_BE(Fo, p)
A_err_BE = A_BE(Fo, p)/A_exact(Fo, p)
print 'Error in A, BE:', A_err_BE.series(Fo, 0, 6)
print latex(A_err_BE.series(Fo, 0, 6))
A_err_CN = A_exact(Fo, p) - A_CN(Fo, p)
A_err_CN = A_CN(Fo, p)/A_exact(Fo, p)
print 'Error in A, CN:', A_err_CN.series(Fo, 0, 6)
print latex(A_err_CN.series(Fo, 0, 6))

raw_input()

show()

"""
doconce combine_images A_Fo20.pdf A_Fo2.pdf diffusion_A_Fo20_Fo2.pdf
doconce combine_images A_Fo20.png A_Fo2.png diffusion_A_Fo20_Fo2.png

doconce combine_images A_Fo05.png A_Fo025.png diffusion_A_Fo05_Fo025.png
doconce combine_images A_Fo05.pdf A_Fo025.pdf diffusion_A_Fo05_Fo025.pdf

doconce combine_images A_Fo01.pdf A_Fo001.pdf diffusion_A_Fo01_Fo001.pdf
doconce combine_images A_Fo01.png A_Fo001.png diffusion_A_Fo01_Fo001.png
"""
