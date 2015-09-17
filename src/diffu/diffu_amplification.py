from numpy import *
#from matplotlib.pyplot import *
from scitools.easyviz.matplotlib_ import *

def A_exact(F, p):
    return exp(-4*F*p**2)

def A_FE(F, p):
    return 1 - 4*F*sin(p)**2

def A_BE(F, p):
    return 1/(1 + 4*F*sin(p)**2)

def A_CN(F, p):
    return (1 - 2*F*sin(p)**2)/(1 + 2*F*sin(p)**2)

def compare_plot(F, p):
    figure()
    plot(p, A_BE(F, p),
         p, A_exact(F, p),
         p, A_CN(F, p),
         p, A_FE(F, p),)
    legend(['BE', 'exact', 'CN', 'FE'])
    title('F=%g' % F)
    print 'F:', F
    if 0.2 >= F > 0.02:
        axis([p[0], p[-1], 0.3, 1])
    elif F <= 0.02:
        axis([p[0], p[-1], 0.75, 1])
    else:
        axis([p[0], p[-1], -1.2, 1])
    xlabel('$p=k\Delta x$')
    savefig('A_F%s.pdf' % (str(F).replace('.', '')))
    savefig('A_F%s.png' % (str(F).replace('.', '')))


p = linspace(0, pi/2, 101)
#for F in 20, 5, 2, 0.5, 0.25, 0.1, 0.01:
#    compare_plot(F, p)

#import sys; sys.exit(0)
from sympy import *
F, p, dx, dt, k = symbols('F p dx dt k')
A_err_FE = A_exact(F, p) - A_FE(F, p)
A_err_FE = A_FE(F, p)/A_exact(F, p)
#print 'Error in A, FE:', A_err_FE.series(F, 0, 6)
A_err_FE = A_err_FE.subs(F, dt/dx**2).subs(sin(p), 1).subs(p, k*dx/2)
print 'Error in A, FE:', A_err_FE.series(dt, 0, 3)
print latex(A_err_FE.series(F, 0, 6))
A_err_BE = A_exact(F, p) - A_BE(F, p)
A_err_BE = A_BE(F, p)/A_exact(F, p)
print 'Error in A, BE:', A_err_BE.series(F, 0, 6)
print latex(A_err_BE.series(F, 0, 6))
A_err_CN = A_exact(F, p) - A_CN(F, p)
A_err_CN = A_CN(F, p)/A_exact(F, p)
print 'Error in A, CN:', A_err_CN.series(F, 0, 6)
print latex(A_err_CN.series(F, 0, 6))

raw_input()

show()

f = open('tmp.sh', 'w')
f.write("""#!/bin/sh
doconce combine_images A_F20.pdf A_F2.pdf diffusion_A_F20_F2.pdf
doconce combine_images A_F20.png A_F2.png diffusion_A_F20_F2.png

doconce combine_images A_F05.png A_F025.png diffusion_A_F05_F025.png
doconce combine_images A_F05.pdf A_F025.pdf diffusion_A_F05_F025.pdf

doconce combine_images A_F01.pdf A_F001.pdf diffusion_A_F01_F001.pdf
doconce combine_images A_F01.png A_F001.png diffusion_A_F01_F001.png
""")
f.close()
import os
os.system('sh -x tmp.sh')
